PROJECT (atus-pro)

cmake_minimum_required(VERSION 2.8)

find_package(MPI REQUIRED)
find_package(Doxygen)

set(CMAKE_CXX_FLAGS_RELEASE "-std=c++0x -march=native -O3 -mieee-fp -funroll-loops -ftree-vectorize -w -s")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++0x -O2 -g -mieee-fp -Wall")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug" )
endif(NOT CMAKE_BUILD_TYPE)

option( BUILD_3D "3D" OFF )
option( BUILD_HTRAP "Build the source with the harmonic trap instead of GOST" OFF )
option( BUILD_VARIANT2 "Variant 2 - Choice of the method how the reference point is computed" ON )
option( BUILD_NEHARI "Disables the initial point in function space being set to the Nehari value" ON )
option( BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

set( HOME $ENV{HOME} CACHE STRING INTERNAL )
set( DIR_INC ${PROJECT_SOURCE_DIR}/include CACHE STRING INTERNAL )

if(BUILD_DOCUMENTATION)
    if(NOT DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()

    set(doxyfile ${PROJECT_SOURCE_DIR}/Doxyfile)

    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif()

### GSL library paths
find_library( GSL_LIBRARY_1 NAMES gsl  PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT GSL_LIBRARY_1 )
  message( FATAL_ERROR "Could not find the gsl library." )
endif()
find_library( GSL_LIBRARY_2 NAMES gslcblas PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT GSL_LIBRARY_2 )
  message( FATAL_ERROR "Could not find the gsl library." )
endif()
get_filename_component( TMP ${GSL_LIBRARY_1} PATH )
get_filename_component( TMP ${TMP} PATH )
set( GSL_INC_DIR ${TMP}/include CACHE STRING INTERNAL )

### p4est library paths
find_library( P4EST_LIBRARY NAMES p4est PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT P4EST_LIBRARY )
  message( FATAL_ERROR "Could not find the p4est library." )
endif()
get_filename_component( TMP ${P4EST_LIBRARY} PATH )
get_filename_component( TMP ${TMP} PATH )
set( P4EST_INC_DIR ${TMP}/include CACHE STRING INTERNAL )

### PETSc library paths
find_library( PETSC_LIBRARY NAMES petsc  PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT PETSC_LIBRARY )
  message( FATAL_ERROR "Could not find the PETSc library." )
endif()
get_filename_component( TMP ${PETSC_LIBRARY} PATH )
get_filename_component( TMP ${TMP} PATH )
set( PETSC_INC_DIR ${TMP}/include CACHE STRING INTERNAL )

### deal.ii library paths
find_library( DEAL_II_LIBRARY NAMES deal_II  PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT DEAL_II_LIBRARY )
  message( FATAL_ERROR "Could not find the deal.ii library." )
endif()
get_filename_component( TMP ${DEAL_II_LIBRARY} PATH )
get_filename_component( TMP ${TMP} PATH )
set( DEAL_INC_DIR ${TMP}/include CACHE STRING INTERNAL )

set( EXECUTABLE_OUTPUT_PATH  ${HOME}/bin )
set( LIBRARY_OUTPUT_PATH  ${HOME}/bin/lib )
set( CMAKE_INSTALL_PREFIX ${HOME}/bin CACHE PATH "Prefix prepended to install directories" FORCE ) 

if( BUILD_HTRAP )
  add_definitions( -DPOTENTIAL=1 )
else()
  add_definitions( -DPOTENTIAL=2 )
endif()

if( BUILD_3D )
  add_definitions( -DDIMENSION=3 )
else()
  add_definitions( -DDIMENSION=2 )
endif()

if( BUILD_VARIANT2 )
  add_definitions( -D__variant_2__ )
else()
  add_definitions( -D__variant_1__ )
endif()

if( BUILD_NEHARI )
  add_definitions( -DNEHARI )
endif()

include_directories( ${MPI_CXX_INCLUDE_PATH} ${DIR_INC} ${DEAL_INC_DIR} ${DEAL_INC_DIR}/deal.II/bundled ${GSL_INC_DIR} ${P4EST_INC_DIR} ${PETSC_INC_DIR} )

add_subdirectory( src )
