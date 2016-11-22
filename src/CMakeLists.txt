
### parameter generators

ADD_EXECUTABLE( gen_params gen_params.cpp MyParameterHandler.cpp tinyxml2.cpp )
TARGET_LINK_LIBRARIES( gen_params ${MPI_CXX_LIBRARIES} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( gen_params_cs gen_params_cs.cpp MyParameterHandler.cpp  tinyxml2.cpp )
TARGET_LINK_LIBRARIES( gen_params_cs ${MPI_CXX_LIBRARIES} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

### stationary solver my newton

ADD_EXECUTABLE( breed_1 breed_1.cpp MyParameterHandler.cpp tinyxml2.cpp my_table.cpp )
TARGET_LINK_LIBRARIES( breed_1 ${MPI_CXX_LIBRARIES} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( breed breed.cpp MyParameterHandler.cpp tinyxml2.cpp my_table.cpp )
TARGET_LINK_LIBRARIES( breed ${MPI_CXX_LIBRARIES} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( breedc breedc.cpp MyParameterHandler.cpp tinyxml2.cpp my_table.cpp )
TARGET_LINK_LIBRARIES( breedc ${MPI_CXX_LIBRARIES} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( breed_cs breed_cs.cpp MyParameterHandler.cpp tinyxml2.cpp my_table.cpp )
TARGET_LINK_LIBRARIES( breed_cs ${MPI_CXX_LIBRARIES} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

### stationary solver sobolev gradient

ADD_EXECUTABLE( breed_sob_1 breed_sob_1.cpp MyParameterHandler.cpp tinyxml2.cpp my_table.cpp )
TARGET_LINK_LIBRARIES( breed_sob_1 ${MPI_CXX_LIBRARIES} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( breed_sob breed_sob.cpp MyParameterHandler.cpp tinyxml2.cpp my_table.cpp )
TARGET_LINK_LIBRARIES( breed_sob ${MPI_CXX_LIBRARIES} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

### real time propagation

ADD_EXECUTABLE( rtprop_1 rtprop_1.cpp MyParameterHandler.cpp tinyxml2.cpp my_table.cpp )
TARGET_LINK_LIBRARIES( rtprop_1 ${MPI_CXX_LIBRARIES} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( rtprop rtprop.cpp MyParameterHandler.cpp tinyxml2.cpp my_table.cpp )
TARGET_LINK_LIBRARIES( rtprop ${MPI_CXX_LIBRARIES} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( rtprop_cs rtprop_cs.cpp MyParameterHandler.cpp tinyxml2.cpp my_table.cpp )
TARGET_LINK_LIBRARIES( rtprop_cs ${MPI_CXX_LIBRARIES} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )
