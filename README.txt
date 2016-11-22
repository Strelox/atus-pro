-------------------------------------------------------------------------------
content of the include directory

file   : breed_ortho_funcs.h
purpose: GSL functions for the reference point. 

file   : CBase.h
purpose: Base class for the multi processor breed solvers.

file   : CBase_no_mpi.h
purpose: Base class for breed_1.cpp

file   : config.h
purpose: The global degree of the FE_Q polynomials can be changed here.

file   : ef.h
purpose: Contains arrays of function pointers for the AiryAi functions, Hermite 
         and Legendre polynome. 

file   : functions_1.h
purpose: Contains the eigenfunctions and potentials for breed_1.cpp. 

file   : functions_cs.h
purpose: Contains the eigenfunctions and potentials for breed_cs.cpp. 

file   : functions.h
purpose: Contains the eigenfunctions and potentials for breed.cpp. 

file   : my_table.h
purpose: header file for my_table class

file   : MyParameterHandler.h
purpose: header file for the MyParameterHandler class

file   : ref_pt_list.h
purpose: A special class that manages the reference points. In general a 
         reference point may not be unique. 

file   : strtk.hpp
purpose: Third party template class for extendend std::string handling.

file   : tinyxml2.h
purpose: Include file of third party c++ classes for handling xml files.

-------------------------------------------------------------------------------
content of the src directory

file   : breed_1.cpp
purpose: Single core solver for the one dimensional 
         time independent Gross-Piteavskii equation for 
         real solutions. Interpolation from real to complex
         solutions is conducted afer each newly found solution.

file   : breedc.cpp
purpose: Multi processor solver for the time independent 
         Gross-Piteavskii equation for two or three 
         spatial dimensions in cartesian coordinates for 
         complex solutions.

file   : breed.cpp
purpose: Multi processor solver for the time independent 
         Gross-Piteavskii equation for two or three 
         spatial dimensions in cartesian coordinates for
         real solutions. Interpolation from real to complex
         solutions is conducted afer each newly found solution.

file   : breed_cs.cpp
purpose: Multi processor solver for the time independent 
         Gross-Piteavskii equation in cylinder symmetry for real
         solutions. Interpolation from real to complex
         solutions is conducted afer each newly found solution.

file   : breed_sob_1.cpp
purpose: Single core solver for the one dimensional time independent 
         Gross-Piteavskii equation for real solutions using the Sobolev 
         gradient method. Interpolation from real to complex solutions 
         is conducted afer a new solution is found.

file   : breed_sob.cpp
purpose: Multi processor solver for the time independent 
         Gross-Piteavskii equation for two or three 
         spatial dimensions in cartesian coordinates
         for real solutions using the Sobolev gradient method. 
         Interpolation from real to complex solutions is conducted 
         afer a new solution is found.

file   : CMakeLists.txt
purpose: Required by the CMake build system. 

file   : gen_params.cpp
purpose: Generates a set of parameter files located in sub folders named 
         by the quantum numbers of the initial guess. The parameter files
         are supposed to be used by breed. 

file   : gen_params_cs.cpp
purpose: Same as above, but for breed_cs. 

file   : my_table.cpp
purpose: Contains a simple table class for double values.

file   : MyParameterHandler.cpp
purpose: Contains the class, which handles the parameter files.

file   : rtprop_1.cpp
purpose: Real time propagation of the one dimensional Gross-Piteavskii
         equation. 

file   : rtprop.cpp
purpose: Real time propagation of the Gross-Piteavskii equation for 
         two or three spatial dimensions in cartesian coordinates.

file   : rtprop_cs.cpp
purpose: Real time propagation of the Gross-Piteavskii equation for 
         two or three spatial dimensions in cartesian coordinates.

file   : tinyxml2.cpp
purpose: Implementation file for the third party tinyxml2 classes.
