matlab_add_mex(NAME DoSARAverageAtCPU SRC DoSARAverageAtCPU.cpp ${CMAKE_SOURCE_DIR}/Matlabdef.def LINK_TO ${Matlab_LIBRARIES})

install(TARGETS DoSARAverageAtCPU DESTINATION ../bin)