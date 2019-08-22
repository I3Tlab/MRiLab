matlab_add_mex(NAME DoGridding SRC DoGridding.cpp ${CMAKE_SOURCE_DIR}/Matlabdef.def LINK_TO ${Matlab_LIBRARIES})

install(TARGETS DoGridding DESTINATION ../bin)
