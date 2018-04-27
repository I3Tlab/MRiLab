matlab_add_mex(NAME DoKSpaceTrajVTK SRC DoKSpaceTrajVTK.cpp ${CMAKE_SOURCE_DIR}/Matlabdef.def LINK_TO ${Matlab_LIBRARIES})

if(VTK_LIBRARIES)
    target_link_libraries(DoKSpaceTrajVTK ${VTK_LIBRARIES})
else(VTK_LIBRARIES)
    target_link_libraries(DoKSpaceTrajVTK vtkHybrid)
endif(VTK_LIBRARIES)

install(TARGETS DoKSpaceTrajVTK DESTINATION ../bin)