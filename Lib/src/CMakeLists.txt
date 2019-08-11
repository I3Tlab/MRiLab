
cmake_minimum_required(VERSION 2.8)
project(MRiLab)

# Ubuntu Linux (local user install)
SET(ENV{IPP_ROOT} /home.old/arwillis/intel)
SET(ENV{MATLAB_ROOT} /usr/local/bin/matlab/R2018a)
# On Mac OS X
#SET(ENV{IPP_ROOT} /opt/intel)
#SET(ENV{MATLAB_ROOT} /Applications/MATLAB_R2016b.app)

SET(CMAKE_VERBOSE_MAKEFILE ON)
#MESSAGE( STATUS "WIN32: " ${WIN32} )
#MESSAGE( STATUS "UNIX: " ${UNIX} )
#MESSAGE( STATUS "APPLE: " ${APPLE} )

if (APPLE)
SET(LLVM_DIR /usr/local/opt/llvm/share/cmake/modules/)
SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
SET(MACOSX_RPATH 1)
SET(CC /usr/local/opt/llvm/bin/clang)
SET(CXX /usr/local/opt/llvm/bin/clang++)
SET(CMAKE_CXX_COMPILER g++-7)
endif(APPLE)


#set(CC gcc)
#set(CXX g++)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_SIZEOF_VOID_P 8) # use 64-bit by default
set(CMAKE_BUILD_TYPE Release) # use Release version by default

add_definitions(/DMATLAB_MEX_FILE) #define matlab macros
add_definitions(/DMX_COMPAT_32)

find_package(Matlab REQUIRED)
find_package(IPP)
find_package(Framewave)

if (APPLE)
  set(LDFLAGS -L/usr/local/opt/llvm/lib/)
  find_package(LLVM)
else(APPLE)
  find_package(OpenMP)
endif(APPLE)

find_package(CUDA)
find_package(VTK)
find_package(Ismrmrd)
find_package(Boost COMPONENTS thread system)

IF(MATLAB_FOUND)
    add_subdirectory(engine)
    add_subdirectory(recon)
    add_subdirectory(utilities)
    add_subdirectory(renderer)
    add_subdirectory(interface)
ELSE(MATLAB_FOUND)
    MESSAGE("MATLAB not found...nothing will be built.")
ENDIF(MATLAB_FOUND)



