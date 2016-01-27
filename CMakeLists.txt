cmake_minimum_required(VERSION 3.2.2)
project(salvus)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
SET(CMAKE_C_COMPILER /usr/local/bin/mpicc)
SET(CMAKE_CXX_COMPILER /usr/local/bin/mpicxx)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake-modules")

SET(PETSC_DIR $ENV{PETSC_DIR})
SET(PETSC_ARCH $ENV{PETSC_ARCH})
find_package(PETSc REQUIRED CXX)
include_directories("$ENV{PETSC_DIR}/include")
include_directories("$ENV{PETSC_DIR}/${PETSC_ARCH}/include")

include_directories("./inc")
add_executable(salvus src/cxx/salvus.cpp src/cxx/Solver.cpp src/cxx/Solver.h src/cxx/Mesh.cpp src/cxx/Mesh.h)
target_link_libraries(salvus ${PETSC_LIBRARY_SINGLE})
