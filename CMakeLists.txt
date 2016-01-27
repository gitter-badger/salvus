cmake_minimum_required(VERSION 3.2.2)
project(salvus)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
SET(CMAKE_C_COMPILER /usr/local/bin/mpicc)
SET(CMAKE_CXX_COMPILER /usr/local/bin/mpicxx)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake-modules")

SET(PETSC_DIR "/opt/petsc")
find_package(PETSc REQUIRED)
include_directories(${PETSC_DIR}/include)

include_directories("./inc")
add_executable(salvus src/cxx/salvus.cpp)
target_link_libraries(salvus ${PETSC_LIBRARY_SINGLE})