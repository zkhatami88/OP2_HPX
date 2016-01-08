cmake_minimum_required(VERSION 2.8)
project(air_seq)

set(HPX_IGNORE_COMPILER_COMPATIBILITY 1)
set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(OP2_C_COMPILER "/usr/bin/gcc")
set(OP2_CXX_COMPILER "/usr/bin/g++")
include_directories(/home/zahra/OP2-Common/op2/c/include)
link_directories(/home/zahra/OP2-Common/op2/c/src)

include(common.cmake)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -fopenmp -I/home/zahra/OP2-Common/op2/c/include -L/home/zahra/OP2-Common/op2/c/lib -lop2_seq")
set(HPX_DIR "/home/zahra/Projects/HPX/build/hpx")
set(PROJECT_BINARY_DIR "/home/zahra/OP2-Common/apps/c/airfoil/myexample/FirstHPX")
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

set(OP2_INC "/home/zahra/OP2-Common/op2/c/include")
set(OP2_LIB "/home/zahra/OP2-Common/op2/c/lib")

find_package(HPX REQUIRED)

find_package(OP2)

include_directories(/home/zahra/OP2-Common/op2/c/include)
link_directories(/home/zahra/OP2-Common/op2/c/src)


set(SOURCE_FILES airfoil.cpp /home/zahra/OP2-Common/op2/c/src/core/op_rt_support.c /home/zahra/OP2-Common/op2/c/src/core/op_lib_core.c /home/zahra/OP2-Common/op2/c/src/core/op_dummy_singlenode.c  /home/zahra/OP2-Common/op2/c/src/sequential/op_seq.c)
#  /home/zahra/OP2-Common/op2/c/src/openmp/op_openmp_decl.c) remove se and add this for openmp version
add_hpx_executable(air_seq
        ESSENTIAL
        SOURCES ${SOURCE_FILES}
        COMPONENT_DEPENDENCIES iostreams)

include_directories(${OP2_SOURCE_DIR}/include)
