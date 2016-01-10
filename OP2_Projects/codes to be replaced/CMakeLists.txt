cmake_minimum_required(VERSION 2.8)
project(air_openmp)

set(HPX_IGNORE_COMPILER_COMPATIBILITY 1)
set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(OP2_C_COMPILER "/usr/bin/gcc")
set(OP2_CXX_COMPILER "/usr/bin/g++")
include_directories(/home/zahra/OP2-Common/op2/c/include)
link_directories(/home/zahra/OP2-Common/op2/c/src)

SET(VAR	“-DOP_PART_SIZE_1=512 -DOP_PART_SIZE_2=1024 -DOP_PART_SIZE_3=64 
#-DOP_BLOCK_SIZE_0=512 -DOP_BLOCK_SIZE_1=64 -DOP_BLOCK_SIZE_2=64 -DOP_BLOCK_SIZE_3=64 -#DOP_BLOCK_SIZE_4=64
”)

SET(SEQ_FLAGS "-std=c++11 g++ -g -fPIC -DUNIX -Wall /home/zahra/OP2-Common/op2/c/include /home/zahra/OP2-Common/op2/c/lib -lop2_seq")

SET(OPENMP_FLAGS "-std=c++11 g++ $(VAR) -g -fPIC -DUNIX -Wall -fopenmp /home/zahra/OP2-Common/op2/c/include /home/zahra/OP2-Common/op2/c/lib -lm -lop2_openmp")


include(common.cmake)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11
# -g -fopenmp -I/home/zahra/OP2-Common/op2/c/include
# -L/home/zahra/OP2-Common/op2/c/lib 
#-lop2_openmp") #if seq:-lop2_seq 

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMP_FLAGS}”) 
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQ_FLAGS}”) 

set(HPX_DIR "/home/zahra/Projects/HPX/build/hpx")
set(PROJECT_BINARY_DIR "/home/zahra/OP2-Common/apps/c/airfoil/myexample/FirstHPX")
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

set(OP2_INC "/home/zahra/OP2-Common/op2/c/include")
set(OP2_LIB "/home/zahra/OP2-Common/op2/c/lib")

find_package(HPX REQUIRED)

find_package(OP2)


include_directories(/home/zahra/OP2-Common/op2/c/include)
link_directories(/home/zahra/OP2-Common/op2/c/src)


set(SOURCE_FILES
airfoil_op.cpp airfoil_kernels.cpp
save_soln_kernel.cpp
adt_calc_kernel.cpp
res_calc_kernel.cpp
bres_calc_kernel.cpp
update_kernel.cpp 
# airfoil.cpp 
/home/zahra/OP2-Common/apps/c/airfoil/myexample/SecondHPX/save_soln.h
/home/zahra/OP2-Common/apps/c/airfoil/myexample/SecondHPX/adt_calc.h
/home/zahra/OP2-Common/apps/c/airfoil/myexample/SecondHPX/res_calc.h
/home/zahra/OP2-Common/apps/c/airfoil/myexample/SecondHPX/bres_calc.h
/home/zahra/OP2-Common/apps/c/airfoil/myexample/SecondHPX/update.h
/home/zahra/OP2-Common/op2/c/src/core/op_rt_support.c 
/home/zahra/OP2-Common/op2/c/src/core/op_lib_core.c 
/home/zahra/OP2-Common/op2/c/src/core/op_dummy_singlenode.c  
# /home/zahra/OP2-Common/op2/c/src/sequential/op_seq.c)
  /home/zahra/OP2-Common/op2/c/src/openmp/op_openmp_decl.c) 
add_hpx_executable(air_openmp
        ESSENTIAL
        SOURCES ${SOURCE_FILES}
        COMPONENT_DEPENDENCIES iostreams)

include_directories(${OP2_SOURCE_DIR}/include)
