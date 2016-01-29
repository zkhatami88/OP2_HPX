cmake_minimum_required(VERSION 2.8)
set(HPX_IGNORE_COMPILER_COMPATIBILITY 1)
include(common.cmake)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -g -fPIC -D_REENTRANT  -DHAVE_PARMETIS -DPARMETIS_VER_4  -DHAVE_PTSCOTCH -DUNIX -Wall -DOP_PART_SIZE_1=512 -DOP_PART_SIZE_2=1024 -DOP_PART_SIZE_3=64 -DOP_BLOCK_SIZE_0=512 -DOP_BLOCK_SIZE_1=64 -DOP_BLOCK_SIZE_2=64 -DOP_BLOCK_SIZE_3=64 -g -fPIC -DUNIX -Wall -fopenmp -I/home/zahra/OP2-Common/op2/c/include -L/home/zahra/OP2-Common/op2/c/lib -I/usr/bin/parmetis-4.0.3/include -L/usr/bin/parmetis-4.0.3/libparmetis -lptscotch -lptscotcherr  -lparmetis -lmetis -lmpi -lop2_mpi") 

set(HPX_INCLUDE_DIR "/home/zahra/Projects/HPX/repo/hpx/include")
set(HPX_LIBRARY_DIR "/home/zahra/Projects/HPX/build/lib")
set(OP2_INC "/home/zahra/OP2-Common/op2/c/include")
set(OP2_LIB "/home/zahra/OP2-Common/op2/c/lib")

find_package(HPX REQUIRED)
find_package(OP2)

include_directories(${HPX_INCLUDE_DIR})
include_directories(${MPI_INCLUDE_PATH})
include_directories(/home/zahra/OP2-Common/op2/c/include)

link_directories(/home/zahra/OP2-Common/op2/c/src  ${MPI_LIBRARIES} ${HPX_LIBRARY_DIR}) 

set(SOURCE_FILES  airfoil_kernels.cpp airfoil_mpi_op.cpp /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/save_soln.h /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/adt_calc.h /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/res_calc.h /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/bres_calc.h /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/update.h /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/op_lib_c.h /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/op_lib_core.h /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/op_lib_cpp.h /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/op_rt_support.h /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/op_mpi_core.h /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/op_cuda_reduction.h /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/op_cuda_rt_support.h /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/op_hdf5.h /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/op_lib_mpi.h /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/op_openmp_rt_support.h /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/uthash.h /home/zahra/OP2-Common/apps/c/airfoil/myexample/airfoil_hpx_MPI/op_util.h  /home/zahra/OP2-Common/op2/c/src/core/op_rt_support.c  /home/zahra/OP2-Common/op2/c/src/core/op_lib_core.c  /home/zahra/OP2-Common/op2/c/src/core/op_dummy_singlenode.c  /home/zahra/OP2-Common/op2/c/src/openmp/op_openmp_decl.c)

add_hpx_executable(air_mpi
         ESSENTIAL
         SOURCES ${SOURCE_FILES}
         COMPONENT_DEPENDENCIES iostreams)




