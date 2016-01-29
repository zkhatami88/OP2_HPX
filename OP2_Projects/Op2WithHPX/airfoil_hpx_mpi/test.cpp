

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <hpx/hpx_init.hpp>
#include <hpx/hpx_main.hpp>
#include <mpi.h>

#include  "op_lib_cpp.h"
#include "op_lib_mpi.h"


int main(int argc, char **argv)
{
  int my_rank;
  int comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  return 0;
}
