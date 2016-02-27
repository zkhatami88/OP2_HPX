#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/include/parallel_algorithm.hpp>
#include <hpx/include/parallel_executor_parameters.hpp>
#include <hpx/include/iostreams.hpp>

typedef struct
{
    int         index;  /* index */
    op_dat      dat;    /* dataset */
    op_map      map;    /* indirect mapping */
    int         dim,    /* dimension of data */
    idx,
    size;   /* size (for sequential execution) */
    char       *data,   /* data on host */
    *data_d; /* data on device (for CUDA execution) */
    int        *map_data,   /* data on host */
    *map_data_d; /* data on device (for CUDA execution) */
    char const *type;   /* datatype */
    op_access   acc;
    op_arg_type argtype;
    int         sent;   /* flag to indicate if this argument has
                         data in flight under non-blocking MPI comms*/
    int         opt;    /* flag to indicate if this argument is in use */
} op_arg;
