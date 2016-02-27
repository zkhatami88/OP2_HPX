
//for direct loops, for indirect use op_arg_dat

#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/include/parallel_algorithm.hpp>
#include <hpx/include/parallel_executor_parameters.hpp>
#include <hpx/include/iostreams.hpp>

hpx::future<op_arg>
op_arg_dat1 ( p_dat dat, int idx, op_map map, int dim, const char * typ, op_access acc )
{
    
    using hpx::lcos::local::dataflow;
    using hpx::util::unwrapped;
    
    return dataflow(unwrapped([idx,map,dim,&typ,acc,dat](){
    op_arg arg;
    
    /* index is not used for now */
    arg.index = -1;
    arg.opt = 1;
    arg.argtype = OP_ARG_DAT;
    
    arg.dat = dat;
    arg.map = map;
    arg.dim = dim;
    arg.idx = idx;
    
    if ( dat != NULL )
    {
        arg.size = dat->size;
        arg.data = dat->data;
        arg.data_d = dat->data_d;
        arg.map_data_d = (map == NULL ? NULL : map->map_d);
        arg.map_data = (map == NULL ? NULL : map->map);
    }
    else
    {
        /* set default values */
        arg.size = -1;
        arg.data = NULL;
        arg.data_d = NULL;
        arg.map_data_d = NULL;
        arg.map_data = NULL;
    }
    
    arg.type = typ;
    arg.acc = acc;
    
    /*initialize to 0 states no-mpi messages inflight for this arg*/
    arg.sent = 0;
    
    return arg;

    });
}