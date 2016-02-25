//
// auto-generated by op2.py
//

//user function
#include "save_soln.h"
#include <vector>
#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/include/parallel_algorithm.hpp>
#include <hpx/include/parallel_executor_parameters.hpp>
#include <hpx/include/iostreams.hpp>

// host stub function
hpx::future<void> op_par_loop_save_soln(hpx::future<char const *> name, hpx::future<op_set> set,
  hpx::future<op_arg> arg0,
  hpx::future<op_arg> arg1)
{

   using hpx::lcos::local::dataflow;
   using hpx::util::unwrapped;

  return dataflow(unwrapped([](char const * name, op_set set, op_arg arg0, op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(0);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  save_soln");
  }

  op_mpi_halo_exchanges(set, nargs, args);
  // set number of threads
  #ifdef _OPENMP
    int nthreads = omp_get_max_threads();
  #else
    int nthreads = 1;
  #endif

  hpx::parallel::dynamic_chunk_size dcs(500);
  if (set->size >0) {

    // execute plan
    auto r=boost::irange(0, nthreads);

    hpx::parallel::for_each(hpx::parallel::par.with(dcs),r.begin(), r.end(),[&](std::size_t thr){
    int start  = (set->size* thr)/nthreads;
    int finish = (set->size*(thr+1))/nthreads;
    for ( int n=start; n<finish; n++ ){

      save_soln(
        &((double*)arg0.data)[4*n],
        &((double*)arg1.data)[4*n]);

    }
    });


  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[0].name      = name;
  OP_kernels[0].count    += 1;
  OP_kernels[0].time     += wall_t2 - wall_t1;
  OP_kernels[0].transfer += (float)set->size * arg0.size;
  OP_kernels[0].transfer += (float)set->size * arg1.size;


}),name,set,arg0,arg1);
}
