//
// header for sequential and MPI+sequentional execution
//

#include "op_lib_cpp.h"

//#define op_malloc malloc

static int op2_stride = 1;
#define OP2_STRIDE(arr, idx) arr[idx]

// scratch space to use for double counting in indirect reduction
static int blank_args_size = 512;
static char* blank_args = (char *)op_malloc(blank_args_size);

inline void op_arg_set(int n, op_arg arg, char **p_arg, int halo){
    *p_arg = arg.data;
    
    if (arg.argtype==OP_ARG_GBL) {
        if (halo && (arg.acc != OP_READ)) *p_arg = blank_args;
    }
    else {
        if (arg.map==NULL || arg.opt==0)         // identity mapping
            *p_arg += arg.size*n;
        else                       // standard pointers
            *p_arg += arg.size*arg.map->map[arg.idx+n*arg.map->dim];
    }
}

inline void op_arg_copy_in(int n, op_arg arg, char **p_arg) {
    for (int i = 0; i < -1*arg.idx; ++i)
        p_arg[i] = arg.data + arg.map->map[i+n*arg.map->dim]*arg.size;
}

inline void op_args_check(op_set set, int nargs, op_arg *args,
                          int *ninds, const char *name) {
    for (int n=0; n<nargs; n++)
        op_arg_check(set,n,args[n],ninds,name);
}

//
//op_par_loop routine for 1 arguments
//
template <class T0>
std::vector<std::vector<std::vector<hpx::shared_future<double>>>> op_par_loop(
//void op_par_loop(
                 void (*kernel)(T0*),
                 char const * name, op_set set,
                 op_arg arg0){
    
    char *p_a[1] = {0};
    op_arg args[1] = {arg0};
    if(arg0.idx < -1) {
        p_a[0] = (char *)op_malloc(-1*args[0].idx*sizeof(T0));
    }
    
    //allocate scratch mememory to do double counting in indirect reduction
    for (int i = 0; i<1;i++)
        if(args[i].argtype == OP_ARG_GBL && args[i].size > blank_args_size )
        {
            blank_args_size = args[i].size;
            blank_args = (char *)op_malloc(blank_args_size);
        }
    // consistency checks
    int ninds = 0;
    if (OP_diags>0) op_args_check(set,1,args,&ninds,name);
    
    if (OP_diags>2) {
        if (ninds==0)
            printf(" kernel routine w/o indirection:  %s\n",name);
        else
            printf(" kernel routine with indirection: %s\n",name);
    }
    // initialise timers
    double cpu_t1, cpu_t2, wall_t1, wall_t2;
    op_timers_core(&cpu_t1, &wall_t1);
    
    // MPI halo exchange and dirty bit setting, if needed
    int n_upper = op_mpi_halo_exchanges(set, 1, args);
    std::vector<std::vector<std::vector<hpx::shared_future<double>>>> new_data(n_upper); //**
    
    // loop over set elements
    int halo = 0;
    
    for (int n=0; n<n_upper; n++) {
        if (n==set->core_size) op_mpi_wait_all(1,args);
        if (n==set->size) halo = 1;
        if (args[0].idx < -1) op_arg_copy_in(n,args[0], (char **)p_a[0]);
        else op_arg_set(n,args[0], &p_a[0],halo);
        
        new_data[i] = hpx::async(kernel( (T0 *)p_a[0]));
    }
    if ( n_upper == set->core_size || n_upper == 0 )
        op_mpi_wait_all (1,args);
    
    //set dirty bit on datasets touched
    op_mpi_set_dirtybit(1, args);
    
    //global reduction for MPI execution, if needed
    //p_a simply used to determine type for MPI reduction
    op_mpi_reduce(&arg0,(T0 *)p_a[0]);
    
    // update timer record
    op_timers_core(&cpu_t2, &wall_t2);
#ifdef COMM_PERF
    void *k_i = op_mpi_perf_time(name, wall_t2 - wall_t1);
    op_mpi_perf_comms(k_i, 1, args);
#else
    op_mpi_perf_time(name, wall_t2 - wall_t1);
#endif
    
    if(arg0.idx < -1) {
        free(p_a[0]);
    }
    
    return new_data;
}
//
//op_par_loop routine for 2 arguments
//
template <class T0,class T1>
std::vector<std::vector<std::vector<hpx::shared_future<double>>>> op_par_loop( //**
//void op_par_loop(
                 void (*kernel)(T0*, T1*),
                 char const * name, op_set set,
                 op_arg arg0, op_arg arg1){
    
    char *p_a[2] = {0,0};
    op_arg args[2] = {arg0, arg1};
    if(arg0.idx < -1) {
        p_a[0] = (char *)op_malloc(-1*args[0].idx*sizeof(T0));
    }
    if(arg1.idx < -1) {
        p_a[1] = (char *)op_malloc(-1*args[1].idx*sizeof(T1));
    }
    
    //allocate scratch mememory to do double counting in indirect reduction
    for (int i = 0; i<2;i++)
        if(args[i].argtype == OP_ARG_GBL && args[i].size > blank_args_size )
        {
            blank_args_size = args[i].size;
            blank_args = (char *)op_malloc(blank_args_size);
        }
    // consistency checks
    int ninds = 0;
    if (OP_diags>0) op_args_check(set,2,args,&ninds,name);
    
    if (OP_diags>2) {
        if (ninds==0)
            printf(" kernel routine w/o indirection:  %s\n",name);
        else
            printf(" kernel routine with indirection: %s\n",name);
    }
    // initialise timers
    double cpu_t1, cpu_t2, wall_t1, wall_t2;
    op_timers_core(&cpu_t1, &wall_t1);
    
    // MPI halo exchange and dirty bit setting, if needed
    int n_upper = op_mpi_halo_exchanges(set, 2, args);
    std::vector<std::vector<std::vector<hpx::shared_future<double>>>> new_data(n_upper); //**
    
    // loop over set elements
    int halo = 0;
    
    for (int n=0; n<n_upper; n++) {
        if (n==set->core_size) op_mpi_wait_all(2,args);
        if (n==set->size) halo = 1;
        if (args[0].idx < -1) op_arg_copy_in(n,args[0], (char **)p_a[0]);
        else op_arg_set(n,args[0], &p_a[0],halo);
        if (args[1].idx < -1) op_arg_copy_in(n,args[1], (char **)p_a[1]);
        else op_arg_set(n,args[1], &p_a[1],halo);
        
        new_data[i] = hpx::async(kernel( (T0 *)p_a[0], (T1 *)p_a[1]));
    }
    if ( n_upper == set->core_size || n_upper == 0 )
        op_mpi_wait_all (2,args);
    
    //set dirty bit on datasets touched
    op_mpi_set_dirtybit(2, args);
    
    //global reduction for MPI execution, if needed
    //p_a simply used to determine type for MPI reduction
    op_mpi_reduce(&arg0,(T0 *)p_a[0]);
    op_mpi_reduce(&arg1,(T1 *)p_a[1]);
    
    // update timer record
    op_timers_core(&cpu_t2, &wall_t2);
#ifdef COMM_PERF
    void *k_i = op_mpi_perf_time(name, wall_t2 - wall_t1);
    op_mpi_perf_comms(k_i, 2, args);
#else
    op_mpi_perf_time(name, wall_t2 - wall_t1);
#endif
    
    if(arg0.idx < -1) {
        free(p_a[0]);
    }
    if(arg1.idx < -1) {
        free(p_a[1]);
    }
    
    return new_data;
}
//
//op_par_loop routine for 3 arguments
//
template <class T0,class T1,class T2>
std::vector<std::vector<std::vector<hpx::shared_future<double>>>> op_par_loop(
//void op_par_loop(
                 void (*kernel)(T0*, T1*, T2*),
                 char const * name, op_set set,
                 op_arg arg0, op_arg arg1, op_arg arg2){
    
    char *p_a[3] = {0,0,0};
    op_arg args[3] = {arg0, arg1, arg2};
    if(arg0.idx < -1) {
        p_a[0] = (char *)op_malloc(-1*args[0].idx*sizeof(T0));
    }
    if(arg1.idx < -1) {
        p_a[1] = (char *)op_malloc(-1*args[1].idx*sizeof(T1));
    }
    if(arg2.idx < -1) {
        p_a[2] = (char *)op_malloc(-1*args[2].idx*sizeof(T2));
    }
    
    //allocate scratch mememory to do double counting in indirect reduction
    for (int i = 0; i<3;i++)
        if(args[i].argtype == OP_ARG_GBL && args[i].size > blank_args_size )
        {
            blank_args_size = args[i].size;
            blank_args = (char *)op_malloc(blank_args_size);
        }
    // consistency checks
    int ninds = 0;
    if (OP_diags>0) op_args_check(set,3,args,&ninds,name);
    
    if (OP_diags>2) {
        if (ninds==0)
            printf(" kernel routine w/o indirection:  %s\n",name);
        else
            printf(" kernel routine with indirection: %s\n",name);
    }
    // initialise timers
    double cpu_t1, cpu_t2, wall_t1, wall_t2;
    op_timers_core(&cpu_t1, &wall_t1);
    
    // MPI halo exchange and dirty bit setting, if needed
    int n_upper = op_mpi_halo_exchanges(set, 3, args);
    std::vector<std::vector<std::vector<hpx::shared_future<double>>>> new_data(n_upper); //**
    
    // loop over set elements
    int halo = 0;
    
    for (int n=0; n<n_upper; n++) {
        if (n==set->core_size) op_mpi_wait_all(3,args);
        if (n==set->size) halo = 1;
        if (args[0].idx < -1) op_arg_copy_in(n,args[0], (char **)p_a[0]);
        else op_arg_set(n,args[0], &p_a[0],halo);
        if (args[1].idx < -1) op_arg_copy_in(n,args[1], (char **)p_a[1]);
        else op_arg_set(n,args[1], &p_a[1],halo);
        if (args[2].idx < -1) op_arg_copy_in(n,args[2], (char **)p_a[2]);
        else op_arg_set(n,args[2], &p_a[2],halo);
        
        new_data[i] = hpx::async(kernel( (T0 *)p_a[0], (T1 *)p_a[1], (T2 *)p_a[2]));
    }
    if ( n_upper == set->core_size || n_upper == 0 )
        op_mpi_wait_all (3,args);
    
    //set dirty bit on datasets touched
    op_mpi_set_dirtybit(3, args);
    
    //global reduction for MPI execution, if needed
    //p_a simply used to determine type for MPI reduction
    op_mpi_reduce(&arg0,(T0 *)p_a[0]);
    op_mpi_reduce(&arg1,(T1 *)p_a[1]);
    op_mpi_reduce(&arg2,(T2 *)p_a[2]);
    
    // update timer record
    op_timers_core(&cpu_t2, &wall_t2);
#ifdef COMM_PERF
    void *k_i = op_mpi_perf_time(name, wall_t2 - wall_t1);
    op_mpi_perf_comms(k_i, 3, args);
#else
    op_mpi_perf_time(name, wall_t2 - wall_t1);
#endif
    
    if(arg0.idx < -1) {
        free(p_a[0]);
    }
    if(arg1.idx < -1) {
        free(p_a[1]);
    }
    if(arg2.idx < -1) {
        free(p_a[2]);
    }
    
    return new_data;
}
//
//op_par_loop routine for 4 arguments
//
template <class T0,class T1,class T2,class T3>
std::vector<std::vector<std::vector<hpx::shared_future<double>>>> op_par_loop(
//void op_par_loop(
                 void (*kernel)(T0*, T1*, T2*, T3*),
                 char const * name, op_set set,
                 op_arg arg0, op_arg arg1, op_arg arg2, op_arg arg3){
    
    char *p_a[4] = {0,0,0,0};
    op_arg args[4] = {arg0, arg1, arg2, arg3};
    if(arg0.idx < -1) {
        p_a[0] = (char *)op_malloc(-1*args[0].idx*sizeof(T0));
    }
    if(arg1.idx < -1) {
        p_a[1] = (char *)op_malloc(-1*args[1].idx*sizeof(T1));
    }
    if(arg2.idx < -1) {
        p_a[2] = (char *)op_malloc(-1*args[2].idx*sizeof(T2));
    }
    if(arg3.idx < -1) {
        p_a[3] = (char *)op_malloc(-1*args[3].idx*sizeof(T3));
    }
    
    //allocate scratch mememory to do double counting in indirect reduction
    for (int i = 0; i<4;i++)
        if(args[i].argtype == OP_ARG_GBL && args[i].size > blank_args_size )
        {
            blank_args_size = args[i].size;
            blank_args = (char *)op_malloc(blank_args_size);
        }
    // consistency checks
    int ninds = 0;
    if (OP_diags>0) op_args_check(set,4,args,&ninds,name);
    
    if (OP_diags>2) {
        if (ninds==0)
            printf(" kernel routine w/o indirection:  %s\n",name);
        else
            printf(" kernel routine with indirection: %s\n",name);
    }
    // initialise timers
    double cpu_t1, cpu_t2, wall_t1, wall_t2;
    op_timers_core(&cpu_t1, &wall_t1);
    
    // MPI halo exchange and dirty bit setting, if needed
    int n_upper = op_mpi_halo_exchanges(set, 4, args);
    std::vector<std::vector<std::vector<hpx::shared_future<double>>>> new_data(n_upper); //**
    
    // loop over set elements
    int halo = 0;
    
    for (int n=0; n<n_upper; n++) {
        if (n==set->core_size) op_mpi_wait_all(4,args);
        if (n==set->size) halo = 1;
        if (args[0].idx < -1) op_arg_copy_in(n,args[0], (char **)p_a[0]);
        else op_arg_set(n,args[0], &p_a[0],halo);
        if (args[1].idx < -1) op_arg_copy_in(n,args[1], (char **)p_a[1]);
        else op_arg_set(n,args[1], &p_a[1],halo);
        if (args[2].idx < -1) op_arg_copy_in(n,args[2], (char **)p_a[2]);
        else op_arg_set(n,args[2], &p_a[2],halo);
        if (args[3].idx < -1) op_arg_copy_in(n,args[3], (char **)p_a[3]);
        else op_arg_set(n,args[3], &p_a[3],halo);
        
        new_data[i] = hpx::async(kernel( (T0 *)p_a[0], (T1 *)p_a[1], (T2 *)p_a[2], (T3 *)p_a[3]));
    }
    if ( n_upper == set->core_size || n_upper == 0 )
        op_mpi_wait_all (4,args);
    
    //set dirty bit on datasets touched
    op_mpi_set_dirtybit(4, args);
    
    //global reduction for MPI execution, if needed
    //p_a simply used to determine type for MPI reduction
    op_mpi_reduce(&arg0,(T0 *)p_a[0]);
    op_mpi_reduce(&arg1,(T1 *)p_a[1]);
    op_mpi_reduce(&arg2,(T2 *)p_a[2]);
    op_mpi_reduce(&arg3,(T3 *)p_a[3]);
    
    // update timer record
    op_timers_core(&cpu_t2, &wall_t2);
#ifdef COMM_PERF
    void *k_i = op_mpi_perf_time(name, wall_t2 - wall_t1);
    op_mpi_perf_comms(k_i, 4, args);
#else
    op_mpi_perf_time(name, wall_t2 - wall_t1);
#endif
    
    if(arg0.idx < -1) {
        free(p_a[0]);
    }
    if(arg1.idx < -1) {
        free(p_a[1]);
    }
    if(arg2.idx < -1) {
        free(p_a[2]);
    }
    if(arg3.idx < -1) {
        free(p_a[3]);
    }
    
    return new_data;
}
//
//op_par_loop routine for 5 arguments
//
template <class T0,class T1,class T2,class T3,
class T4>
std::vector<std::vector<std::vector<hpx::shared_future<double>>>> op_par_loop(
//void op_par_loop(
                 void (*kernel)(T0*, T1*, T2*, T3*,
                                T4*),
                 char const * name, op_set set,
                 op_arg arg0, op_arg arg1, op_arg arg2, op_arg arg3,
                 op_arg arg4){
    
    char *p_a[5] = {0,0,0,0,0};
    op_arg args[5] = {arg0, arg1, arg2, arg3,
        arg4};
    if(arg0.idx < -1) {
        p_a[0] = (char *)op_malloc(-1*args[0].idx*sizeof(T0));
    }
    if(arg1.idx < -1) {
        p_a[1] = (char *)op_malloc(-1*args[1].idx*sizeof(T1));
    }
    if(arg2.idx < -1) {
        p_a[2] = (char *)op_malloc(-1*args[2].idx*sizeof(T2));
    }
    if(arg3.idx < -1) {
        p_a[3] = (char *)op_malloc(-1*args[3].idx*sizeof(T3));
    }
    if(arg4.idx < -1) {
        p_a[4] = (char *)op_malloc(-1*args[4].idx*sizeof(T4));
    }
    
    //allocate scratch mememory to do double counting in indirect reduction
    for (int i = 0; i<5;i++)
        if(args[i].argtype == OP_ARG_GBL && args[i].size > blank_args_size )
        {
            blank_args_size = args[i].size;
            blank_args = (char *)op_malloc(blank_args_size);
        }
    // consistency checks
    int ninds = 0;
    if (OP_diags>0) op_args_check(set,5,args,&ninds,name);
    
    if (OP_diags>2) {
        if (ninds==0)
            printf(" kernel routine w/o indirection:  %s\n",name);
        else
            printf(" kernel routine with indirection: %s\n",name);
    }
    // initialise timers
    double cpu_t1, cpu_t2, wall_t1, wall_t2;
    op_timers_core(&cpu_t1, &wall_t1);
    
    // MPI halo exchange and dirty bit setting, if needed
    int n_upper = op_mpi_halo_exchanges(set, 5, args);
    std::vector<std::vector<std::vector<hpx::shared_future<double>>>> new_data(n_upper); //**
    
    // loop over set elements
    int halo = 0;
    
    for (int n=0; n<n_upper; n++) {
        if (n==set->core_size) op_mpi_wait_all(5,args);
        if (n==set->size) halo = 1;
        if (args[0].idx < -1) op_arg_copy_in(n,args[0], (char **)p_a[0]);
        else op_arg_set(n,args[0], &p_a[0],halo);
        if (args[1].idx < -1) op_arg_copy_in(n,args[1], (char **)p_a[1]);
        else op_arg_set(n,args[1], &p_a[1],halo);
        if (args[2].idx < -1) op_arg_copy_in(n,args[2], (char **)p_a[2]);
        else op_arg_set(n,args[2], &p_a[2],halo);
        if (args[3].idx < -1) op_arg_copy_in(n,args[3], (char **)p_a[3]);
        else op_arg_set(n,args[3], &p_a[3],halo);
        if (args[4].idx < -1) op_arg_copy_in(n,args[4], (char **)p_a[4]);
        else op_arg_set(n,args[4], &p_a[4],halo);
        
        new_data[i] = hpx::async(kernel( (T0 *)p_a[0], (T1 *)p_a[1], (T2 *)p_a[2], (T3 *)p_a[3],
               (T4 *)p_a[4]));
    }
    if ( n_upper == set->core_size || n_upper == 0 )
        op_mpi_wait_all (5,args);
    
    //set dirty bit on datasets touched
    op_mpi_set_dirtybit(5, args);
    
    //global reduction for MPI execution, if needed
    //p_a simply used to determine type for MPI reduction
    op_mpi_reduce(&arg0,(T0 *)p_a[0]);
    op_mpi_reduce(&arg1,(T1 *)p_a[1]);
    op_mpi_reduce(&arg2,(T2 *)p_a[2]);
    op_mpi_reduce(&arg3,(T3 *)p_a[3]);
    op_mpi_reduce(&arg4,(T4 *)p_a[4]);
    
    // update timer record
    op_timers_core(&cpu_t2, &wall_t2);
#ifdef COMM_PERF
    void *k_i = op_mpi_perf_time(name, wall_t2 - wall_t1);
    op_mpi_perf_comms(k_i, 5, args);
#else
    op_mpi_perf_time(name, wall_t2 - wall_t1);
#endif
    
    if(arg0.idx < -1) {
        free(p_a[0]);
    }
    if(arg1.idx < -1) {
        free(p_a[1]);
    }
    if(arg2.idx < -1) {
        free(p_a[2]);
    }
    if(arg3.idx < -1) {
        free(p_a[3]);
    }
    if(arg4.idx < -1) {
        free(p_a[4]);
    }
    
    return new_data;
}
//
//op_par_loop routine for 6 arguments
//
template <class T0,class T1,class T2,class T3,
class T4,class T5>
std::vector<std::vector<std::vector<hpx::shared_future<double>>>> op_par_loop(
//void op_par_loop(
                 void (*kernel)(T0*, T1*, T2*, T3*,
                                T4*, T5*),
                 char const * name, op_set set,
                 op_arg arg0, op_arg arg1, op_arg arg2, op_arg arg3,
                 op_arg arg4, op_arg arg5){
    
    char *p_a[6] = {0,0,0,0,0,0};
    op_arg args[6] = {arg0, arg1, arg2, arg3,
        arg4, arg5};
    if(arg0.idx < -1) {
        p_a[0] = (char *)op_malloc(-1*args[0].idx*sizeof(T0));
    }
    if(arg1.idx < -1) {
        p_a[1] = (char *)op_malloc(-1*args[1].idx*sizeof(T1));
    }
    if(arg2.idx < -1) {
        p_a[2] = (char *)op_malloc(-1*args[2].idx*sizeof(T2));
    }
    if(arg3.idx < -1) {
        p_a[3] = (char *)op_malloc(-1*args[3].idx*sizeof(T3));
    }
    if(arg4.idx < -1) {
        p_a[4] = (char *)op_malloc(-1*args[4].idx*sizeof(T4));
    }
    if(arg5.idx < -1) {
        p_a[5] = (char *)op_malloc(-1*args[5].idx*sizeof(T5));
    }
    
    //allocate scratch mememory to do double counting in indirect reduction
    for (int i = 0; i<6;i++)
        if(args[i].argtype == OP_ARG_GBL && args[i].size > blank_args_size )
        {
            blank_args_size = args[i].size;
            blank_args = (char *)op_malloc(blank_args_size);
        }
    // consistency checks
    int ninds = 0;
    if (OP_diags>0) op_args_check(set,6,args,&ninds,name);
    
    if (OP_diags>2) {
        if (ninds==0)
            printf(" kernel routine w/o indirection:  %s\n",name);
        else
            printf(" kernel routine with indirection: %s\n",name);
    }
    // initialise timers
    double cpu_t1, cpu_t2, wall_t1, wall_t2;
    op_timers_core(&cpu_t1, &wall_t1);
    
    // MPI halo exchange and dirty bit setting, if needed
    int n_upper = op_mpi_halo_exchanges(set, 6, args);
    std::vector<std::vector<std::vector<hpx::shared_future<double>>>> new_data(n_upper); //**
    
    // loop over set elements
    int halo = 0;
    
    for (int n=0; n<n_upper; n++) {
        if (n==set->core_size) op_mpi_wait_all(6,args);
        if (n==set->size) halo = 1;
        if (args[0].idx < -1) op_arg_copy_in(n,args[0], (char **)p_a[0]);
        else op_arg_set(n,args[0], &p_a[0],halo);
        if (args[1].idx < -1) op_arg_copy_in(n,args[1], (char **)p_a[1]);
        else op_arg_set(n,args[1], &p_a[1],halo);
        if (args[2].idx < -1) op_arg_copy_in(n,args[2], (char **)p_a[2]);
        else op_arg_set(n,args[2], &p_a[2],halo);
        if (args[3].idx < -1) op_arg_copy_in(n,args[3], (char **)p_a[3]);
        else op_arg_set(n,args[3], &p_a[3],halo);
        if (args[4].idx < -1) op_arg_copy_in(n,args[4], (char **)p_a[4]);
        else op_arg_set(n,args[4], &p_a[4],halo);
        if (args[5].idx < -1) op_arg_copy_in(n,args[5], (char **)p_a[5]);
        else op_arg_set(n,args[5], &p_a[5],halo);
        
        new_data[i] = hpx::async(kernel( (T0 *)p_a[0], (T1 *)p_a[1], (T2 *)p_a[2], (T3 *)p_a[3],
               (T4 *)p_a[4], (T5 *)p_a[5]));
    }
    if ( n_upper == set->core_size || n_upper == 0 )
        op_mpi_wait_all (6,args);
    
    //set dirty bit on datasets touched
    op_mpi_set_dirtybit(6, args);
    
    //global reduction for MPI execution, if needed
    //p_a simply used to determine type for MPI reduction
    op_mpi_reduce(&arg0,(T0 *)p_a[0]);
    op_mpi_reduce(&arg1,(T1 *)p_a[1]);
    op_mpi_reduce(&arg2,(T2 *)p_a[2]);
    op_mpi_reduce(&arg3,(T3 *)p_a[3]);
    op_mpi_reduce(&arg4,(T4 *)p_a[4]);
    op_mpi_reduce(&arg5,(T5 *)p_a[5]);
    
    // update timer record
    op_timers_core(&cpu_t2, &wall_t2);
#ifdef COMM_PERF
    void *k_i = op_mpi_perf_time(name, wall_t2 - wall_t1);
    op_mpi_perf_comms(k_i, 6, args);
#else
    op_mpi_perf_time(name, wall_t2 - wall_t1);
#endif
    
    if(arg0.idx < -1) {
        free(p_a[0]);
    }
    if(arg1.idx < -1) {
        free(p_a[1]);
    }
    if(arg2.idx < -1) {
        free(p_a[2]);
    }
    if(arg3.idx < -1) {
        free(p_a[3]);
    }
    if(arg4.idx < -1) {
        free(p_a[4]);
    }
    if(arg5.idx < -1) {
        free(p_a[5]);
    }
    
    return new_data;
}