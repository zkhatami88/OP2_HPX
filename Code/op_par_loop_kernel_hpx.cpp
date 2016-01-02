
//op_seq.h
// Making future output of kernel
//op_par_loop routine for 2 arguments
//
template <class T0,class T1>
std::vector<std::vector<std::vector<hpx::shared_future<double>>>> op_par_loop( //****add time here
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
    }c
    // initialise timers
    double cpu_t1, cpu_t2, wall_t1, wall_t2;
    op_timers_core(&cpu_t1, &wall_t1);
    
    // MPI halo exchange and dirty bit setting, if needed
    int n_upper = op_mpi_halo_exchanges(set, 2, args);
    
    std::vector<std::vector<std::vector<hpx::shared_future<double>>>> new_data(n_upper); //****
    
    // loop over set elements
    int halo = 0;
    
    for (int n=0; n<n_upper; n++) { // I think it can be parallelized
        if (n==set->core_size) op_mpi_wait_all(2,args);
        if (n==set->size) halo = 1;
        if (args[0].idx < -1) op_arg_copy_in(n,args[0], (char **)p_a[0]);
        else op_arg_set(n,args[0], &p_a[0],halo);
        if (args[1].idx < -1) op_arg_copy_in(n,args[1], (char **)p_a[1]);
        else op_arg_set(n,args[1], &p_a[1],halo);
        
        // kernel( (T0 *)p_a[0], (T1 *)p_a[1]);
        
        new_data[i] = hpx::async(kernel, (T0 *)p_a[0], (T1 *)p_a[1]); //****

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