//op_mpi_rt_support.c

void op_partition(const char* lib_name, const char* lib_routine,
                  op_set prime_set, op_map prime_map, op_dat coords )
{
    partition(lib_name, lib_routine, prime_set, prime_map, coords );
}

//////////////////////////////////////////////////////////////////////////
// op_mpi_core.c

int op_mpi_halo_exchanges(op_set set, int nargs, op_arg *args) {
    int size = set->size;
    int direct_flag = 1;
    
    if (OP_diags>0) {
        int dummy;
        for (int n=0; n<nargs; n++)
            op_arg_check(set,n,args[n],&dummy,"");
    }
    
    if (OP_hybrid_gpu) {
        for (int n=0; n<nargs; n++)
            if(args[n].opt && args[n].argtype == OP_ARG_DAT && args[n].dat->dirty_hd == 2) {
                op_download_dat(args[n].dat);
                args[n].dat->dirty_hd = 0;
            }
    }
    
    //check if this is a direct loop
    for (int n=0; n<nargs; n++)
        if(args[n].opt && args[n].argtype == OP_ARG_DAT && args[n].idx != -1)
            direct_flag = 0;
    
    if (direct_flag == 1) return size;
    
    //not a direct loop ...
    int exec_flag = 0;
    for (int n=0; n<nargs; n++) {
        if(args[n].opt && args[n].idx != -1 && args[n].acc != OP_READ) {
            size = set->size + set->exec_size;
            exec_flag = 1;
        }
    }
    op_timers_core(&c1, &t1);
    for (int n=0; n<nargs; n++) {
        if(args[n].opt && args[n].argtype == OP_ARG_DAT) {
            if (args[n].map == OP_ID) {
                op_exchange_halo(&args[n], exec_flag);///???????
            } else {
                //Check if dat-map combination was already done or if there is a mismatch (same dat, diff map)
                int found = 0;
                int fallback = 0;
                for (int m = 0; m < nargs; m++) {
                    if (m < n && args[n].dat == args[m].dat && args[n].map == args[m].map) found = 1;
                    else if (args[n].dat == args[m].dat && args[n].map != args[m].map) fallback = 1;
                }
                //If there was a map mismatch with other argument, do full halo exchange
                if (fallback) op_exchange_halo(&args[n], exec_flag);
                else if (!found) { //Otherwise, if partial halo exchange is enabled for this map, do it
                    if (OP_map_partial_exchange[args[n].map->index]) op_exchange_halo_partial(&args[n], exec_flag);
                    else op_exchange_halo(&args[n], exec_flag);
                }
            }
        }
    }
    op_timers_core(&c2, &t2);
    if (OP_kern_max>0) OP_kernels[OP_kern_curr].mpi_time += t2-t1;
    return size;
}

//////////////////////////////////////////////////////////////////////////////////
// op_mpi_rt_support.c

void op_exchange_halo(op_arg* arg, int exec_flag)
{
    op_dat dat = arg->dat;
    
    if (arg->opt == 0) return;
    
    if(arg->sent == 1)
    {
        printf("Error: Halo exchange already in flight for dat %s\n", dat->name);
        fflush(stdout);
        MPI_Abort(OP_MPI_WORLD, 2);
    }
    if (exec_flag == 0 && arg->idx == -1) return;
    
    //For a directly accessed op_dat do not do halo exchanges if not executing over
    //redundant compute block
    if (exec_flag == 0 && arg->idx == -1) return;
    
    arg->sent = 0; //reset flag
    
    //need to exchange both direct and indirect data sets if they are dirty
    if((arg->acc == OP_READ || arg->acc == OP_RW /* good for debug || arg->acc == OP_INC*/) &&
       (dat->dirtybit == 1))
    {
        printf("Exchanging Halo of data array %10s\n",dat->name); //************
        halo_list imp_exec_list = OP_import_exec_list[dat->set->index];
        halo_list imp_nonexec_list = OP_import_nonexec_list[dat->set->index];
        
        halo_list exp_exec_list = OP_export_exec_list[dat->set->index];
        halo_list exp_nonexec_list = OP_export_nonexec_list[dat->set->index];
        
        //-------first exchange exec elements related to this data array--------
        
        //sanity checks
        if(compare_sets(imp_exec_list->set,dat->set) == 0)
        {
            printf("Error: Import list and set mismatch\n");
            MPI_Abort(OP_MPI_WORLD, 2);
        }
        if(compare_sets(exp_exec_list->set,dat->set) == 0)
        {
            printf("Error: Export list and set mismatch\n");
            MPI_Abort(OP_MPI_WORLD, 2);
        }
        
        int set_elem_index;
        for(int i=0; i<exp_exec_list->ranks_size; i++) {
            for(int j = 0; j < exp_exec_list->sizes[i]; j++)
            {
                set_elem_index = exp_exec_list->list[exp_exec_list->disps[i]+j];
                memcpy(&((op_mpi_buffer)(dat->mpi_buffer))->
                       buf_exec[exp_exec_list->disps[i]*dat->size+j*dat->size],
                       (void *)&dat->data[dat->size*(set_elem_index)],dat->size);
            }
            printf("export from %d to %d data %10s, number of elements of size %d | sending:\n ", //************
                      my_rank, exp_exec_list->ranks[i], dat->name,exp_exec_list->sizes[i]);
            MPI_Isend(&((op_mpi_buffer)(dat->mpi_buffer))->
                      buf_exec[exp_exec_list->disps[i]*dat->size],
                      dat->size*exp_exec_list->sizes[i],
                      MPI_CHAR, exp_exec_list->ranks[i],
                      dat->index, OP_MPI_WORLD,
                      &((op_mpi_buffer)(dat->mpi_buffer))->
                      s_req[((op_mpi_buffer)(dat->mpi_buffer))->s_num_req++]);
        }
        
        
        int init = dat->set->size*dat->size;
        for(int i=0; i < imp_exec_list->ranks_size; i++) {
            // printf("import on to %d from %d data %10s, number of elements of size %d | recieving:\n ",
            //       my_rank, imp_exec_list->ranks[i], dat->name, imp_exec_list->sizes[i]);
            MPI_Irecv(&(dat->data[init+imp_exec_list->disps[i]*dat->size]),
                      dat->size*imp_exec_list->sizes[i],
                      MPI_CHAR, imp_exec_list->ranks[i],
                      dat->index, OP_MPI_WORLD,
                      &((op_mpi_buffer)(dat->mpi_buffer))->
                      r_req[((op_mpi_buffer)(dat->mpi_buffer))->r_num_req++]);
        }
        
        //-----second exchange nonexec elements related to this data array------
        //sanity checks
        if(compare_sets(imp_nonexec_list->set,dat->set) == 0)
        {
            printf("Error: Non-Import list and set mismatch");
            MPI_Abort(OP_MPI_WORLD, 2);
        }
        if(compare_sets(exp_nonexec_list->set,dat->set)==0)
        {
            printf("Error: Non-Export list and set mismatch");
            MPI_Abort(OP_MPI_WORLD, 2);
        }
        
        for(int i=0; i<exp_nonexec_list->ranks_size; i++) {
            for(int j = 0; j < exp_nonexec_list->sizes[i]; j++)
            {
                set_elem_index = exp_nonexec_list->list[exp_nonexec_list->disps[i]+j];
                memcpy(&((op_mpi_buffer)(dat->mpi_buffer))->
                       buf_nonexec[exp_nonexec_list->disps[i]*dat->size+j*dat->size],
                       (void *)&dat->data[dat->size*(set_elem_index)],dat->size);
            }
            //printf("export from %d to %d data %10s, number of elements of size %d | sending:\n ",
            //          my_rank, exp_nonexec_list->ranks[i], dat->name,exp_nonexec_list->sizes[i]);
            MPI_Isend(&((op_mpi_buffer)(dat->mpi_buffer))->
                      buf_nonexec[exp_nonexec_list->disps[i]*dat->size],
                      dat->size*exp_nonexec_list->sizes[i],
                      MPI_CHAR, exp_nonexec_list->ranks[i],
                      dat->index, OP_MPI_WORLD,
                      &((op_mpi_buffer)(dat->mpi_buffer))->
                      s_req[((op_mpi_buffer)(dat->mpi_buffer))->s_num_req++]);
        }
        
        int nonexec_init = (dat->set->size+imp_exec_list->size)*dat->size;
        for(int i=0; i<imp_nonexec_list->ranks_size; i++) {
            //printf("import on to %d from %d data %10s, number of elements of size %d | recieving:\n ",
            //      my_rank, imp_nonexec_list->ranks[i], dat->name, imp_nonexec_list->sizes[i]);
            MPI_Irecv(&(dat->data[nonexec_init+imp_nonexec_list->disps[i]*dat->size]),
                      dat->size*imp_nonexec_list->sizes[i],
                      MPI_CHAR, imp_nonexec_list->ranks[i],
                      dat->index, OP_MPI_WORLD,
                      &((op_mpi_buffer)(dat->mpi_buffer))->
                      r_req[((op_mpi_buffer)(dat->mpi_buffer))->r_num_req++]);
        }
        //clear dirty bit
        dat->dirtybit = 0;
        arg->sent = 1;
    }
}

void op_exchange_halo_partial(op_arg* arg, int exec_flag)
{
    op_dat dat = arg->dat;
    
    if (arg->opt == 0) return;
    
    if(arg->sent == 1)
    {
        printf("Error: Halo exchange already in flight for dat %s\n", dat->name);
        fflush(stdout);
        MPI_Abort(OP_MPI_WORLD, 2);
    }
    arg->sent = 0; //reset flag
    
    //need to exchange indirect data sets if they are dirty
    if((arg->acc == OP_READ || arg->acc == OP_RW /* good for debug || arg->acc == OP_INC*/) &&
       (dat->dirtybit == 1))
    {
        halo_list imp_nonexec_list = OP_import_nonexec_permap[arg->map->index];
        halo_list exp_nonexec_list = OP_export_nonexec_permap[arg->map->index];
        //-------exchange nonexec elements related to this data array and map--------
        
        //sanity checks
        if(compare_sets(imp_nonexec_list->set,dat->set) == 0)
        {
            printf("Error: Import list and set mismatch\n");
            MPI_Abort(OP_MPI_WORLD, 2);
        }
        if(compare_sets(exp_nonexec_list->set,dat->set) == 0)
        {
            printf("Error: Export list and set mismatch\n");
            MPI_Abort(OP_MPI_WORLD, 2);
        }
        
        int set_elem_index;
        for(int i=0; i<exp_nonexec_list->ranks_size; i++) {
            for(int j = 0; j < exp_nonexec_list->sizes[i]; j++)
            {
                set_elem_index = exp_nonexec_list->list[exp_nonexec_list->disps[i]+j];
                memcpy(&((op_mpi_buffer)(dat->mpi_buffer))->
                       buf_nonexec[exp_nonexec_list->disps[i]*dat->size+j*dat->size],
                       (void *)&dat->data[dat->size*(set_elem_index)],dat->size);
            }
            MPI_Isend(&((op_mpi_buffer)(dat->mpi_buffer))->
                      buf_nonexec[exp_nonexec_list->disps[i]*dat->size],
                      dat->size*exp_nonexec_list->sizes[i],
                      MPI_CHAR, exp_nonexec_list->ranks[i],
                      dat->index, OP_MPI_WORLD,
                      &((op_mpi_buffer)(dat->mpi_buffer))->
                      s_req[((op_mpi_buffer)(dat->mpi_buffer))->s_num_req++]);
        }
        
        int init = exp_nonexec_list->size;
        for(int i=0; i < imp_nonexec_list->ranks_size; i++) {
            MPI_Irecv(&((op_mpi_buffer)(dat->mpi_buffer))->
                      buf_nonexec[(init+imp_nonexec_list->disps[i])*dat->size],
                      dat->size*imp_nonexec_list->sizes[i],
                      MPI_CHAR, imp_nonexec_list->ranks[i],
                      dat->index, OP_MPI_WORLD,
                      &((op_mpi_buffer)(dat->mpi_buffer))->
                      r_req[((op_mpi_buffer)(dat->mpi_buffer))->r_num_req++]);
        }
        
        //note that we are not settinging the dirtybit to 0, since it's not a full exchange
        arg->sent = 1;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

