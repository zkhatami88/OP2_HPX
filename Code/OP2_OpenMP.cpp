
#include <op_lib_c.h>
#include <op_rt_support.h>

/*
 * Routines called by user code and kernels
 * these wrappers are used by non-CUDA versions
 * op_lib.cu provides wrappers for CUDA version
 */

void
op_init ( int argc, char ** argv, int diags )
{
    op_init_core ( argc, argv, diags );
}

op_dat
op_decl_dat_char ( op_set set, int dim, char const * type, int size, char * data, char const *name )
{
    return op_decl_dat_core ( set, dim, type, size, data, name );
}

op_dat
op_decl_dat_temp_char (op_set set, int dim, char const * type, int size, char const *name )
{
    char* data = NULL;
    op_dat dat = op_decl_dat_temp_core ( set, dim, type, size, data, name );
    
    dat->data = (char*) calloc(set->size*dim*size, 1); //initialize data bits to 0
    dat-> user_managed = 0;
    return dat;
}

int op_free_dat_temp_char ( op_dat dat )
{
    return op_free_dat_temp_core (dat);
}

void op_upload_all ()
{
}

void
op_fetch_data_char ( op_dat dat, char * usr_ptr )
{
    //need to copy data into memory pointed to by usr_ptr
    memcpy((void *)usr_ptr, (void *)dat->data, dat->set->size*dat->size);
}

void
op_fetch_data_hdf5_char ( op_dat dat, char * usr_ptr, int low, int high )
{
    if(low < 0 || high > dat->set->size -1)
    {
        printf("op_fetch_data_hdf5: Indices not within range of elements held in %s\n",
               dat->name);
        exit(2);
    }
    //need to copy data into memory pointed to by usr_ptr
    memcpy((void *)usr_ptr, (void *)&dat->data[low*dat->size],
           (high+1)*dat->size);
}

/*
 * No specific action is required for constants in OpenMP
 */

void
op_decl_const_char ( int dim, char const * type, int typeSize, char * data, char const * name )
{
    (void)dim;
    (void)type;
    (void)typeSize;
    (void)data;
    (void)name;
}

op_plan *
op_plan_get ( char const * name, op_set set, int part_size,
             int nargs, op_arg * args, int ninds, int *inds )
{
    return op_plan_get_stage ( name, set, part_size, nargs, args, ninds, inds, OP_STAGE_ALL );
}

op_plan *
op_plan_get_stage ( char const * name, op_set set, int part_size,
                   int nargs, op_arg * args, int ninds, int *inds, int staging )
{
    return op_plan_core ( name, set, part_size, nargs, args, ninds, inds, staging );
}


int op_get_size(op_set set)
{
    return set->size;
}

void op_printf(const char* format, ...)
{
    va_list argptr;
    va_start(argptr, format);
    vprintf(format, argptr);
    va_end(argptr);
}

void op_print(const char* line)
{
    printf("%s\n",line);
}

void
op_exit (  )
{
    op_rt_exit (  );
    op_exit_core (  );
}

/*
 * Wrappers of core lib
 */

op_set
op_decl_set ( int size, char const *name )
{
    return op_decl_set_core ( size, name );
}

op_map
op_decl_map ( op_set from, op_set to, int dim, int * imap, char const * name )
{
    return op_decl_map_core ( from, to, dim, imap, name );
}

op_arg
op_arg_dat ( op_dat dat, int idx, op_map map, int dim, char const * type, op_access acc )
{
    return op_arg_dat_core ( dat, idx, map, dim, type, acc );
}

op_arg
op_opt_arg_dat ( int opt, op_dat dat, int idx, op_map map, int dim, char const * type, op_access acc )
{
    return op_opt_arg_dat_core ( opt, dat, idx, map, dim, type, acc );
}

void op_timers(double * cpu, double * et)
{
    op_timers_core(cpu,et);
}

op_arg
op_arg_gbl_char ( char * data, int dim, const char *type, int size, op_access acc )
{
    return op_arg_gbl_core ( data, dim, type, size, acc );
}

void op_timing_output()
{
    op_timing_output_core();
    printf("Total plan time: %8.4f\n", OP_plan_time);
}

void op_print_dat_to_binfile(op_dat dat, const char *file_name)
{
    op_print_dat_to_binfile_core(dat, file_name);
}

void op_print_dat_to_txtfile(op_dat dat, const char *file_name)
{
    op_print_dat_to_txtfile_core(dat, file_name);
}