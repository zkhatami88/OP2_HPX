
#348
    code('int blockId = blkmap[blockIdx + block_offset];')
    code('nelem    = nelems[blockId];')
    code('offset_b = offset[blockId];')
    code('')
    code('new_data[blockIdx].resize(nelem);')
    code('typedef boost::counting_iterator<std::size_t> iterator;')
    #std::vector<std::vector<hpx::shared_future<double>> new_data(nelem,
        #std::vector<hpx::shared_future<double>(arg4.data->size));
    code('using namespace hpx::parallel;')
    
#644 kernel call for indirect version
    #
    if ninds>0:
        code('op_plan *Plan = op_plan_get(name,set,part_size,nargs,args,ninds,inds);')
        code('')
        comm(' execute plan')
        code('int block_offset = 0;')
        FOR('col','0','Plan->ncolors')
        IF('col==Plan->ncolors_core')
        code('op_mpi_wait_all(nargs, args);')
        ENDIF()
        code('int nblocks = Plan->ncolblk[col];')
        code('')
        
        comm(' hpx parallel')
        comm(' // new_data[i][0 ,1, 2] which i is blockId')
        code('std::vector<std::vector<hpx::shared_future<double> > > new_data(2);')
        FOR('blockIdx','0','nblocks')
        code('op_x86_'+name+'( blockIdx,')
            
        for m in range(1,ninds+1):
            g_m = invinds[m-1]
            code('(TYP *)ARG.data,')
                                                                        
        code('Plan->ind_map,')
        code('Plan->loc_map,')