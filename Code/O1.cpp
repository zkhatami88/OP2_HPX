


double * bres_calc(const double *x1, const double *x2, const double *q1,
                      const double *adt1, double *res1, const int *bound) {
    
    double dx,dy,mu, ri, p1,vol1, p2,vol2, f;
    
    dx = x1[0] - x2[0];
    dy = x1[1] - x2[1];
    
    ri = 1.0f/q1[0];
    p1 = gm1*(q1[3]-0.5f*ri*(q1[1]*q1[1]+q1[2]*q1[2]));
    
    if (*bound==1) {
        res1[1] += + p1*dy;
        res1[2] += - p1*dx;
    }
    else {
        vol1 =  ri*(q1[1]*dy - q1[2]*dx);
        
        ri   = 1.0f/qinf[0];
        p2   = gm1*(qinf[3]-0.5f*ri*(qinf[1]*qinf[1]+qinf[2]*qinf[2]));
        vol2 =  ri*(qinf[1]*dy - qinf[2]*dx);
        
        mu = (*adt1)*eps;
        
        f = 0.5f*(vol1* q1[0]         + vol2* qinf[0]        ) + mu*(q1[0]-qinf[0]);
        res1[0] += f;
        f = 0.5f*(vol1* q1[1] + p1*dy + vol2* qinf[1] + p2*dy) + mu*(q1[1]-qinf[1]);
        res1[1] += f;
        f = 0.5f*(vol1* q1[2] - p1*dx + vol2* qinf[2] - p2*dx) + mu*(q1[2]-qinf[2]);
        res1[2] += f;
        f = 0.5f*(vol1*(q1[3]+p1)     + vol2*(qinf[3]+p2)    ) + mu*(q1[3]-qinf[3]);
        res1[3] += f;
    }
    
    return res1;
}




for(int iter=1; iter<=niter; iter++)
{
    
    op_plan *Plan = op_plan_get(name,set,part_size,nargs,args,ninds,inds);
    
    for ( int blockIdx=0; blockIdx<nblocks; blockIdx++ ){
    
    int blockId  = Plan->blkmap[blockIdx + block_offset];
    int nelem    = Plan->nelems[blockId];
    int offset_b = Plan->offset[blockId];
    
    for ( int n=offset_b; n<offset_b+nelem; n++ ){
        
        int map0idx = arg0.map_data[n * arg0.map->dim + 0];
        int map1idx = arg0.map_data[n * arg0.map->dim + 1];
        int map2idx = arg2.map_data[n * arg2.map->dim + 0];
        
        next[i] = dataflow(
                           hpx::launch::async, Op,
                           current[idx(i-1, np)], current[i], current[idx(i+1, np)]
                           );
        
        arg4.data = dataflow(
                             hpx::launch::async, bres_calc,
                             current[idx(i-1, np)], current[i], current[idx(i+1, np)]
                             );
        

        bres_calc(
                  &((hpx::future<double*>)arg0.data)[2 * map0idx],
                  &((hpx::future<double*>)arg0.data)[2 * map1idx],
                  &((hpx::future<double*>)arg2.data)[4 * map2idx],
                  &((hpx::future<double*>)arg3.data)[1 * map2idx],
                  &((hpx::future<double*>)arg4.data)[4 * map2idx],
                  &((hpx::future<int*>)arg5.data)[1 * n]);
        
        
        /*arg4.data = dataflow(
                        hpx::launch::async,
                        unwrapped(
                                  [(hpx::future<double*>)arg0.data)[2 * map0idx],
                                   (hpx::future<double*>)arg0.data)[2 * map1idx],
                                   (hpx::future<double*>)arg2.data)[4 * map2idx],
                                   (hpx::future<double*>)arg3.data)[1 * map2idx],
                                   (hpx::future<double*>)arg4.data)[4 * map2idx],
                                   (hpx::future<int*>)arg5.data)[1 * n])]
                                  (double *data0, double *data1, double *data2,
                                   double *data3, double *data4, int *data5) -> hpx::future<double*>
                                  {
                                  
                                  
                                  
                                  }
                                  ),
                             
                             arg0.data[2 * map0idx].get,
                             arg0.data[2 * map1idx].get,
                             arg2.data[4 * map2idx].get,
                             arg3.data[1 * map2idx].get,
                             arg4.data[4 * map2idx].get,
                             arg5.data)[1 * n].get
        );*/
        
        hpx::future<double*> data0 = arg0.data[2 * map0idx];
        hpx::future<double*> data1 = arg0.data[2 * map1idx];
        hpx::future<double*> data2 = arg2.data[4 * map2idx];
        hpx::future<double*> data3 = arg3.data[1 * map2idx];
        hpx::future<double*> data4 = arg4.data)[4 * map2idx];
        hpx::future<int*> data5 = arg5.data[1 * n];

        arg4.data = dataflow(unwrapped(bres_calc), (hpx::future<double*>arg0.data)[2 * map0idx],
                             (hpx::future<double*>arg0.data)[2 * map1idx],
                             (hpx::future<double*>arg2.data)[4 * map2idx],
                             (hpx::future<double*>arg3.data)[1 * map2idx],
                             (hpx::future<double*>arg4.data)[4 * map2idx],
                             (hpx::future<int*>arg5.data)[1 * n]);
    }

}




for ( int blockIdx=0; blockIdx<nblocks; blockIdx++ ){
    int blockId  = Plan->blkmap[blockIdx + block_offset];
    int nelem    = Plan->nelems[blockId];
    int offset_b = Plan->offset[blockId];
    for ( int n=offset_b; n<offset_b+nelem; n++ ){
        int map0idx = arg0.map_data[n * arg0.map->dim + 0];
        int map1idx = arg0.map_data[n * arg0.map->dim + 1];
        int map2idx = arg2.map_data[n * arg2.map->dim + 0];
        
        bres_calc(
                  &((double*)arg0.data)[2 * map0idx],
                  &((double*)arg0.data)[2 * map1idx],
                  &((double*)arg2.data)[4 * map2idx],
                  &((double*)arg3.data)[1 * map2idx],
                  &((double*)arg4.data)[4 * map2idx],
                  &((int*)arg5.data)[1 * n]);
    }
}
    
    
    partition_data boost::shared_array<double>
    
    typedef hpx::shared_future<partition_data> partition;
    typedef std::vector<partition> space;
    
    std::vector<hpx::shared_future<double>> space;
    or
    std::vector<hpx::shared_future<boost::shared_array<double>>> space; //vector is for each time & in each time we have future<array<double>>
    

    static partition_data heat_part(partition_data const& left,
                                    partition_data const& middle, partition_data const& right)
    {
        std::size_t size = middle.size();
        partition_data next(size);
        
        typedef boost::counting_iterator<std::size_t> iterator;
        
        next[0] = heat(left[size-1], middle[0], middle[1]);
        
        using namespace hpx::parallel;
        for_each(par, iterator(1), iterator(size-1),
                 [&next, &middle](std::size_t i)
                 {
                     next[i] = heat(middle[i-1], middle[i], middle[i+1]);
                 });
        
        next[size-1] = heat(middle[size-2], middle[size-1], right[0]);
        
        return next;
    }
    
    
    for (std::size_t t = 0; t != nt; ++t)
    {
        space const& current = U[t % 2];
        space& next = U[(t + 1) % 2];
        
        typedef boost::counting_iterator<std::size_t> iterator;
        
        for_each(par, iterator(0), iterator(np),
                 [&next, &current, np, &Op](std::size_t i)
                 {
                     next[i] = dataflow(
                                        hpx::launch::async, Op,
                                        current[idx(i-1, np)], current[i], current[idx(i+1, np)]
                                        );
                 });
    }




