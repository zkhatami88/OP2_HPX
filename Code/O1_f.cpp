

std::vector<double> bres_calc(const double *x1, const double *x2, const double *q1,
                      const double *adt1, double *res1, const int *bound) {
    
    double dx,dy,mu, ri, p1,vol1, p2,vol2, f;
    std::vector<double> result(3,0);
    
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
    
    result=res1;
    
    return result;

}

//should be produced
std::vector<std::vector<double>> work(int start, int fin, op_arg data0, op_arg data1, op_arg data2, op_arg data3,
                         op_arg data4, op_arg data5){

    typedef boost::counting_iterator<std::size_t> iterator;
    std::vector<std::vector<double>> new_data(nelem,data4.data->size);
    
    using namespace hpx::parallel;
    for_each(par, iterator(start), iterator(finish-1),
              [&&new_data, &data0, &data1, &data2, &data3, &data4, &data5](std::size_t i)
             {
                 int map0idx = data0.map_data[i * arg0.map->dim + 0];
                 int map1idx = data0.map_data[i * arg0.map->dim + 1];
                 int map2idx = data2.map_data[i * arg2.map->dim + 0];
    
                 new_data[i]=bres_calc(
                           &((double*)data0.data)[2 * map0idx],
                           &((double*)data0.data)[2 * map1idx],
                           &((double*)data2.data)[4 * map2idx],
                           &((double*)data3.data)[1 * map2idx],
                           &((double*)data4.data)[4 * map2idx],
                           &((int*)data5.data)[1 * n]);
             });
    
    return new_data; //how to return two ar.data??
}

//should be produced
int main_hpx(){
    
    op_plan *Plan = op_plan_get(name,set,part_size,nargs,args,ninds,inds);
    
    hpx::future<op_arg> data0 = hpx::make_ready_future(arg0);
    hpx::future<op_arg> data1 = hpx::make_ready_future(arg1);
    hpx::future<op_arg> data2 = hpx::make_ready_future(arg2);
    hpx::future<op_arg> data3 = hpx::make_ready_future(arg3);
    hpx::future<op_arg> data4 = hpx::make_ready_future(arg4);
    hpx::future<op_arg> data5 = hpx::make_ready_future(arg5);
    
    typedef boost::counting_iterator<std::size_t> blockIdx;
    std::vector<std::vector<hpx::future<double>>> new_data(max_nelem,arg4.data->size); //max_nelem should be computed before following loop then the size of each should be changed

    for_each(par, blockIdx(0), blockIdx(nblocks),
                [&Plan, &unwrapped(bres_calc), &new_data,
                &data0, &data1, &data2, &data3, &data4, &data5](std::size_t i)
                {
                    int blockId  = Plan->blkmap[i + block_offset];
                    int nelem    = Plan->nelems[i];
                    int offset_b = Plan->offset[i];
                    
                    new_data = dataflow(
                                        hpx::launch::async, unwrapped(work),
                                        hpx::make_ready_future(offset_b),
                                        hpx::make_ready_future(offset_b+nelem),
                                        data0, data1, data2, data3, data4, data5); //should be future
                });
    
    
    return hpx::when_all(new_data); // which new_data is arg4.data

}


