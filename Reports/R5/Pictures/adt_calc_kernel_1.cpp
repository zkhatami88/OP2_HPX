

      #pragma omp parallel for
      for ( int blockIdx=0; blockIdx<nblocks; blockIdx++ ){
        int blockId  = Plan->blkmap[blockIdx + block_offset];
        int nelem    = Plan->nelems[blockId];
        int offset_b = Plan->offset[blockId];
        for ( int n=offset_b; n<offset_b+nelem; n++ ){
          int map0idx = arg0.map_data[n * arg0.map->dim + 0];
          int map1idx = arg0.map_data[n * arg0.map->dim + 1];
          int map2idx = arg0.map_data[n * arg0.map->dim + 2];
          int map3idx = arg0.map_data[n * arg0.map->dim + 3];

          adt_calc(
            &((double*)arg0.data)[2 * map0idx],
            &((double*)arg0.data)[2 * map1idx],
            &((double*)arg0.data)[2 * map2idx],
            &((double*)arg0.data)[2 * map3idx],
            &((double*)arg4.data)[4 * n],
            &((double*)arg5.data)[1 * n]);
        }
      }




