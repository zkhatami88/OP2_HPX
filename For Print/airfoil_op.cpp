

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/lcos/gather.hpp>
#include <hpx/runtime/serialization/vector.hpp>
#include <hpx/runtime/serialization/serialize.hpp>
#include <hpx/include/async.hpp>
// global constants

double gam, gm1, cfl, eps, mach, alpha, qinf[4];

//
// OP header file
//

#include  "op_lib_cpp.h"


#define STRIDE(x,y) x
int nodes_stride = 1;
int edges_stride = 1;
int bedges_stride = 1;
int cells_stride = 1;
//
// op_par_loop declarations
//



hpx::future<op_arg>
op_arg_dat1 ( hpx::shared_future<op_dat> dat, int idx, op_map map, int dim, const char * typ, op_access acc )
{
    
    using hpx::lcos::local::dataflow;
    using hpx::util::unwrapped;
  
    return dataflow(unwrapped([idx,map,dim,&typ,acc](op_dat dat){
    op_arg arg;
    
    
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
        
        arg.size = -1;
        arg.data = NULL;
        arg.data_d = NULL;
        arg.map_data_d = NULL;
        arg.map_data = NULL;
    }
    
    arg.type = typ;
    arg.acc = acc;
    
    
    arg.sent = 0;
    
    return arg;

    }),dat);

}

hpx::shared_future<op_dat> op_par_loop_save_soln(char const *, op_set,
  hpx::future<op_arg>,
  hpx::future<op_arg> );

hpx::shared_future<op_dat> op_par_loop_adt_calc(char const *, op_set,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg> );

hpx::shared_future<op_dat> op_par_loop_res_calc(char const *, op_set,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg> );

hpx::shared_future<op_dat> op_par_loop_bres_calc(char const *, op_set,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg> );

hpx::shared_future<op_dat> op_par_loop_update1(char const *, op_set,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  op_arg );

hpx::shared_future<op_dat> op_par_loop_update2(char const *, op_set,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  hpx::future<op_arg>,
  op_arg );


//
// kernel routines for parallel loops
//

#include "save_soln.h"
#include "adt_calc.h"
#include "res_calc.h"
#include "bres_calc.h"
#include "update.h"

// main program

int hpx_main(int argc, char** argv){



  // OP initialisation
  op_init(argc,argv,2);

  int    *becell, *ecell,  *bound, *bedge, *edge, *cell;
  double  *x, *q, *qold, *adt, *res;

  int    nnode,ncell,nedge,nbedge,niter;
  double  rms;

  //timer
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  hpx::future<void> myvar;
  // read in grid

  op_printf("reading in grid \n");

  FILE *fp;
  //if ( (fp = fopen("./new_grid.dat","r")) == NULL) {  
  if((fp=fopen("../../FirstHPX/new_grid.dat","r"))==NULL){
    op_printf("can't open file new_grid.dat\n"); exit(-1);
  }

  if (fscanf(fp,"%d %d %d %d \n",&nnode, &ncell, &nedge, &nbedge) != 4) {
    op_printf("error reading from new_grid.dat\n"); exit(-1);
  }

  cell   = (int *) malloc(4*ncell*sizeof(int));
  edge   = (int *) malloc(2*nedge*sizeof(int));
  ecell  = (int *) malloc(2*nedge*sizeof(int));
  bedge  = (int *) malloc(2*nbedge*sizeof(int));
  becell = (int *) malloc(  nbedge*sizeof(int));
  bound  = (int *) malloc(  nbedge*sizeof(int));

  x      = (double *) malloc(2*nnode*sizeof(double));
  q      = (double *) malloc(4*ncell*sizeof(double));
  qold   = (double *) malloc(4*ncell*sizeof(double));
  res    = (double *) malloc(4*ncell*sizeof(double));
  adt    = (double *) malloc(  ncell*sizeof(double));

  for (int n=0; n<nnode; n++) {
    if (fscanf(fp,"%lf %lf \n",&x[2*n], &x[2*n+1]) != 2) {
      op_printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  for (int n=0; n<ncell; n++) {
    if (fscanf(fp,"%d %d %d %d \n",&cell[4*n  ], &cell[4*n+1],
                                   &cell[4*n+2], &cell[4*n+3]) != 4) {
      op_printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  for (int n=0; n<nedge; n++) {
    if (fscanf(fp,"%d %d %d %d \n",&edge[2*n], &edge[2*n+1],
                                   &ecell[2*n],&ecell[2*n+1]) != 4) {
      op_printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  for (int n=0; n<nbedge; n++) {
    if (fscanf(fp,"%d %d %d %d \n",&bedge[2*n],&bedge[2*n+1],
                                   &becell[n], &bound[n]) != 4) {
      op_printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  fclose(fp);

  // set constants and initialise flow field and residual

  op_printf("initialising flow field \n");

  gam = 1.4f;
  gm1 = gam - 1.0f;
  cfl = 0.9f;
  eps = 0.05f;

  double mach  = 0.4f;
  double alpha = 3.0f*atan(1.0f)/45.0f;
  double p     = 1.0f;
  double r     = 1.0f;
  double u     = sqrt(gam*p/r)*mach;
  double e     = p/(r*gm1) + 0.5f*u*u;

  qinf[0] = r;
  qinf[1] = r*u;
  qinf[2] = 0.0f;
  qinf[3] = r*e;

  for (int n=0; n<ncell; n++) {
    for (int m=0; m<4; m++) {
        q[4*n+m] = qinf[m];
      res[4*n+m] = 0.0f;
    }
  }

  // declare sets, pointers, datasets and global constants

  op_set nodes  = op_decl_set(nnode,  "nodes");
  op_set edges  = op_decl_set(nedge,  "edges");
  op_set bedges = op_decl_set(nbedge, "bedges");
  op_set cells  = op_decl_set(ncell,  "cells");

  op_map pedge   = op_decl_map(edges, nodes,2,edge,  "pedge");
  op_map pecell  = op_decl_map(edges, cells,2,ecell, "pecell");
  op_map pbedge  = op_decl_map(bedges,nodes,2,bedge, "pbedge");
  op_map pbecell = op_decl_map(bedges,cells,1,becell,"pbecell");
  op_map pcell   = op_decl_map(cells, nodes,4,cell,  "pcell");

  op_dat p_bound = op_decl_dat(bedges,1,"int"  ,bound,"p_bound");
  op_dat p_x     = op_decl_dat(nodes ,2,"double",x    ,"p_x");
  op_dat p_q     = op_decl_dat(cells ,4,"double",q    ,"p_q");
  op_dat p_qold  = op_decl_dat(cells ,4,"double",qold ,"p_qold");
  op_dat p_adt   = op_decl_dat(cells ,1,"double",adt  ,"p_adt");
  op_dat p_res   = op_decl_dat(cells ,4,"double",res  ,"p_res");

  op_decl_const2("gam",1,"double",&gam);
  op_decl_const2("gm1",1,"double",&gm1);
  op_decl_const2("cfl",1,"double",&cfl);
  op_decl_const2("eps",1,"double",&eps);
  op_decl_const2("mach",1,"double",&mach);
  op_decl_const2("alpha",1,"double",&alpha);
  op_decl_const2("qinf",4,"double",qinf);

  op_diagnostic_output();

  op_timers(&cpu_t1, &wall_t1);

  niter = 1000;

    hpx::shared_future<op_dat> dat0=hpx::make_ready_future(p_q);
    hpx::shared_future<op_dat> dat1=hpx::make_ready_future(p_qold);
    hpx::shared_future<op_dat> dat2=hpx::make_ready_future(p_adt);
    hpx::shared_future<op_dat> dat3=hpx::make_ready_future(p_res);
    hpx::shared_future<op_dat> dat4=hpx::make_ready_future(p_x);
    hpx::shared_future<op_dat> dat5=hpx::make_ready_future(p_bound);

    std::vector<hpx::shared_future<op_dat>> v0;
    std::vector<hpx::shared_future<op_dat>> v1;
    std::vector<hpx::shared_future<op_dat>> v2;
    std::vector<hpx::shared_future<op_dat>> v3;
    std::vector<hpx::shared_future<op_dat>> v4;
    std::vector<hpx::shared_future<op_dat>> v5;
    
    v0.resize(niter);
    v1.resize(niter);
    v2.resize(niter);
    v3.resize(niter);
    v4.resize(niter);
    v5.resize(niter);

  for(int i=0; i<niter; ++i){    
    v0[i]=dat0;
    v1[i]=dat1;
    v2[i]=dat2;
    v3[i]=dat3;
    v4[i]=dat4;
    v5[i]=dat5;
}

boost::uint64_t t = hpx::util::high_resolution_clock::now();
   
  for(int iter=1; iter<niter-1; iter++) {

    v1[iter]=op_par_loop_save_soln("save_soln", cells, 
      op_arg_dat1(v0[iter-1], -1,OP_ID, 4,"double",OP_READ ),
      op_arg_dat1(v1[iter-1], -1,OP_ID, 4,"double",OP_WRITE));

      // calculate area/timstep

      v2[iter]=op_par_loop_adt_calc("adt_calc",cells, 
          op_arg_dat1(v4[0],   0,pcell, 2,"double",OP_READ ),
          op_arg_dat1(v4[0],   1,pcell, 2,"double",OP_READ ),
          op_arg_dat1(v4[0],   2,pcell, 2,"double",OP_READ ),
          op_arg_dat1(v4[0],   3,pcell, 2,"double",OP_READ ),
          op_arg_dat1(v0[iter-1],  -1,OP_ID, 4,"double",OP_READ ),
          op_arg_dat1(v2[iter-1],-1,OP_ID, 1,"double",OP_WRITE));

      v3[iter]=op_par_loop_res_calc("res_calc",edges, 
          op_arg_dat1(v4[0],    0,pedge, 2,"double",OP_READ),
          op_arg_dat1(v4[0],    1,pedge, 2,"double",OP_READ),
          op_arg_dat1(v0[iter-1],    0,pecell,4,"double",OP_READ),
          op_arg_dat1(v0[iter-1],    1,pecell,4,"double",OP_READ),
          op_arg_dat1(v2[iter],  0,pecell,1,"double",OP_READ),
          op_arg_dat1(v2[iter],  1,pecell,1,"double",OP_READ),
          op_arg_dat1(v3[iter-1],  0,pecell,4,"double",OP_INC ),
          op_arg_dat1(v3[iter-1],  1,pecell,4,"double",OP_INC ));



      v3[iter]=op_par_loop_bres_calc("bres_calc",bedges,
          op_arg_dat1(v4[0],     0,pbedge, 2,"double",OP_READ),
          op_arg_dat1(v4[0],     1,pbedge, 2,"double",OP_READ),
          op_arg_dat1(v0[iter-1],     0,pbecell,4,"double",OP_READ),
          op_arg_dat1(v2[iter],   0,pbecell,1,"double",OP_READ),
          op_arg_dat1(v3[iter],   0,pbecell,4,"double",OP_INC ), 
          op_arg_dat1(v5[0],-1,OP_ID  ,1,"int",  OP_READ));


     rms = 0.0;

      v3[iter]=op_par_loop_update1("update1",cells,
          op_arg_dat1(v1[iter],-1,OP_ID, 4,"double",OP_READ ),
          op_arg_dat1(v0[iter-1],   -1,OP_ID, 4,"double",OP_WRITE),
          op_arg_dat1(v3[iter], -1,OP_ID, 4,"double",OP_RW   ),
          op_arg_dat1(v2[iter], -1,OP_ID, 1,"double",OP_READ ), 
          op_arg_gbl(&rms,1,"double",OP_INC));

        
     v0[iter]=op_par_loop_update2("update2",cells,
          op_arg_dat1(v1[iter],-1,OP_ID, 4,"double",OP_READ ),
          op_arg_dat1(v0[iter-1],   -1,OP_ID, 4,"double",OP_WRITE),
          op_arg_dat1(v3[iter], -1,OP_ID, 4,"double",OP_RW   ),
          op_arg_dat1(v2[iter], -1,OP_ID, 1,"double",OP_READ ), 
          op_arg_gbl(&rms,1,"double",OP_INC));

  }


when_all(v0).get();
when_all(v1).get();
when_all(v2).get();
when_all(v3).get();
when_all(v4).get();
when_all(v5).get();

boost::uint64_t elapsed = hpx::util::high_resolution_clock::now() - t;
std::cout<<(boost::format("%.14g")%(elapsed/ 1e9)) <<std::flush<<std::endl;

  op_timers(&cpu_t2, &wall_t2);

  //output the result dat array to files
  op_print_dat_to_txtfile(p_q, "out_grid_seq.dat"); //ASCI
  op_print_dat_to_binfile(p_q, "out_grid_seq.bin"); //Binary

  op_timing_output();
  op_printf("Max total runtime = \n%f\n",wall_t2-wall_t1);

  op_exit();

  free(cell);
  free(edge);
  free(ecell);
  free(bedge);
  free(becell);
  free(bound);
  free(x);
  free(q);
  free(qold);
  free(res);
  free(adt);

  std::cout<<"Yayyyy"<<std::endl;
  std::cout<<"----------------- \n";

  hpx::finalize();
  return 0;

};

int main(int argc, char** argv){
    return hpx::init(argc, argv);
}
