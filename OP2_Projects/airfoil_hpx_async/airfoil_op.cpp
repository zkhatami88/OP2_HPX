
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/lcos/gather.hpp>
#include <hpx/include/async.hpp>
#include <hpx/lcos/when_all.hpp>
#include <hpx/runtime/serialization/serialize.hpp>

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

void op_par_loop_save_soln(char const *, op_set,
  op_arg,
  op_arg );


std::vector<hpx::future<void>> op_par_loop_adt_calc(char const *, op_set,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg );


std::vector<hpx::future<void>> op_par_loop_res_calc(char const *, op_set,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg );


std::vector<hpx::future<void>> op_par_loop_bres_calc(char const *, op_set,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg );



void op_par_loop_update(char const *, op_set,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
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

int hpx_main(int argc, char **argv)
{
  // OP initialisation
  op_init(argc,argv,2);

  int    *becell, *ecell,  *bound, *bedge, *edge, *cell;
  double  *x, *q, *qold, *adt, *res;

  int    nnode,ncell,nedge,nbedge,niter;
  double  rms;

  //timer
  double cpu_t1, cpu_t2, wall_t1, wall_t2;

  // read in grid

  op_printf("reading in grid \n");

  FILE *fp;
  if((fp=fopen("../new_grid.dat","r"))==NULL){
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

  //initialise timers for total execution wall time
  op_timers(&cpu_t1, &wall_t1);

  // main time-marching loop

  niter = 500;
   
  boost::uint64_t t=hpx::util::high_resolution_clock::now();

    



  for(int iter=1; iter<=niter; iter++) {

  std::vector<std::vector<hpx::future<void>>> new_data1B;

    // save old flow solution

    op_par_loop_save_soln("save_soln",cells,
                op_arg_dat(p_q,-1,OP_ID,4,"double",OP_READ),
                op_arg_dat(p_qold,-1,OP_ID,4,"double",OP_WRITE));

    // predictor/corrector update loop

   for(int k=0; k<2; k++) {

      // calculate area/timstep
     new_data1B.push_back(op_par_loop_adt_calc("adt_calc",cells,
                  op_arg_dat(p_x,0,pcell,2,"double",OP_READ),
                  op_arg_dat(p_x,1,pcell,2,"double",OP_READ),
                  op_arg_dat(p_x,2,pcell,2,"double",OP_READ),
                  op_arg_dat(p_x,3,pcell,2,"double",OP_READ),
                  op_arg_dat(p_q,-1,OP_ID,4,"double",OP_READ),
                  op_arg_dat(p_adt,-1,OP_ID,1,"double",OP_WRITE)));

      // calculate flux residual
     new_data1B.push_back(op_par_loop_res_calc("res_calc",edges,
                  op_arg_dat(p_x,0,pedge,2,"double",OP_READ),
                  op_arg_dat(p_x,1,pedge,2,"double",OP_READ),
                  op_arg_dat(p_q,0,pecell,4,"double",OP_READ),
                  op_arg_dat(p_q,1,pecell,4,"double",OP_READ),
                  op_arg_dat(p_adt,0,pecell,1,"double",OP_READ),
                  op_arg_dat(p_adt,1,pecell,1,"double",OP_READ), 
                  op_arg_dat(p_res,0,pecell,4,"double",OP_INC),  
                  op_arg_dat(p_res,1,pecell,4,"double",OP_INC)));  

      new_data1B.push_back(op_par_loop_bres_calc("bres_calc",bedges,
                  op_arg_dat(p_x,0,pbedge,2,"double",OP_READ),
                  op_arg_dat(p_x,1,pbedge,2,"double",OP_READ),
                  op_arg_dat(p_q,0,pbecell,4,"double",OP_READ),
                  op_arg_dat(p_adt,0,pbecell,1,"double",OP_READ),
                  op_arg_dat(p_res,0,pbecell,4,"double",OP_INC),
                  op_arg_dat(p_bound,-1,OP_ID,1,"int",OP_READ)));


	for(int i=0; i<new_data1B.size();++i)
			for(int j=0; j<new_data1B[i].size();++j)
				new_data1B[i][j].wait();


      // update flow field
      rms = 0.0;
     op_par_loop_update("update",cells,
                  	op_arg_dat(p_qold,-1,OP_ID,4,"double",OP_READ),
                  	op_arg_dat(p_q,-1,OP_ID,4,"double",OP_WRITE),
                  	op_arg_dat(p_res,-1,OP_ID,4,"double",OP_RW),
                  	op_arg_dat(p_adt,-1,OP_ID,1,"double",OP_READ),
                  	op_arg_gbl(&rms,1,"double",OP_INC));

    }

    // print iteration history
    rms = sqrt(rms/(double) op_get_size(cells));
    if (iter%100 == 0)
      op_printf(" %d  %10.5e \n",iter,rms);
  }


  boost::uint64_t elapsed=hpx::util::high_resolution_clock::now()-t;
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

	return hpx::finalize();
}


int main(){
	return hpx::init();
}




