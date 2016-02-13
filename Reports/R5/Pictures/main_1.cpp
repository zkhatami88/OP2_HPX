

    op_par_loop_save_soln("save_soln",cells,
                op_arg_dat(p_q,-1,OP_ID,4,"double",OP_READ),
                op_arg_dat(p_qold,-1,OP_ID,4,"double",OP_WRITE));

    for(int k=0; k<2; k++) {

      op_par_loop_adt_calc("adt_calc",cells,
                  op_arg_dat(p_x,0,pcell,2,"double",OP_READ),
                  op_arg_dat(p_x,1,pcell,2,"double",OP_READ),
                  op_arg_dat(p_x,2,pcell,2,"double",OP_READ),
                  op_arg_dat(p_x,3,pcell,2,"double",OP_READ),
                  op_arg_dat(p_q,-1,OP_ID,4,"double",OP_READ),
                  op_arg_dat(p_adt,-1,OP_ID,1,"double",OP_WRITE));

      op_par_loop_res_calc("res_calc",edges,
                  op_arg_dat(p_x,0,pedge,2,"double",OP_READ),
                  op_arg_dat(p_x,1,pedge,2,"double",OP_READ),
                  op_arg_dat(p_q,0,pecell,4,"double",OP_READ),
                  op_arg_dat(p_q,1,pecell,4,"double",OP_READ),
                  op_arg_dat(p_adt,0,pecell,1,"double",OP_READ),
                  op_arg_dat(p_adt,1,pecell,1,"double",OP_READ),
                  op_arg_dat(p_res,0,pecell,4,"double",OP_INC),
                  op_arg_dat(p_res,1,pecell,4,"double",OP_INC));

      op_par_loop_bres_calc("bres_calc",bedges,
                  op_arg_dat(p_x,0,pbedge,2,"double",OP_READ),
                  op_arg_dat(p_x,1,pbedge,2,"double",OP_READ),
                  op_arg_dat(p_q,0,pbecell,4,"double",OP_READ),
                  op_arg_dat(p_adt,0,pbecell,1,"double",OP_READ),
                  op_arg_dat(p_res,0,pbecell,4,"double",OP_INC),
                  op_arg_dat(p_bound,-1,OP_ID,1,"int",OP_READ));

      op_par_loop_update("update",cells,
                  op_arg_dat(p_qold,-1,OP_ID,4,"double",OP_READ),
                  op_arg_dat(p_q,-1,OP_ID,4,"double",OP_WRITE),
                  op_arg_dat(p_res,-1,OP_ID,4,"double",OP_RW),
                  op_arg_dat(p_adt,-1,OP_ID,1,"double",OP_READ),
                  op_arg_gbl(&rms,1,"double",OP_INC));

    }



