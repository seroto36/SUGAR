#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include "common.h"
#include "define.h"
#include "cosmo_tools.h"

int main(int argc, char *argv[]){

  read_argc(argc, argv);
  timer(0);
  ini_conf_values();
  read_conf_file();
  read_param_file();

  if(!strcmp(CODE_OPT, "LGC")){
    lightcone lc_cat;
    inter interv;

    init_lightcone_struct(&lc_cat);
    r_to_z();

    if(READ_LC) load_lightcone(&lc_cat);
    else{
      box *box_cat;
      init_inter_struct(&interv);
      set_redshift_intervals(&interv);
      interval_distances(&interv);
      init_boxes_struct(interv,&box_cat);
      load_boxes(interv,&box_cat);
      construct_lightcone(interv, box_cat, &lc_cat);
      free_boxes(interv.ns_in,box_cat);
      free_inter(interv);
    }
    lightcone_selection(lc_cat);
    if(RANDOM){
      ran_lc lc_ran;
      construct_random(&lc_ran,lc_cat);
      free_random(lc_ran);}
    write_output_info(interv,lc_cat.np);
    free_lightcone(lc_cat);
    end_r_to_z();
  }

  if(!strcmp(CODE_OPT, "BOX")){
    box box_cat;
    ROTBOX = 1;
    init_box_struct(&box_cat);
    load_box(&box_cat);
    box_selection(box_cat);
    free_box(box_cat);
  }
  
  return 0;
}
