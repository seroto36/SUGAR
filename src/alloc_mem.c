#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "cosmo_tools.h"
#include "common.h"
#include "define.h"


long n_halos_ascii(char *file){

   FILE     *fd;
   char     buffer[CHAR_SIZE];
   long     lines=0;

   if(!(fd = fopen(file, "r")))
     error_open_file(file);
   
   while(fgets(buffer,CHAR_SIZE,fd) != NULL)lines++;
   
   fclose(fd);

   return(lines);
}

void init_lightcone_struct(lightcone *light_cat){
  (*light_cat).np = 0;
  (*light_cat).lc = NULL;
}

void init_lightcone_rand(ran_lc *light_cat){
  (*light_cat).np = 0;
  (*light_cat).lc = NULL;
}

void init_box_struct(box *box_cat){
    (*box_cat).np = 0;
    (*box_cat).bx = NULL;
}

void init_boxes_struct(inter interv,box **box_cat){
  int ii=0;
  (*box_cat) = NULL;
  (*box_cat) = (box *)malloc(interv.ns_in*sizeof(box));
  if((*box_cat)==NULL) error_malloc_allocation("box_cat");
  for(ii=0; ii<interv.ns_in; ii++){
    (*box_cat)[ii].np = 0;
    (*box_cat)[ii].bx = NULL;    
  }
}


void init_inter_struct(inter *interval){
    (*interval).ns_in    = 0;
    (*interval).ns_tot   = 0;
    (*interval).n_lim    = 0;
    (*interval).first_sn = 0;
    (*interval).last_sn  = 0;
    (*interval).snap_num = NULL;
    (*interval).snap_rds = NULL;
    (*interval).lim_rds  = NULL;
    (*interval).lim_dis  = NULL;
}

void init_mass_sel(mass_sel **sm_sel){
  int ii=0;
  (*sm_sel) = NULL;
  (*sm_sel) = (mass_sel *)calloc(NPART,sizeof(mass_sel));
  for(ii=0; ii<NPART; ii++){
    (*sm_sel)[ii].i_obj   = NULL;    
    (*sm_sel)[ii].f_obj   = NULL;    
    (*sm_sel)[ii].tot_obj = NULL;    
    (*sm_sel)[ii].prob    = NULL;
    (*sm_sel)[ii].tot_obj = (int *)calloc(MASS_NBIN,sizeof(int));
    (*sm_sel)[ii].i_obj   = (int *)calloc(MASS_NBIN,sizeof(int));
    (*sm_sel)[ii].f_obj   = (float *)calloc(MASS_NBIN,sizeof(float));
    (*sm_sel)[ii].prob    = (float *)calloc(MASS_NBIN,sizeof(float));
  
    if((*sm_sel)[ii].i_obj==NULL) error_malloc_allocation("sm_sel.i_obj");
    if((*sm_sel)[ii].f_obj==NULL) error_malloc_allocation("sm_sel.f_obj");
    if((*sm_sel)[ii].tot_obj==NULL) error_malloc_allocation("sm_sel.tot_obj");
    if((*sm_sel)[ii].prob==NULL) error_malloc_allocation("sm_sel.prob_obj");
  }
}

void init_sheels_ham(hl_proxy **ham_px){
  int ii=0;
  (*ham_px) = NULL;
  (*ham_px) = (hl_proxy *)malloc(NPART*sizeof(hl_proxy));
  if((*ham_px)==NULL) error_malloc_allocation("ham_px");
    for(ii=0; ii<NPART; ii++){
    (*ham_px)[ii].nsel = 0;
    (*ham_px)[ii].np   = 0;
    (*ham_px)[ii].vol  = 0.0;
    (*ham_px)[ii].prx  = NULL;    
  }
}

void allocate_sheels_ham(hl_proxy **ham_px,lightcone lc_cat){
  long ii=0;
  float BINSIZE = (MAXZ-MINZ)/(float)NPART;
  float INV_BINSIZE = (float)NPART/(MAXZ-MINZ);

  for(ii=0; ii<lc_cat.np; ii++){
    int itemp=(int)((lc_cat.lc[ii].rds-MINZ)*INV_BINSIZE);
    if(itemp<0 || itemp>=NPART)
      continue;	
    (*ham_px)[itemp].np++;
  }
  
  read_radial_nbar(FALSE,NULL);
  for(ii=0; ii<NPART; ii++){
    float MIN_L = MINZ + (float)ii*BINSIZE;
    float MAX_L = MINZ + (float)(ii+1)*BINSIZE;
    float MIN_R = com_distance(MIN_L);
    float MAX_R = com_distance(MAX_L);
    float z = (MAX_L+MIN_L)*0.5;
    float vol = Volume(MIN_R,MAX_R);
    long size = (*ham_px)[ii].np;
    (*ham_px)[ii].vol  = vol;
    (*ham_px)[ii].nsel = (int)(number_density(z)*vol);
    (*ham_px)[ii].prx  = (proxy *)malloc(size*sizeof(proxy));
  }
}

void allocate_mem_lightcone(long size,lightcone *light_cat){
  (*light_cat).np = size;
  (*light_cat).lc = (lc_data *)malloc(size*sizeof(lc_data));

  if((*light_cat).lc==NULL) error_malloc_allocation("light_cat.lc");
}

void allocate_mem_random(long size,ran_lc *light_cat){
  (*light_cat).np = size;
  (*light_cat).lc = (ran_data *)malloc(size*sizeof(ran_data));
  
  if((*light_cat).lc==NULL) error_malloc_allocation("light_ran.lc");
}

void allocate_mem_box(long size,box *box_cat){
   
  (*box_cat).np = size;
  (*box_cat).bx =  (box_data *)malloc(size*sizeof(box_data));
  
  if((*box_cat).bx==NULL) error_malloc_allocation("box_cat.bx");

}

void load_lightcone(lightcone *light_cat){

  fprintf(stderr,"Reading input lightcone... ");

  if(READ_HDF5==TRUE){
    read_hdf5_lightcone(&(*light_cat));
  }
  else{
    long crows=0;
  
    crows = n_halos_ascii(INFILE);
    allocate_mem_lightcone(crows,light_cat);    
    read_lightcone_ascii(*light_cat);
  }
  timer(1);
}

void load_box(box *box_cat){

  fprintf(stderr,"Reading input boxes... ");

  if(READ_HDF5==TRUE){
    if(ROTBOX==1){Xh5="x";Yh5="y";Zh5="z";Vxh5="vx";Vyh5="vy";Vzh5="vz";}
    if(ROTBOX==2){Zh5="x";Xh5="y";Yh5="z";Vzh5="vx";Vxh5="vy";Vyh5="vz";}
    if(ROTBOX==3){Yh5="x";Zh5="y";Xh5="z";Vyh5="vx";Vzh5="vy";Vxh5="vz";}
    if(ROTBOX==4){Xh5="x";Yh5="y";Zh5="z";Vxh5="vx";Vyh5="vy";Vzh5="vz";}
    if(ROTBOX==5){Zh5="x";Xh5="y";Yh5="z";Vzh5="vx";Vxh5="vy";Vyh5="vz";}
    if(ROTBOX==6){Yh5="x";Zh5="y";Xh5="z";Vyh5="vx";Vzh5="vy";Vxh5="vz";}
    read_hdf5_box(INFILE,&(*box_cat));
  }
  else{
    long crows=0;
  
    crows = n_halos_ascii(INFILE);
    allocate_mem_box(crows,&(*box_cat));
    read_box_ascii(INFILE,*box_cat);
  }
  timer(1);
}

void load_boxes(inter interv,box **box_cat){

  int ii=0;
  
  fprintf(stderr,"Reading input boxes... ");

  if(READ_HDF5==TRUE){
    if(ROTBOX==1){Xh5="x";Yh5="y";Zh5="z";Vxh5="vx";Vyh5="vy";Vzh5="vz";}
    if(ROTBOX==2){Zh5="x";Xh5="y";Yh5="z";Vzh5="vx";Vxh5="vy";Vyh5="vz";}
    if(ROTBOX==3){Yh5="x";Zh5="y";Xh5="z";Vyh5="vx";Vzh5="vy";Vxh5="vz";}
    if(ROTBOX==4){Xh5="x";Yh5="y";Zh5="z";Vxh5="vx";Vyh5="vy";Vzh5="vz";}
    if(ROTBOX==5){Zh5="x";Xh5="y";Yh5="z";Vzh5="vx";Vxh5="vy";Vyh5="vz";}
    if(ROTBOX==6){Yh5="x";Zh5="y";Xh5="z";Vyh5="vx";Vzh5="vy";Vxh5="vz";}
    for(ii = 0; ii < interv.ns_in; ii++ ){
      char file[CHAR_SIZE];
      int SNAP_NUM=interv.snap_num[ii+interv.first_sn];
      
      sprintf(file,"%s%02d%s",BASEFILE,SNAP_NUM,ENDFILE);
      read_hdf5_box(file,&(*box_cat)[ii]);
    }
  }
  else{
    
#pragma omp parallel for
    for(ii = 0; ii < interv.ns_in; ii++ ){
      char file[CHAR_SIZE];
      long crows=0;
      int SNAP_NUM=interv.snap_num[ii+interv.first_sn];
      
      sprintf(file,"%s%02d%s",BASEFILE,SNAP_NUM,ENDFILE);
      
      crows = n_halos_ascii(file);
      allocate_mem_box(crows,&(*box_cat)[ii]);
      read_box_ascii(file,(*box_cat)[ii]);
    }    
  }
  timer(1);
}

void free_inter(inter interv){
  if(interv.ns_tot>0){
    if(interv.snap_num){free(interv.snap_num); interv.snap_num=NULL;}
    if(interv.snap_rds){free(interv.snap_rds); interv.snap_rds=NULL;}
  }
  if(interv.ns_in>0){
    if(interv.lim_rds){free(interv.lim_rds); interv.lim_rds=NULL;}
    if(interv.lim_dis){free(interv.lim_dis); interv.lim_dis=NULL;}
  }
}

void free_boxes(int num,box *cat){
  int ii=0;
  for(ii=0; ii<num; ii++)
    free_box(cat[ii]);
  if(cat){free(cat); cat=NULL;}
}

void free_box(box cat){
  if(cat.np>0) {
    if(cat.bx){free(cat.bx); cat.bx=NULL;}
  }
}

void free_lightcone(lightcone cat){
  if(cat.np>0) {
    if(cat.lc){free(cat.lc); cat.lc=NULL;}
  }
}


void free_random(ran_lc cat){
  if(cat.np>0) {
    if(cat.lc){free(cat.lc); cat.lc=NULL;}
  }
}

void free_gsl_nbar(){
  gsl_interp_accel_free(Acc_nbar);
  gsl_spline_free(spline_nbar);
}

void free_gsl_NORMnbar(){
  gsl_interp_accel_free(Acc_NORMnbar);
  gsl_spline_free(spline_NORMnbar);
}

void free_gsl_CMF(){
  gsl_interp_accel_free(Acc_CMF);
  gsl_spline_free(spline_CMF);
}
