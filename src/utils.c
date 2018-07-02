#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <rng_gsl.c>
#include <gsl/gsl_integration.h>
#include "hdf5.h"
#include "common.h"
#include "define.h"
#include "cosmo_tools.h"

static time_t absbeg,absend;

void ini_conf_values(){
  int ii=0;
  for(ii=0;ii<N_VAR_LC;ii++)
    PAR_VAL_LC[ii]=-1;
  for(ii=0;ii<N_VAR_BX;ii++)
    PAR_VAL_BX[ii]=-1;
}

void ini_var_conf_lc(){
  int ii;
  for(ii=0;ii<N_VAR_LC;ii++){
    if(!strcmp(PAR_NAME_LC[ii],"MASK")) MASK=PAR_VAL_LC[ii];
    if(!strcmp(PAR_NAME_LC[ii],"RANDOM")) RANDOM=PAR_VAL_LC[ii];
    if(!strcmp(PAR_NAME_LC[ii],"READ_HDF5")) READ_HDF5=PAR_VAL_LC[ii];
    if(!strcmp(PAR_NAME_LC[ii],"WRITE_HDF5")) WRITE_HDF5=PAR_VAL_LC[ii];
    if(!strcmp(PAR_NAME_LC[ii],"READ_LC")) READ_LC=PAR_VAL_LC[ii];
    if(!strcmp(PAR_NAME_LC[ii],"PART_CAT")) PART_CAT=PAR_VAL_LC[ii];
    if(!strcmp(PAR_NAME_LC[ii],"HAM")) HAM=PAR_VAL_LC[ii];
    if(!strcmp(PAR_NAME_LC[ii],"HOD")) HOD=PAR_VAL_LC[ii];
    if(!strcmp(PAR_NAME_LC[ii],"MODIFIED_HAM")) MODIFIED_HAM=PAR_VAL_LC[ii];
    if(!strcmp(PAR_NAME_LC[ii],"DOWNSAMPLE")) DOWNSAMPLE=PAR_VAL_LC[ii];
    if(!strcmp(PAR_NAME_LC[ii],"COR_FUN")) COR_FUN=PAR_VAL_LC[ii];
    if(!strcmp(PAR_NAME_LC[ii],"FIND_PARM")) FIND_PARM=PAR_VAL_LC[ii];
    if(!strcmp(PAR_NAME_LC[ii],"VAR_CUT")) VAR_CUT=PAR_VAL_LC[ii];
  }
}

void ini_var_conf_box(){
  int ii;
  for(ii=0;ii<N_VAR_BX;ii++){
    if(!strcmp(PAR_NAME_BX[ii],"WRITE_HDF5")) WRITE_HDF5=PAR_VAL_BX[ii];
    if(!strcmp(PAR_NAME_BX[ii],"READ_HDF5")) READ_HDF5=PAR_VAL_BX[ii];
    if(!strcmp(PAR_NAME_BX[ii],"HAM")) HAM=PAR_VAL_BX[ii];
    if(!strcmp(PAR_NAME_BX[ii],"HOD")) HOD=PAR_VAL_BX[ii];
    if(!strcmp(PAR_NAME_BX[ii],"MODIFIED_HAM")) MODIFIED_HAM=PAR_VAL_BX[ii];
    if(!strcmp(PAR_NAME_BX[ii],"DOWNSAMPLE")) DOWNSAMPLE=PAR_VAL_BX[ii];
    if(!strcmp(PAR_NAME_BX[ii],"COR_FUN")) COR_FUN=PAR_VAL_BX[ii];
    if(!strcmp(PAR_NAME_BX[ii],"FIND_PARM")) FIND_PARM=PAR_VAL_BX[ii];
    if(!strcmp(PAR_NAME_BX[ii],"VAR_CUT")) VAR_CUT=PAR_VAL_BX[ii];
  }
}

void ini_var_param(){
  int ii;
  for(ii=0;ii<N_VAR;ii++){
    if(!strcmp(VAR_NAME[ii],"OMEGA_M") && VAR_VALUE[ii]!=NULL)
      OMEGA_M=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"OMEGA_L") && VAR_VALUE[ii]!=NULL)
      OMEGA_L=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"HUBBLE_PAR") && VAR_VALUE[ii]!=NULL)
      HUBBLE_PAR=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"MINZ") && VAR_VALUE[ii]!=NULL)
      MINZ=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"MAXZ") && VAR_VALUE[ii]!=NULL)
      MAXZ=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"RA_LOS") && VAR_VALUE[ii]!=NULL)
      RA_LOS=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"DEC_LOS") && VAR_VALUE[ii]!=NULL)
      DEC_LOS=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"RA_SIZE") && VAR_VALUE[ii]!=NULL)
      RA_SIZE=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"DEC_SIZE") && VAR_VALUE[ii]!=NULL)
      DEC_SIZE=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"AREA") && VAR_VALUE[ii]!=NULL)
      AREA=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"NPART") && VAR_VALUE[ii]!=NULL)
      NPART=atoi(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"X_OBS") && VAR_VALUE[ii]!=NULL)
      X_OBS=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"Y_OBS") && VAR_VALUE[ii]!=NULL)
      Y_OBS=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"Z_OBS") && VAR_VALUE[ii]!=NULL)
      Z_OBS=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"LBOX") && VAR_VALUE[ii]!=NULL)
      LBOX=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"FRAN") && VAR_VALUE[ii]!=NULL)
      FRAN=atoi(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"ROTBOX") && VAR_VALUE[ii]!=NULL)
      ROTBOX=atoi(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"MIN_COMP") && VAR_VALUE[ii]!=NULL)
      MIN_COMP=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"SCAT_HAM") && VAR_VALUE[ii]!=NULL)
      SCAT_HAM=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"ZBOX") && VAR_VALUE[ii]!=NULL)
      ZBOX=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"NBAR") && VAR_VALUE[ii]!=NULL)
      NBAR=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"FSAT") && VAR_VALUE[ii]!=NULL)
      FSAT=atof(VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"ANG_DOWN") && VAR_VALUE[ii]!=NULL){
      if(!strcmp(VAR_VALUE[ii],"yes")) ANG_DOWN=TRUE;
      if(strcmp(VAR_VALUE[ii],"yes") && strcmp(VAR_VALUE[ii],"no")){
	fprintf(stderr,"ANG_DOWN only accepts \'yes\' or \'no\'\n"
		 "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
    }
    if(!strcmp(VAR_NAME[ii],"SHUFFLE") && VAR_VALUE[ii]!=NULL){
      if(!strcmp(VAR_VALUE[ii],"yes")) SHUFFLE=TRUE;
      if(strcmp(VAR_VALUE[ii],"yes") && strcmp(VAR_VALUE[ii],"no")){
	fprintf(stderr,"SHUFFLE only accepts \'yes\' or \'no\'\n"
		 "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
    }
    if(!strcmp(VAR_NAME[ii],"CSMF") && VAR_VALUE[ii]!=NULL){
      if(!strcmp(VAR_VALUE[ii],"yes")) CSMF=TRUE;
      if(strcmp(VAR_VALUE[ii],"yes") && strcmp(VAR_VALUE[ii],"no")){
	fprintf(stderr,"CSMF only accepts \'yes\' or \'no\'\n"
		 "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
    }
    
    if(!strcmp(VAR_NAME[ii],"CONF_BOX") && VAR_VALUE[ii]!=NULL)
      strcpy(CONF_BOX,VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"INFILE") && VAR_VALUE[ii]!=NULL)
      strcpy(INFILE,VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"BASEFILE") && VAR_VALUE[ii]!=NULL)
      strcpy(BASEFILE,VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"ENDFILE") && VAR_VALUE[ii]!=NULL)
      strcpy(ENDFILE,VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"OUTFILE") && VAR_VALUE[ii]!=NULL)
      strcpy(OUTFILE,VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"RANFILE") && VAR_VALUE[ii]!=NULL)
      strcpy(RANFILE,VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"NBARFILE") && VAR_VALUE[ii]!=NULL)
      strcpy(NBARFILE,VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"F_RDS_SIM") && VAR_VALUE[ii]!=NULL)
      strcpy(F_RDS_SIM,VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"MASKFILE") && VAR_VALUE[ii]!=NULL)
      strcpy(MASKFILE,VAR_VALUE[ii]);
    if(!strcmp(VAR_NAME[ii],"SMFCAT") && VAR_VALUE[ii]!=NULL)
      strcpy(SMFCAT,VAR_VALUE[ii]);
  }
}
 
void n_good_sanpshots(inter *interv){

  if((*interv).ns_tot>1){
    int ii=0,nn=0,nz=0;
    double step1,step2;
    
    while((*interv).snap_rds[ii] >= MAXZ) ii++;
    if(ii==0){
      nz++;
      (*interv).first_sn=0;
    }
    else{
      step1=(*interv).snap_rds[ii-1]-MAXZ;
      step2=MAXZ-(*interv).snap_rds[ii];
      if(step1<step2){
	nz=nz+2;
	(*interv).first_sn=ii-1;
      }
      else{
	nz++;
	(*interv).first_sn=ii;
      }
    }

    nn=ii+1;
    while((*interv).snap_rds[nn] > MINZ){
      nz++;
      (*interv).last_sn=nn;
      nn++;
    }

    step1=(*interv).snap_rds[nn-1]-MINZ;
    step2=MINZ-(*interv).snap_rds[nn];
    if(step1>step2){
      nz++;
      (*interv).last_sn=nn;
    }
    else (*interv).last_sn=nn-1;
    
    (*interv).ns_in=nz;
    if((*interv).last_sn<(*interv).first_sn)
      (*interv).last_sn=(*interv).first_sn;
  }
  else
    (*interv).ns_in=1;
}

void select_intervals(inter *interv){
   
  int ii=0,jj=1;
  
  (*interv).lim_rds[0]=MAXZ;
  for(ii=(*interv).first_sn;ii < (*interv).last_sn; ii++){
    (*interv).lim_rds[jj]
      =((*interv).snap_rds[ii]+(*interv).snap_rds[ii+1])*0.5;
    jj++;
  }
  (*interv).lim_rds[jj]=MINZ;
  
}

void interval_distances(inter *interv){
   
   int ii=0;
   for(ii=0; ii<=(*interv).ns_in; ii++){
     (*interv).lim_dis[ii]=com_distance((*interv).lim_rds[ii]);
   } 
}

void set_redshift_intervals(inter *interv){

  int lines=(int)n_halos_ascii(F_RDS_SIM);
  (*interv).ns_tot = lines;
  (*interv).snap_num = (int *)malloc(lines*sizeof(int));
  (*interv).snap_rds = (float *)malloc(lines*sizeof(float));
  if((*interv).snap_num==NULL) error_malloc_allocation("interv.snap_num");
  if((*interv).snap_rds==NULL) error_malloc_allocation("interv.snap_rds");
  
  read_redshift((*interv));
  n_good_sanpshots(interv);
  
  int tot_lims=(*interv).ns_in+1;
  (*interv).n_lim = tot_lims;
  (*interv).lim_rds = (float *)malloc(tot_lims*sizeof(float));
  (*interv).lim_dis = (float *)malloc(tot_lims*sizeof(float));
  
  select_intervals(interv);
}

/*float fast_atof (char *p){
  int frac;
  float sign, value, scale;
  
  while ((*p) == ' ' || (*p) == '\t' ) p += 1;
  
  sign = 1.0;
  if (*p == '-'){
    sign = -1.0;
    p += 1;
  }
  else if (*p == '+')
    p += 1;

  for (value = 0.0; ((*p)>='0'&&(*p)<='9'); p += 1) 
    value = value * 10.0 + (*p - '0');
  
  if (*p == '.') {
    float pow10 = 10.0;
    p += 1;
    while ((*p)>='0'&& (*p)<='9') {
      value += (*p - '0') / pow10;
      pow10 *= 10.0;
      p += 1;
    }
  }

  frac = 0;
  scale = 1.0;
  if ((*p == 'e') || (*p == 'E')) {
    unsigned int expon;

    p += 1;
    if (*p == '-') {
      frac = 1;
      p += 1;
    }
    else if (*p == '+')
      p += 1;

    for (expon = 0; ((*p)>='0'&&(*p)<='9'); p += 1) 
      expon = expon * 10 + (*p - '0');
    
    if (expon > 308) expon = 308;
    
    while (expon >= 50) { scale *= 1E50; expon -= 50; }
    while (expon >=  8) { scale *= 1E8;  expon -=  8; }
    while (expon >   0) { scale *= 10.0; expon -=  1; }
  }
  
  return sign * (frac ? (value / scale) : (value * scale));
  }*/

int comp_var(const void *p1, const void *p2){

   int j=1000000;
   proxy* pvel1=(proxy*)p1;
   proxy* pvel2=(proxy*)p2;
   return (int) (pvel2->var*j-pvel1->var*j);
}

void fill_sheels(hl_proxy *ham_px,lightcone lc_cat){
  int ii=0;
  int idx[NPART];
  float INV_BINSIZE = (float)NPART/(MAXZ-MINZ);
  unsigned int seed = (int)time(NULL);
  void *rng;
  rng = rng_init(seed);
  memset(idx,0,NPART*sizeof(int));
 
  for(ii=0; ii<lc_cat.np; ii++){
    int itemp=(int)((lc_cat.lc[ii].rds-MINZ)*INV_BINSIZE);
    if(itemp<0 || itemp>=NPART)
      continue;	
    float GAUSSRAN = gsl_ran_gaussian (rng, SCAT_HAM);
    ham_px[itemp].prx[idx[itemp]].pst = ii;
    ham_px[itemp].prx[idx[itemp]].var = lc_cat.lc[ii].Vpeak*(1.0+GAUSSRAN);
    idx[itemp]++;
  }
  
#pragma omp parallel for schedule(dynamic)
  for(ii=0; ii<NPART; ii++){
    qsort(ham_px[ii].prx,ham_px[ii].np,sizeof(proxy),comp_var);
  }
  
  rng_kill( rng );
}

void fill_proxy_box(proxy *ham_px,box bx_cat){
  int ii=0;
  unsigned int seed = (int)time(NULL);
  void *rng;
  rng = rng_init(seed);

  for(ii=0; ii<bx_cat.np; ii++){
    float GAUSSRAN = gsl_ran_gaussian (rng, SCAT_HAM);
    ham_px[ii].pst = ii;
    ham_px[ii].var = bx_cat.bx[ii].Vpeak*(1.0+GAUSSRAN);
  }

  qsort(ham_px,bx_cat.np,sizeof(proxy),comp_var);
  
  rng_kill( rng );
}

void timer(int ii){
  // timer(0) -> initialize absolute clock
  // timer(1) -> read absolute clock
  if(ii==0)
    absbeg=omp_get_wtime();
  else if(ii==1) {
    absend=omp_get_wtime();
    fprintf(stderr,"  %2dm%2ds\n",(int)((absend-absbeg)/60),(int)(absend-absbeg)%60);
  }
}

void write_hdf5(lightcone lc_cat){
  char        FILE[CHAR_SIZE];
  hid_t       file, memtype, strtype, space, dset;
  hsize_t     dims[1] = {lc_cat.np};
  int ii=0;
  hdf5_data   *h5temp;   
  h5temp  =  (hdf5_data *)malloc(sizeof(hdf5_data)*lc_cat.np);
#pragma omp parallel
  {
#pragma omp for nowait schedule(dynamic)
    for(ii=0;ii<lc_cat.np;ii++){
      
      h5temp[ii].id      = lc_cat.lc[ii].id;
      h5temp[ii].ra      = lc_cat.lc[ii].ra;
      h5temp[ii].dec     = lc_cat.lc[ii].dec;
      h5temp[ii].rds     = lc_cat.lc[ii].rds;
      h5temp[ii].zreal   = lc_cat.lc[ii].zreal;
      h5temp[ii].Mv      = lc_cat.lc[ii].Mv;
      h5temp[ii].Rv      = lc_cat.lc[ii].Rv;
      h5temp[ii].M200    = lc_cat.lc[ii].M200;
      h5temp[ii].Vmax    = lc_cat.lc[ii].Vmax;
      h5temp[ii].Vpeak   = lc_cat.lc[ii].Vpeak;
      h5temp[ii].str_id  = lc_cat.lc[ii].str_id;
      h5temp[ii].lc_id   = lc_cat.lc[ii].lc_id;
      h5temp[ii].boxnum  = lc_cat.lc[ii].boxnum;
      h5temp[ii].snapnum = lc_cat.lc[ii].snapnum;
    }
  }
  
  sprintf(FILE,"%s%s%.2f%s%.2f%s",OUTFILE,".minz.",MINZ,".maxz.",MAXZ,".h5");
  file = H5Fcreate (FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  space = H5Screate_simple (1, dims, NULL);

  strtype = H5Tcopy (H5T_C_S1);
  H5Tset_size (strtype, H5T_VARIABLE);
  
  memtype = H5Tcreate (H5T_COMPOUND, sizeof (hdf5_data));
  H5Tinsert (memtype, "ra",HOFFSET(hdf5_data,ra), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "dec",HOFFSET (hdf5_data,dec), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "rds",HOFFSET (hdf5_data,rds), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "zreal",HOFFSET (hdf5_data,zreal), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Mv",HOFFSET (hdf5_data,Mv), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Rv",HOFFSET (hdf5_data,Rv), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "M200",HOFFSET (hdf5_data,M200), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Vmax",HOFFSET (hdf5_data,Vmax), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Vpeak",HOFFSET (hdf5_data,Vpeak), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "str_id",HOFFSET (hdf5_data,str_id), H5T_NATIVE_LONG);
  H5Tinsert (memtype, "id",HOFFSET (hdf5_data,id), H5T_NATIVE_LONG);
  H5Tinsert (memtype, "lc_id",HOFFSET (hdf5_data,lc_id), H5T_NATIVE_LONG);
  H5Tinsert (memtype, "boxnum",HOFFSET (hdf5_data,boxnum), H5T_NATIVE_LONG);
  H5Tinsert (memtype, "snapnum",HOFFSET (hdf5_data,snapnum), H5T_NATIVE_INT);

  dset = H5Dcreate (file, DATASET, memtype, space, H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);

  H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,h5temp);

  H5Dclose (dset);
  H5Sclose (space);
  H5Tclose (memtype);
  H5Fclose (file);
  H5Tclose (strtype);
  free(h5temp);
}

void write_ascii(lightcone dat){
  FILE *fd;
  int  ii=0;
  char fname[CHAR_SIZE];

  sprintf(fname,"%s%s%.2f%s%.2f%s",OUTFILE,".minz.",MINZ,".maxz.",MAXZ,".dat");
  fd = fopen(fname,"w");
  for(ii=0; ii<dat.np; ii++){
    fprintf(fd,"%.6f %.6f %.6f %.6f %.5e %.5e %.6f %.6f %.6f %ld %ld %d %ld\n", \
	    dat.lc[ii].ra,dat.lc[ii].dec,dat.lc[ii].rds,dat.lc[ii].zreal,
	    dat.lc[ii].Mv,dat.lc[ii].M200,dat.lc[ii].Rv,dat.lc[ii].Vmax,
	    dat.lc[ii].Vpeak,dat.lc[ii].str_id,dat.lc[ii].id,dat.lc[ii].snapnum,
	    dat.lc[ii].boxnum); 
  }

  
  /*for(ii=0; ii<dat.np; ii++){
    fprintf(fd,"%.6f %.6f %.6f %.6f\n", \
	    dat.lc[ii].ra,dat.lc[ii].dec,dat.lc[ii].zreal,dat.lc[ii].rds);
  }*/

  fclose(fd);  
}

void write_hdf5_sel(lightcone lc_cat){
  char        FILE[CHAR_SIZE];
  hid_t       file, memtype, strtype, space, dset;
  hsize_t     dims[1] = {NSEL};
  int ii=0,jj=0;
  hdf5_data   *h5temp;   
  h5temp  =  (hdf5_data *)malloc(sizeof(hdf5_data)*NSEL);
  for(ii=0;ii<lc_cat.np;ii++){
    if(lc_cat.lc[ii].sel){
      h5temp[jj].id      = lc_cat.lc[ii].id;
      h5temp[jj].ra      = lc_cat.lc[ii].ra;
      h5temp[jj].dec     = lc_cat.lc[ii].dec;
      h5temp[jj].rds     = lc_cat.lc[ii].rds;
      h5temp[jj].zreal   = lc_cat.lc[ii].zreal;
      h5temp[jj].wfkp    = lc_cat.lc[ii].wfkp;
      h5temp[jj].Mv      = lc_cat.lc[ii].Mv;
      h5temp[jj].Mstar   = lc_cat.lc[ii].Mstar;
      h5temp[jj].Rv      = lc_cat.lc[ii].Rv;
      h5temp[jj].M200    = lc_cat.lc[ii].M200;
      h5temp[jj].Vmax    = lc_cat.lc[ii].Vmax;
      h5temp[jj].Vpeak   = lc_cat.lc[ii].Vpeak;
      h5temp[jj].str_id  = lc_cat.lc[ii].str_id;
      h5temp[jj].lc_id   = lc_cat.lc[ii].lc_id;
      h5temp[jj].boxnum  = lc_cat.lc[ii].boxnum;
      h5temp[jj].snapnum = lc_cat.lc[ii].snapnum;
      jj++;
    }
  }
  
  sprintf(FILE,"%s%s%.2f%s%.2f%s",OUTFILE,".minz.",MINZ,".maxz.",MAXZ,".h5");
  file = H5Fcreate (FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  space = H5Screate_simple (1, dims, NULL);

  strtype = H5Tcopy (H5T_C_S1);
  H5Tset_size (strtype, H5T_VARIABLE);
  
  memtype = H5Tcreate (H5T_COMPOUND, sizeof (hdf5_data));
  H5Tinsert (memtype, "ra",HOFFSET(hdf5_data,ra), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "dec",HOFFSET (hdf5_data,dec), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "rds",HOFFSET (hdf5_data,rds), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "zreal",HOFFSET (hdf5_data,zreal), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Mv",HOFFSET (hdf5_data,Mv), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Rv",HOFFSET (hdf5_data,Rv), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "M200",HOFFSET (hdf5_data,M200), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Vmax",HOFFSET (hdf5_data,Vmax), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Vpeak",HOFFSET (hdf5_data,Vpeak), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "str_id",HOFFSET (hdf5_data,str_id), H5T_NATIVE_LONG);
  H5Tinsert (memtype, "id",HOFFSET (hdf5_data,id), H5T_NATIVE_LONG);
  H5Tinsert (memtype, "lc_id",HOFFSET (hdf5_data,lc_id), H5T_NATIVE_LONG);
  H5Tinsert (memtype, "boxnum",HOFFSET (hdf5_data,boxnum), H5T_NATIVE_LONG);
  H5Tinsert (memtype, "snapnum",HOFFSET (hdf5_data,snapnum), H5T_NATIVE_INT);
  H5Tinsert (memtype, "wfkp",HOFFSET (hdf5_data,wfkp), H5T_NATIVE_FLOAT);
  if(CSMF || !NO_SEL)
    H5Tinsert (memtype, "Mstar",HOFFSET (hdf5_data,Mstar), H5T_NATIVE_FLOAT);

  dset = H5Dcreate (file, DATASET, memtype, space, H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);

  H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,h5temp);

  H5Dclose (dset);
  H5Sclose (space);
  H5Tclose (memtype);
  H5Fclose (file);
  H5Tclose (strtype);
  free(h5temp);
}

void write_ascii_sel(lightcone dat){
  FILE *fd;
  int  ii=0;
  char fname[CHAR_SIZE];

  sprintf(fname,"%s%s%.2f%s%.2f%s",OUTFILE,".minz.",MINZ,".maxz.",MAXZ,".dat");
  fd = fopen(fname,"w");
  if(!CSMF || !HAM)
    for(ii=0; ii<dat.np; ii++)
      if(dat.lc[ii].sel)
	fprintf(fd,"%.6f %.6f %.6f %.6f %.6f %.5e %.5e %.6f %.6f %.6f %ld"
		" %ld %d %ld\n",dat.lc[ii].ra,dat.lc[ii].dec,dat.lc[ii].rds,
		dat.lc[ii].wfkp,dat.lc[ii].zreal,dat.lc[ii].Mv,dat.lc[ii].M200,
		dat.lc[ii].Rv,dat.lc[ii].Vmax,dat.lc[ii].Vpeak,dat.lc[ii].str_id,
		dat.lc[ii].id,dat.lc[ii].snapnum,dat.lc[ii].boxnum);
  

  if(CSMF && HAM)
    for(ii=0; ii<dat.np; ii++)
      if(dat.lc[ii].sel)
	fprintf(fd,"%.6f %.6f %.6f %.6f %.6f %.6f %.5e %.5e %.6f %.6f %.6f %ld"
		" %ld %d %ld\n",dat.lc[ii].ra,dat.lc[ii].dec,dat.lc[ii].rds,
		dat.lc[ii].wfkp,dat.lc[ii].Mstar,dat.lc[ii].zreal,dat.lc[ii].Mv,
		dat.lc[ii].M200,dat.lc[ii].Rv,dat.lc[ii].Vmax,dat.lc[ii].Vpeak,
		dat.lc[ii].str_id,dat.lc[ii].id,dat.lc[ii].snapnum,dat.lc[ii].boxnum);

  fclose(fd);  
}

void write_hdf5_sel_box(box box_cat){
  /*char        FILE[CHAR_SIZE];
  hid_t       file, memtype, strtype, space, dset;
  hsize_t     dims[1] = {NSEL};
  int ii=0,jj=0;
  hdf5_box_w   *h5temp;   
  h5temp  =  (hdf5_data *)malloc(sizeof(hdf5_data)*NSEL);
  for(ii=0;ii<lc_cat.np;ii++){
    if(box_cat.bx[ii].sel){
      h5temp[jj].id      = box_cat.bx[ii].id;
      h5temp[jj].x       = box_cat.bx[ii].x;
      h5temp[jj].y       = box_cat.bx[ii].y;
      h5temp[jj].z       = box_cat.bx[ii].z;
      h5temp[jj].vx      = box_cat.bx[ii].vx;
      h5temp[jj].vy      = box_cat.bx[ii].vy;
      h5temp[jj].vz      = box_cat.bx[ii].vz;
      h5temp[jj].Mv      = box_cat.bx[ii].Mv;
      h5temp[jj].Mstar   = box_cat.bx[ii].Mstar;
      h5temp[jj].Rv      = box_cat.bx[ii].Rv;
      h5temp[jj].M200    = box_cat.bx[ii].M200;
      h5temp[jj].Vmax    = box_cat.bx[ii].Vmax;
      h5temp[jj].Vpeak   = box_cat.bx[ii].Vpeak;
      h5temp[jj].str_id  = box_cat.bx[ii].str_id;
      h5temp[jj].box_id  = box_cat.bx[ii].box_id;
      jj++;
    }
  }
   
  sprintf(FILE,"%s%s",OUTFILE,".h5");
  file = H5Fcreate (FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  space = H5Screate_simple (1, dims, NULL);

  strtype = H5Tcopy (H5T_C_S1);
  H5Tset_size (strtype, H5T_VARIABLE);
  
  memtype = H5Tcreate (H5T_COMPOUND, sizeof (hdf5_box_w));
  H5Tinsert (memtype, "x",HOFFSET(hdf5_box_w,x), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "y",HOFFSET (hdf5_box_w,y), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "z",HOFFSET (hdf5_box_w,z), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "vx",HOFFSET (hdf5_box_w,vx), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "vy",HOFFSET (hdf5_box_w,vy), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "vz",HOFFSET (hdf5_box_w,vz), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Mv",HOFFSET (hdf5_box_w,Mv), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Rv",HOFFSET (hdf5_box_w,Rv), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "M200",HOFFSET (hdf5_box_w,M200), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Vmax",HOFFSET (hdf5_box_w,Vmax), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Vpeak",HOFFSET (hdf5_box_w,Vpeak), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "str_id",HOFFSET (hdf5_box_w,str_id), H5T_NATIVE_LONG);
  H5Tinsert (memtype, "id",HOFFSET (hdf5_box_w,id), H5T_NATIVE_LONG);
  H5Tinsert (memtype, "box_id",HOFFSET (hdf5_box_w,box_id), H5T_NATIVE_LONG);
  if(CSMF || !NO_SEL)
    H5Tinsert (memtype, "Mstar",HOFFSET (hdf5_box_w,Mstar), H5T_NATIVE_FLOAT);

  dset = H5Dcreate (file, DATASET, memtype, space, H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);

  H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,h5temp);

  H5Dclose (dset);
  H5Sclose (space);
  H5Tclose (memtype);
  H5Fclose (file);
  H5Tclose (strtype);
  free(h5temp);*/
}

void write_ascii_sel_box(box dat){
  FILE *fd;
  int  ii=0;
  char fname[CHAR_SIZE];

  sprintf(fname,"%s%s",OUTFILE,".dat");
  fd = fopen(fname,"w");
  if(!CSMF)
    for(ii=0; ii<dat.np; ii++)
      if(dat.bx[ii].sel)
	fprintf(fd,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.5e %.5e %.6f %.6f %.6f"
		" %ld %ld %ld\n",dat.bx[ii].x,dat.bx[ii].y,dat.bx[ii].z,
		dat.bx[ii].zrds,dat.bx[ii].vx,dat.bx[ii].vy,dat.bx[ii].vz,
		dat.bx[ii].Mv,dat.bx[ii].M200,dat.bx[ii].Rv,dat.bx[ii].Vmax,
		dat.bx[ii].Vpeak,dat.bx[ii].str_id,dat.bx[ii].id,
		dat.bx[ii].box_id);
  

  if(CSMF)
    for(ii=0; ii<dat.np; ii++)
      if(dat.bx[ii].sel)
	fprintf(fd,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.5e %.5e %.6f %.6f"
		" %.6f %ld %ld %ld\n",dat.bx[ii].x,dat.bx[ii].y,dat.bx[ii].z,
		dat.bx[ii].zrds,dat.bx[ii].vx,dat.bx[ii].vy,dat.bx[ii].vz,
		dat.bx[ii].Mstar,dat.bx[ii].Mv,dat.bx[ii].M200,dat.bx[ii].Rv,
		dat.bx[ii].Vmax,dat.bx[ii].Vpeak,dat.bx[ii].str_id,dat.bx[ii].id,
		dat.bx[ii].box_id);

		fclose(fd);  
}


void write_output_info(inter interv,int tot){
  FILE *fd;
  char fname[CHAR_SIZE];
  float rs=RA_SIZE*D2R;
  float ds=DEC_SIZE*D2R;
  float ld=DEC_LOS*D2R;
  float theta=sin(ld+ds*0.5)-sin(ld-ds*0.5);
  float area=rs*theta;
  sprintf(fname,"%s%s%.2f%s%.2f%s",OUTFILE,".minz.",MINZ,".maxz.",MAXZ,".info");
  fd = fopen(fname,"w");
  
  fprintf(fd,"***************\nSUGAR Code v1.0\n***************\n\n"
	  "Cosmological Parameters used:\n"
	  "Omega_m = %.6f\n"
	  "Omega_l = %.6f\n"
	  "h = %.6f\n\n",OMEGA_M,OMEGA_L,HUBBLE_PAR);
  
  if(!strcmp(CODE_OPT, "LGC"))
    fprintf(fd,"Light-Cone code\n\n"
	    "Minimum z = %.3f\n"
	    "Maximum z = %.3f\n",MINZ,MAXZ);
  if(NO_SEL)
    fprintf(fd,"\nAny selection method used\n"
	    "Lightcone built with all objects\n"
	    "Area =  %.5f\n",area);
  if(MASK){
    fprintf(fd,"\nApplying geometry from file:\n"
	    "%s\n",MASKFILE);
    if(ANG_DOWN)
      fprintf(fd,"Angular downsampling included\n");
    if(MIN_COMP!=FALSE)
      fprintf(fd,"Excluding regions with completeness smaller than %.3f:\n",
	      MIN_COMP);      
  }

  if(DOWNSAMPLE)
      fprintf(fd,"\nDOWNSAMPLING selection\n"
	    "nbar file: %s\n",NBARFILE);

  
	    fclose(fd);  
  }

/*Press-Schechter SMF used by Rodriguez-Torres et al . 2016*/
double f (double x, void * params) {
  int alpha = *(int *) params;
  double b[3];
  if(alpha==1){
    b[0]=0.0002663; b[1]=11.42;b[2]=2.447;}
  if(alpha==2){
    b[0]=0.004002; b[1]=10.76; b[2]=0.938;}
  double f = log(10)*b[0]*(pow(10,(x-b[1])*(1-b[2])))*exp(-pow(10,x-b[1]));
  return f;
}

double CSMF_integral(float smass){

  double int1=0.0,int2=0.0,error=0.0,total=0.0;
  size_t neval;
  int highm = 1;
  int lowm  = 2;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000);
    
  if(smass<=SMCUT){
    gsl_function G;
    G.function = &f;
    G.params = &lowm;
    gsl_integration_qag(&G,smass,SMCUT,0,1e-3,GSL_INTEG_GAUSS41,50000,w,&int1,&error);

    gsl_function F;
    F.function = &f;
    F.params = &highm;
    gsl_integration_qag(&F,SMCUT,MAX_MASS,0,1e-3,GSL_INTEG_GAUSS41,50000,w,&int2,&error);     
    total=int1+int2;
  }
  else{
    gsl_function F;
    F.function = &f;
    F.params = &highm;      
    gsl_integration_qng (&F,smass,MAX_MASS,0,1e-7,&total,&error,&neval); 
  }
  return total;
}

void CSMF_calculator(){
  int ii=0.0;
  int nP = 1000;
  double BINSIZE=(MAX_MASS-MIN_MASS)/(double)nP;
  double SM[nP],CSMF[nP];
  for(ii=1;ii<=nP;ii++){
    SM[nP-ii] = MIN_MASS+(ii-0.5)*BINSIZE;
    CSMF[nP-ii] = CSMF_integral(SM[nP-ii]);
  }

  Acc_CMF=gsl_interp_accel_alloc();
  spline_CMF=gsl_spline_alloc(gsl_interp_cspline,nP);
  gsl_spline_init(spline_CMF,CSMF,SM,nP);
}
