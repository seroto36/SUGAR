#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "hdf5.h"
#include <time.h>  
#include "cosmo_tools.h"
#include "common.h"
#include "define.h"

int find_maximum(double a[], double rds[], int n){
  int ii=0, index=0;
  double max=0.0;
 
  for (ii = 1; ii < n; ii++)
    if (a[ii] > max && rds[ii] >= MINZ && rds[ii] < MAXZ){
      index = ii;
      max = a[ii];
    }
 
  return index;
}

void split_line(int nline, char *line,char *name,
		char *value, char *file){
  
  char delims[] = " ";
  char *result = NULL; 
  char *var[3];
  int jj = 0;
  result = strtok( line, delims );
  while( result != NULL ) {
    var[jj]=result;
    result = strtok( NULL, delims );
    jj ++;
    if(jj>2){
      if(var[jj-1][0]=='#')
	break;
      else
	error_line_numvar(nline,file);}
  }
  if(jj<2) error_line_numvar(nline,file);
    
  strcpy(name, var[0]);
  strcpy(value, var[1]);
}

void read_redshift(inter interv){
   
   FILE *red;
   int ii=0;
   
   if(!(red = fopen(F_RDS_SIM, "r")))
     error_open_file(F_RDS_SIM);
   
   for(ii = 0; ii < interv.ns_tot; ii++ ){
     fscanf(red,"%d  %f",(interv.snap_num+ii),(interv.snap_rds+ii));
     if(ii>0 && interv.snap_rds[ii]>interv.snap_rds[ii-1]){
       fprintf(stderr,"Error in list of redshift file:\n"
	       "Redshifts have to be in descending order\n"
	       "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
   }
   
   fclose(red);
}

void normalise_nbar(int num,sheels parts[],double rds[],double nbar[]){
  int ii=0;
  float bin_size=(MAXZ-MINZ)/(int)NPART;
  for(ii=0; ii<num; ii++){
    int itemp=(int)((rds[ii]-MINZ)/bin_size);
    if(itemp<0){
      nbar[ii]=nbar[ii]/parts[0].nbar;
      continue;
    }
    if(itemp>=NPART){
      nbar[ii]=nbar[ii]/parts[NPART-1].nbar;
      continue;
    }
    nbar[ii]=nbar[ii]/parts[itemp].nbar;
  }
}

void read_radial_nbar(int NORM,sheels parts[]){

  FILE *fd;
  char buffer[CHAR_SIZE];
  int ii=0,nP=0;
  int id_max=0;
  
  if(!(fd = fopen(NBARFILE, "r")))
    error_open_file(NBARFILE);

  while(fgets(buffer,CHAR_SIZE,fd)!=NULL)nP++;
  rewind(fd);
  double rds[nP],nbar[nP];
  
  for(ii=0; ii<nP; ii++)
    fscanf(fd,"%lf %lf",&rds[ii],&nbar[ii]);

  fclose(fd);

  id_max=find_maximum(nbar,rds, nP);
  MAX_NBAR=nbar[id_max];

  Acc_nbar=gsl_interp_accel_alloc();
  spline_nbar=gsl_spline_alloc(gsl_interp_cspline,nP);
  gsl_spline_init(spline_nbar,rds,nbar,nP);
  
  if(NORM){
    normalise_nbar(nP,parts,rds,nbar);
    Acc_NORMnbar=gsl_interp_accel_alloc();
    spline_NORMnbar=gsl_spline_alloc(gsl_interp_cspline,nP);
    gsl_spline_init(spline_NORMnbar,rds,nbar,nP);
  }
}

void ReadandClass_SM(mass_sel *sm_sel,hl_proxy *ham_px){
  FILE *fd;
  float rds,smass,wgt;
  float tot_rds[NPART];
  float ZBIN_SIZE = (MAXZ-MINZ)/(float)NPART;
  float MBIN_SIZE = (MAX_MASS-MIN_MASS)/(float)MASS_NBIN;
  char buffer[CHAR_SIZE];
  int ii=0,jj=0,nP=0;

  memset(tot_rds,0,NPART*sizeof(int));
  
  if(!(fd = fopen(SMFCAT, "r")))
    error_open_file(SMFCAT);

  while(fgets(buffer,CHAR_SIZE,fd)!=NULL)nP++;
  rewind(fd);
  
  for(ii=0; ii<nP; ii++){
    fscanf(fd,"%f %f %f",&rds,&smass,&wgt);
    if(rds>=MAXZ || rds<MINZ || smass>=MAX_MASS || smass<MIN_MASS) continue;
    else{
      int z_idx = (int)((rds-MINZ)/ZBIN_SIZE);
      int m_idx = (int)((smass-MIN_MASS)/MBIN_SIZE);
      sm_sel[z_idx].f_obj[m_idx]+=wgt;
      tot_rds[z_idx]+=wgt;
    }
  }
  fclose(fd);

  for(ii=0;ii<NPART;ii++)
    for(jj=0;jj<MASS_NBIN;jj++)
      sm_sel[ii].f_obj[jj]=sm_sel[ii].f_obj[jj]*ham_px[ii].nsel/tot_rds[ii];
}

void read_conf_file(){

  FILE *file;
  char line[CHAR_SIZE];
  char name[CHAR_SIZE],value[CHAR_SIZE];
  char YES[]="yes",NO[]="no";
  int nline=1,cline=0;
  
  if((file=fopen(fconf, "r"))==NULL)
    error_open_file(fconf);
  
  /*Read lightcone params CODE_OPT=LGC*/
  if(!strcmp(CODE_OPT, "LGC")){

    fprintf(stderr,"\n**************\nLightcone code\n**************\n\n");

    while(fgets(line, sizeof(line), file) != NULL){
      if ((*line == '#') || (*line == '\n')){nline++; continue;}
      if(line[strlen(line) - 1] == '\n') line[strlen(line) - 1] = '\0';
      
      split_line(nline,line,name,value,fconf);
      for(cline = 0; cline < N_VAR_LC; cline++){
	if(!strcmp(PAR_NAME_LC[cline], name)){
	  if(PAR_VAL_LC[cline]!=DEFAULT)
	    error_repeat_param(nline, fconf, name);
	  if(!strcmp(YES, value)) PAR_VAL_LC[cline]=TRUE;
	  else{
	    if(!strcmp(NO, value)) PAR_VAL_LC[cline]=FALSE;
	    else error_line_cfile(nline);
	  }
	  break;
	}
      }
      if(cline==N_VAR_LC)
	error_input_param(nline,name,fconf);
      nline++;
    }

    for(cline = 0; cline < N_VAR_LC; cline++)
      if(PAR_VAL_LC[cline]==DEFAULT)
	error_missing_param(fconf,PAR_NAME_LC[cline]);

    ini_var_conf_lc(); //save results in configuration variables
    check_conf_lc();
  }

  /*Read box params CODE_OPT=BOX*/
  else{
    
    fprintf(stderr,"\n********\nBOX code\n********\n\n");

    while(fgets(line, sizeof(line), file) != NULL){
      if ((*line == '#') || (*line == '\n')){nline++; continue;}
      if(line[strlen(line) - 1] == '\n') line[strlen(line) - 1] = '\0';
      
      split_line(nline,line,name,value,fconf);
      for(cline = 0; cline < N_VAR_BX; cline++){
	if(!strcmp(PAR_NAME_BX[cline], name)){
	  if(PAR_VAL_BX[cline]!=DEFAULT)
	    error_repeat_param(nline, fconf, name);
	  if(!strcmp(YES, value)) PAR_VAL_BX[cline]=TRUE;
	  else{
	    if(!strcmp(NO, value)) PAR_VAL_BX[cline]=FALSE;
	    else error_line_cfile(nline);
	  }
	  break;
	}
      }
      if(cline==N_VAR_BX)
	error_input_param(nline,name,fconf);
      
      nline++;
    }

    for(cline = 0; cline < N_VAR_BX; cline++)
      if(PAR_VAL_BX[cline]==DEFAULT)
	error_missing_param(fconf,PAR_NAME_BX[cline]);

    ini_var_conf_box(); //save results in configuration variables
    check_conf_bx();
  }

}

void read_argc(int argc, char **argv){

  int n_argc=4;  
  if( argc != n_argc ) {
    fprintf( stderr, "USAGE:\t%s BOX/LGC config_file param_file\n", argv[0] );
    exit(EXIT_FAILURE);
  }
  if((!strcmp(argv[1],"LGC")) || (!strcmp(argv[1],"BOX"))){
    CODE_OPT = argv[1];
    fconf = argv[2];
    fparam = argv[3];}
  else{
    fprintf(stderr,"Unknown option %s\n", argv[1]);
    exit(EXIT_FAILURE);
  }
    
}

void read_param_file(){
  FILE *file;
  char line[CHAR_SIZE];
  char name[CHAR_SIZE],value[CHAR_SIZE];
  int nline=1,cline=0;
  
  if((file=fopen(fparam, "r"))==NULL)
    error_open_file(fparam);

  while(fgets(line, sizeof(line), file) != NULL){
    if ((*line == '#') || (*line == '\n')){nline++; continue;}
    if(line[strlen(line) - 1] == '\n') line[strlen(line) - 1] = '\0';
	  
    split_line(nline,line,name,value,fparam);
    for(cline = 0; cline < N_VAR; cline++){
      if(!strcmp(VAR_NAME[cline], name)){
	if(VAR_VALUE[cline]!=NULL){
	  error_repeat_param(nline, fparam, name);}
	VAR_VALUE[cline] = malloc (CHAR_SIZE);
	strcpy(VAR_VALUE[cline],value);
	break;
      }
    }

    if(cline==N_VAR)
      error_input_param(nline,name,fparam);

    nline++;
    
  }
  ini_var_param();
  check_params();
}

void read_hdf5_box(char *fname,box *box_cat){

  hdf5_box *read_h5data;
  hid_t    file,dset,space,memtype,strtype;
  hsize_t  dims[1];
  long     ii=0;
  float DIS=0.0,sig=0.0;

  if(ROTBOX==1){sig=1.0; DIS=0.0;}
  if(ROTBOX==2){sig=1.0; DIS=0.0;}
  if(ROTBOX==3){sig=1.0; DIS=0.0;}
  if(ROTBOX==4){sig=-1.0; DIS=LBOX;}
  if(ROTBOX==5){sig=-1.0; DIS=LBOX;}
  if(ROTBOX==6){sig=-1.0; DIS=LBOX;}

  strtype = H5Tcopy (H5T_C_S1);
  H5Tset_size (strtype, H5T_VARIABLE);

  memtype = H5Tcreate (H5T_COMPOUND, sizeof (hdf5_box));
  H5Tinsert (memtype, Xh5,HOFFSET(hdf5_box,x), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, Yh5,HOFFSET(hdf5_box,y), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, Zh5,HOFFSET(hdf5_box,z), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, Vxh5,HOFFSET(hdf5_box,vx), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, Vyh5,HOFFSET(hdf5_box,vy), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, Vzh5,HOFFSET(hdf5_box,vz), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Mv",HOFFSET(hdf5_box,Mv), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Rv",HOFFSET(hdf5_box,Rv), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "M200",HOFFSET(hdf5_box,M200), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Vmax",HOFFSET(hdf5_box,Vmax), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Vpeak",HOFFSET(hdf5_box,Vpeak), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "str_id",HOFFSET(hdf5_box,str_id), H5T_NATIVE_LONG);
  H5Tinsert (memtype, "id",HOFFSET(hdf5_box,id), H5T_NATIVE_LONG);
    
  file = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  dset = H5Dopen (file, DATABOX, H5P_DEFAULT);

  space = H5Dget_space (dset);
  H5Sget_simple_extent_dims (space, dims, NULL);
  read_h5data = (hdf5_box *) malloc (dims[0] * sizeof (hdf5_box));
  H5Dread (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, read_h5data);

  long size_cat     = (long)dims[0];
  allocate_mem_box(size_cat,box_cat);  
  
#pragma omp parallel for shared(DIS,sig)
  for(ii=0;ii<size_cat;ii++){
    (*box_cat).bx[ii].x      = DIS+sig*read_h5data[ii].x;
    (*box_cat).bx[ii].y      = DIS+sig*read_h5data[ii].y;
    (*box_cat).bx[ii].z      = DIS+sig*read_h5data[ii].z;
    (*box_cat).bx[ii].vx     = sig*read_h5data[ii].vx;
    (*box_cat).bx[ii].vy     = sig*read_h5data[ii].vy;
    (*box_cat).bx[ii].vz     = sig*read_h5data[ii].vz;
    (*box_cat).bx[ii].Mv     = read_h5data[ii].Mv;
    (*box_cat).bx[ii].M200   = read_h5data[ii].M200;
    (*box_cat).bx[ii].Rv     = read_h5data[ii].Rv;
    (*box_cat).bx[ii].Vpeak  = read_h5data[ii].Vpeak;
    (*box_cat).bx[ii].Vmax   = read_h5data[ii].Vmax;
    (*box_cat).bx[ii].str_id = read_h5data[ii].str_id;
    (*box_cat).bx[ii].id     = read_h5data[ii].id;
    (*box_cat).bx[ii].box_id = ii;
    (*box_cat).bx[ii].sel    = FALSE;
  }

  H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, read_h5data);
  H5Dclose (dset);
  H5Sclose (space);
  H5Tclose (memtype);
  
  H5Fclose (file);
  free (read_h5data);
  H5Tclose (strtype);
}

void read_hdf5_lightcone(lightcone *data){

  hdf5_data_r *read_h5data;
  hid_t       file,dset,space,memtype,strtype;
  hsize_t     dims[1];
  long        ii=0;

  strtype = H5Tcopy (H5T_C_S1);
  H5Tset_size (strtype, H5T_VARIABLE);

  memtype = H5Tcreate (H5T_COMPOUND, sizeof (hdf5_data_r));
  H5Tinsert (memtype, "ra",HOFFSET(hdf5_data_r,ra), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "dec",HOFFSET(hdf5_data_r,dec), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "rds",HOFFSET(hdf5_data_r,rds), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "zreal",HOFFSET(hdf5_data_r,zreal), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Mv",HOFFSET(hdf5_data_r,Mv), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "M200",HOFFSET(hdf5_data_r,M200), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Rv",HOFFSET(hdf5_data_r,Rv), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Vmax",HOFFSET(hdf5_data_r,Vmax), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "Vpeak",HOFFSET(hdf5_data_r,Vpeak), H5T_NATIVE_FLOAT);
  H5Tinsert (memtype, "str_id",HOFFSET(hdf5_data_r,str_id), H5T_NATIVE_LONG);
  H5Tinsert (memtype, "id",HOFFSET(hdf5_data_r,id), H5T_NATIVE_LONG);
  H5Tinsert (memtype, "snapnum",HOFFSET(hdf5_data_r,snapnum), H5T_NATIVE_INT);
  H5Tinsert (memtype, "boxnum",HOFFSET(hdf5_data_r,boxnum), H5T_NATIVE_LONG);
  H5Tinsert (memtype, "lc_id",HOFFSET(hdf5_data_r,lc_id), H5T_NATIVE_LONG);
    
  file = H5Fopen (INFILE, H5F_ACC_RDONLY, H5P_DEFAULT);
  dset = H5Dopen (file, DATASET, H5P_DEFAULT);

  space = H5Dget_space (dset);
  H5Sget_simple_extent_dims (space, dims, NULL);
  read_h5data = (hdf5_data_r *) malloc (dims[0] * sizeof (hdf5_data_r));
  H5Dread (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, read_h5data);

  long size_cat   = (long)dims[0];
  allocate_mem_lightcone(size_cat,data);
    
#pragma omp parallel for 
  for(ii=0;ii<size_cat;ii++){
    
    (*data).lc[ii].ra      = read_h5data[ii].ra;
    (*data).lc[ii].dec     = read_h5data[ii].dec;
    (*data).lc[ii].rds     = read_h5data[ii].rds;
    (*data).lc[ii].zreal   = read_h5data[ii].zreal;
    (*data).lc[ii].Mv      = read_h5data[ii].Mv;
    (*data).lc[ii].M200    = read_h5data[ii].M200;
    (*data).lc[ii].Rv      = read_h5data[ii].Rv;
    (*data).lc[ii].Vmax    = read_h5data[ii].Vmax;
    (*data).lc[ii].Vpeak   = read_h5data[ii].Vpeak;
    (*data).lc[ii].str_id  = read_h5data[ii].str_id;
    (*data).lc[ii].id      = read_h5data[ii].id;
    (*data).lc[ii].lc_id   = read_h5data[ii].lc_id;;
    (*data).lc[ii].boxnum  = read_h5data[ii].boxnum;
    (*data).lc[ii].snapnum = read_h5data[ii].snapnum;
    (*data).lc[ii].sel     = FALSE;
  }
  
  H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, read_h5data);
  H5Dclose (dset);
  H5Sclose (space);
  H5Tclose (memtype);
  
  H5Fclose (file);
  free (read_h5data);
  H5Tclose (strtype);
}

void read_lightcone_ascii(lightcone lc_cat){

   FILE *fd;
   long  n=0;
   char line[2*CHAR_SIZE];
   char delims[] = " ";
   
   if(!(fd = fopen(INFILE, "r")))
     error_open_file(INFILE);

   while(fgets(line, sizeof(line), fd) != NULL){
     int ii=0;
     char *result = NULL;
     char *var[N_COL_LC_ASCII];
     result = strtok( line, delims );
     while( result != NULL ) {
       var[ii]=result;
       result = strtok( NULL, delims );
       ii ++;
     }
     lc_cat.lc[n].ra      = atof(var[0]);
     lc_cat.lc[n].dec     = atof(var[1]);
     lc_cat.lc[n].rds     = atof(var[2]);
     lc_cat.lc[n].zreal   = atof(var[3]);
     lc_cat.lc[n].Mv      = atof(var[4]);
     lc_cat.lc[n].M200    = atof(var[5]);
     lc_cat.lc[n].Rv      = atof(var[6]);
     lc_cat.lc[n].Vmax    = atof(var[7]);
     lc_cat.lc[n].Vpeak   = atof(var[8]);
     lc_cat.lc[n].str_id  = atol(var[9]);
     lc_cat.lc[n].id      = atol(var[10]);
     lc_cat.lc[n].snapnum = atoi(var[11]);
     lc_cat.lc[n].boxnum  = atoi(var[12]);
     lc_cat.lc[n].lc_id   = n;
     lc_cat.lc[n].sel     = FALSE;
     n++;
   }
   
   fclose(fd);
}

void select_rotation(long *i,long *j,long *k,float *sig, float *disp){
  if(ROTBOX==1){*i=0; *j=1; *k=2; *sig=1.0; *disp=0.0;}
  if(ROTBOX==2){*i=1; *j=2; *k=0; *sig=1.0; *disp=0.0;}
  if(ROTBOX==3){*i=3; *j=0; *k=1; *sig=1.0; *disp=0.0;}
  if(ROTBOX==4){*i=0; *j=1; *k=2; *sig=-1.0; *disp=LBOX;}
  if(ROTBOX==5){*i=1; *j=2; *k=0; *sig=-1.0; *disp=LBOX;}
  if(ROTBOX==6){*i=3; *j=0; *k=1; *sig=-1.0; *disp=LBOX;}
}

void read_box_ascii(char *fname,box box_cat){

   FILE *fd;
   long  n=0;
   char line[2*CHAR_SIZE];
   char delims[] = " ";
   char *tmp;
   long i=0,j=0,k=0;
   float DIS=0.0,sig=0.0;
   
   select_rotation(&i,&j,&k,&sig,&DIS);

   fprintf(stderr,"%ld %s\n",box_cat.np,fname);
   if(!(fd = fopen(fname, "r")))
     error_open_file(fname);

   while(fgets(line, sizeof(line), fd) != NULL){
     int ii=0;
     char *result = NULL;
     char *var[N_COL_BOX_ASCII];
     result = strtok_r( line, delims,&tmp);
     while( result != NULL ) {
       var[ii]=result;
       result = strtok_r( NULL, delims,&tmp);
       ii ++;
     }
     box_cat.bx[n].id     = atol(var[0]);
     box_cat.bx[n].x      = DIS+sig*atof(var[i+1]);
     box_cat.bx[n].y      = DIS+sig*atof(var[j+1]);
     box_cat.bx[n].z      = DIS+sig*atof(var[k+1]);
     box_cat.bx[n].vx     = sig*atof(var[i+4]);
     box_cat.bx[n].vy     = sig*atof(var[j+4]);
     box_cat.bx[n].vz     = sig*atof(var[k+4]);
     box_cat.bx[n].Mv     = atof(var[7]);
     box_cat.bx[n].Rv     = atof(var[8]);
     box_cat.bx[n].M200   = atof(var[9]);
     box_cat.bx[n].str_id = atol(var[10]);
     box_cat.bx[n].Vpeak  = atof(var[11]);
     box_cat.bx[n].Vmax   = atof(var[12]);
     box_cat.bx[n].box_id = n;
     box_cat.bx[n].sel    = FALSE;
     /*box_cat.bx[n].x      = DIS+sig*atof(var[i]);
     box_cat.bx[n].y      = DIS+sig*atof(var[j]);
     box_cat.bx[n].z      = DIS+sig*atof(var[k]);
     box_cat.bx[n].vx     = sig*atof(var[i+3]);
     box_cat.bx[n].vy     = sig*atof(var[j+3]);
     box_cat.bx[n].vz     = sig*atof(var[k+3]);
     box_cat.bx[n].sel    = FALSE;*/
     n++;
   }
   
   fclose(fd);
}

void read_angular_random(ran_lc lc_ran){

   FILE *fd;
   long  n=0;
   char line[2*CHAR_SIZE];
   char delims[] = " ";
   
   if(!(fd = fopen(RANFILE, "r")))
     error_open_file(RANFILE);

   while(fgets(line, sizeof(line), fd) != NULL){
     int ii=0;
     char *result = NULL;
     char *var[5];
     result = strtok( line, delims );
     for(ii=0;ii<2;ii++){
       var[ii]=result;
       result = strtok( NULL, delims );
     }
     lc_ran.lc[n].ra      = atof(var[0]);
     lc_ran.lc[n].dec     = atof(var[1]);
     lc_ran.lc[n].sel     = FALSE;
     n++;
     if(n==lc_ran.np) break;
   }
   
   fclose(fd);
}

