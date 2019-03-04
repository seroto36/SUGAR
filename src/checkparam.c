#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"
#include "define.h"


void error_open_file(char *fname){
  fprintf(stderr,"Could not open file: \"%s\"\n"
	  "(PROCESS ABORTED)\n\n",fname);
  exit(EXIT_FAILURE);
}

void error_line_cfile(int nline){
  fprintf(stderr,"Could not read line %d in file \"%s\"\n" 
	  "  Do not include space in the beginnig of the line\n"
	  "  Line format: NAME_VARIBLE yes/no (#comment)\n"
	  "(PROCESS ABORTED)\n\n",nline,fconf);
  exit(EXIT_FAILURE);
}

void error_line_numvar(int nline, char *file){
  fprintf(stderr,"Could not read line %d in file \"%s\"\n" 
	  "  wrong number of columns\n"
	  "  number of columns = 2\n"
	  "(PROCESS ABORTED)\n\n",nline,file);
  exit(EXIT_FAILURE);
}

void error_input_param(int nline,char *name,char *fname){
  fprintf(stderr,"Unknown parameter \"%s\" \n"
	  "file: \"%s\" line %d\n"
	  "(PROCESS ABORTED)\n\n",name,fname,nline);
  exit(EXIT_FAILURE);
}

void error_missing_param(char *fname, char *param){
  fprintf(stderr,"Missing param \"%s\" in file \"%s\"\n"
	  "(PROCESS ABORTED)\n\n",param,fname);
  exit(EXIT_FAILURE);
}

void error_repeat_param(int nline, char *fname, char *param){
  fprintf(stderr,"Repeat param \"%s\" in line %d. File \"%s\"\n"
	  "(PROCESS ABORTED)\n\n",param,nline,fname);
  exit(EXIT_FAILURE);
}

void error_malloc_allocation(char *var){
  fprintf(stderr,"Varible %s unallocated\n"
	  "(PROCESS ABORTED)\n\n",var);
  exit(EXIT_FAILURE);;
}

void check_nbar_file(){

  FILE *fd;
  char buffer[CHAR_SIZE];
  int ii=0,nP=0;
  
  if(!(fd = fopen(NBARFILE, "r")))
    error_open_file(NBARFILE);

  while(fgets(buffer,CHAR_SIZE,fd)!=NULL)nP++;
  rewind(fd);
  double rds[nP],nbar[nP];
  
  for(ii=0; ii<nP; ii++){
    fscanf(fd,"%lf %lf",&rds[ii],&nbar[ii]);
    if(ii>0 && rds[ii]<rds[ii-1]){
      fprintf(stderr,"Error in n(z) file:\n"
	      "Redshifts have to be in ascending order\n"
	      "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  } 
  
  fclose(fd);
}


void check_conf_lc(){

  if(HAM && HOD){
    printf("HAM and HOD cannot be used together\n\n");
    exit(EXIT_FAILURE);
  }
  if(HAM && MODIFIED_HAM){
    printf("HAM and MODIFIED_HAM cannot be used together\n"
	   "(PROCESS ABORTED)\n\n");
    exit(EXIT_FAILURE);
  }
  if(HAM && DOWNSAMPLE){
    printf("HAM and DOWNSAMPLE cannot be used together\n"
	   "(PROCESS ABORTED)\n\n");
    exit(EXIT_FAILURE);
  }
  if(HOD && MODIFIED_HAM){
    printf("HOD and MODIFIED_HAM cannot be used together\n"
	   "(PROCESS ABORTED)\n\n");
    exit(EXIT_FAILURE);
  }
  if(HOD && DOWNSAMPLE){
    printf("HOD and DOWNSAMPLE cannot be used together\n"
	   "(PROCESS ABORTED)\n\n");
    exit(EXIT_FAILURE);
  }
  if(MODIFIED_HAM && DOWNSAMPLE){
    printf("MODIFIED_HAM and DOWNSAMPLE cannot be used together\n"
	   "(PROCESS ABORTED)\n\n");
    exit(EXIT_FAILURE);
  }
  if(VAR_CUT && (DOWNSAMPLE || MODIFIED_HAM)){
    printf("MODIFIED_HAM and DOWNSAMPLE cannot be used with VAR_CUT.\n"
	   "It is necessary to include a n(z) file\n"
	   "(PROCESS ABORTED)\n\n");
    exit(EXIT_FAILURE);
  }
  if(PART_CAT && (HAM || HOD || MODIFIED_HAM)){
    printf("PART_CAT does not work with selec. methods: HAM, HOD, MODIFIED_HAM\n"
	   "(PROCESS ABORTED)\n\n");
    exit(EXIT_FAILURE);
  }
  if(!HAM && !HOD && !MODIFIED_HAM && !DOWNSAMPLE && !VAR_CUT)
    NO_SEL=TRUE;
  if(NO_SEL && READ_LC && !MASK){
    printf("Any action set to be applied on the input lightcone\n"
	   "(PROCESS ABORTED)\n\n");
    exit(EXIT_FAILURE);
  }
}

void check_conf_bx(){
  if(HAM && HOD){
    printf("HAM and HOD cannot be used together\n"
	   "(PROCESS ABORTED)\n\n");
    exit(EXIT_FAILURE);
  }
  if(HAM && MODIFIED_HAM){
    printf("HAM and MODIFIED_HAM cannot be used together\n"
	   "(PROCESS ABORTED)\n\n");
    exit(EXIT_FAILURE);
  }
  if(HAM && DOWNSAMPLE){
    printf("HAM and DOWNSAMPLE cannot be used together\n"
	   "(PROCESS ABORTED)\n\n");
    exit(EXIT_FAILURE);
  }
  if(HOD && MODIFIED_HAM){
    printf("HOD and MODIFIED_HAM cannot be used together\n\n");
    exit(EXIT_FAILURE);
  }
  if(HOD && DOWNSAMPLE){
    printf("HOD and DOWNSAMPLE cannot be used together\n"
	   "(PROCESS ABORTED)\n\n");
    exit(EXIT_FAILURE);
  }
  if(MODIFIED_HAM && DOWNSAMPLE){
    printf("MODIFIED_HAM and DOWNSAMPLE cannot be used together\n"
	   "(PROCESS ABORTED)\n\n");
    exit(EXIT_FAILURE);
  }
  if(PART_CAT==TRUE && (HAM==TRUE || HOD==TRUE || MODIFIED_HAM==TRUE)){
    printf("PART_CAT does not work with selec. methods: HAM, HOD, MODIFIED_HAM\n"
	   "(PROCESS ABORTED)\n\n");
    exit(EXIT_FAILURE);
  }
  if(VAR_CUT==TRUE && (DOWNSAMPLE==TRUE || MODIFIED_HAM==TRUE)){
    printf("MODIFIED_HAM and DOWNSAMPLE cannot be used with VAR_CUT.\n"
	   "It is necessary to provide a constant number density in param file\n"
	   "(PROCESS ABORTED)\n\n");
    exit(EXIT_FAILURE);
  }
}

void check_params(){

  /*Mandatory parametters*/
  if(OMEGA_M==DEFAULT){
    fprintf(stderr,"OMEGA_M: Parameter not found in param file\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  if(OMEGA_L==DEFAULT){
    fprintf(stderr,"OMEGA_L: Parameter not found in param file\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  if(HUBBLE_PAR==DEFAULT){
    fprintf(stderr,"HUBBLE_PAR: Parameter not found in param file\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  if(MINZ==DEFAULT && !strcmp(CODE_OPT, "LGC")){
    fprintf(stderr,"MINZ: Parameter not found in param file\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  if(MAXZ==DEFAULT && !strcmp(CODE_OPT, "LGC")){
    fprintf(stderr,"MAXZ: Parameter not found in param file\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  if(((!strcmp(CODE_OPT, "LGC") && !READ_LC) || !strcmp(CODE_OPT, "BOX"))
     && LBOX==DEFAULT){
    fprintf(stderr,"LBOX: Parameter not found in param file\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  if(!strcmp(CODE_OPT, "LGC") && NO_SEL==DEFAULT && NPART<=0){
    fprintf(stderr,"N_PART: Parameter not found or negative value\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  if(ZBOX==DEFAULT && !strcmp(CODE_OPT, "BOX")){
    fprintf(stderr,"ZBOX: Parameter not found in param file\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  if(NBAR==DEFAULT && !strcmp(CODE_OPT, "BOX")){
    fprintf(stderr,"NBAR: Parameter not found in param file\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  if(FSAT!=DEFAULT && HAM)
    fprintf(stderr,"Using FSAR = %.2f with HAM mode\n",FSAT);	    

  
  /*Optional parameters*/
  if(READ_LC && !strcmp(INFILE, "UNDEFINED")){
    fprintf(stderr,"INFILE: Parameter not found in param file\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  if(!strcmp(INFILE, "UNDEFINED") && !strcmp(CODE_OPT, "BOX")){
    fprintf(stderr,"INFILE: Parameter not found in param file\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  if(!READ_LC && !strcmp(CODE_OPT, "LGC") &&
     !strcmp(F_RDS_SIM,"UNDEFINED")){
       fprintf(stderr,"F_RDS_SIM: Parameter not found in param file\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  if(MASK && !strcmp(CODE_OPT, "LGC") &&!strcmp(MASKFILE,"UNDEFINED")){
       fprintf(stderr,"MASKFILE: Parameter not found in param file\n"
	       "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  if(!READ_LC && !strcmp(CODE_OPT, "LGC") &&
     !strcmp(BASEFILE,"UNDEFINED")){
       fprintf(stderr,"BASEFILE: Parameter not found in param file\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  if(!READ_LC && !strcmp(CODE_OPT, "LGC") &&
     !strcmp(ENDFILE,"UNDEFINED")){
       fprintf(stderr,"ENDFILE: Parameter not found in param file\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  if(!READ_LC && !strcmp(CODE_OPT, "LGC") && (RA_LOS<0.0 || RA_LOS>=360.0)){
       fprintf(stderr,"RA_LOS: Parameter not found or out of range [0,360)\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
  if(!READ_LC && !strcmp(CODE_OPT, "LGC") && (DEC_LOS<-90 || DEC_LOS>90)){
       fprintf(stderr,"DEC_LOS: Parameter not found or out of range (-90,90)\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
    if(!READ_LC && !strcmp(CODE_OPT, "LGC") && RA_SIZE<0){
       fprintf(stderr,"RA_SIZE: Parameter not found in param file\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
    if(!READ_LC && !strcmp(CODE_OPT, "LGC") && DEC_SIZE<0){
       fprintf(stderr,"DEC_SIZE: Parameter not found in param file\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
    if(!READ_LC && !strcmp(CODE_OPT, "LGC") && !strcmp(CONF_BOX, "UNDEFINED")){
      fprintf(stderr,"CONF_BOX: Parameter not found in param file\n"
	      "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
    if(!READ_LC && !strcmp(CODE_OPT, "LGC") && strcmp(CONF_BOX, "SDSS_REMAP")
       && strcmp(CONF_BOX, "NO_REMAP") ){
       fprintf(stderr,"CONF_BOX: Error in parameter -> SDSS_REMAP or NO_REMAP\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
    if(!READ_LC && !strcmp(CODE_OPT, "LGC") &&
       (X_OBS==DEF_COOR || Y_OBS==DEF_COOR || Z_OBS==DEF_COOR) ){
       fprintf(stderr,"_OBS: Observer coordinates error, check X_OBS,Y_OBS,Z_OBS\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
    if(!strcmp(OUTFILE,"UNDEFINED")){
      fprintf(stderr,"OUTFILE: Parameter not found in param file\n"
	    "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
    if(!strcmp(CODE_OPT, "LGC") && !strcmp(NBARFILE, "UNDEFINED") &&
       (DOWNSAMPLE  || HAM)){
      fprintf(stderr,"NBARFILE: Parameter not found in param file\n"
	      "Necessary when DOWNSAMPLE,HAM,MODIFIED_HAM or HOD are selected\n"
	      "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
    if(!strcmp(CODE_OPT, "LGC") && (DOWNSAMPLE || HAM))
      check_nbar_file();
    if(!strcmp(CODE_OPT, "LGC") && READ_LC && AREA==DEFAULT){
      fprintf(stderr,"AREA: Parameter necessary for selection\n"
	      "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
    if(ROTBOX<0 || ROTBOX>6){
      fprintf(stderr,"ROTBOX: takes integer values from 1 to 6\n"
	      "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
    if(!strcmp(CODE_OPT, "LGC") && RANDOM && !strcmp(RANFILE, "UNDEFINED") &&
       (RA_LOS==DEF_COOR || DEC_LOS==DEF_COOR || !RA_SIZE || !DEC_SIZE)){
      fprintf(stderr,"RA_SIZE,DEC_SIZE,RA_LOS,DEC_LOS: Missing for randoms\n"
	      "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
    if(HAM && SCAT_HAM==DEFAULT){
      fprintf(stderr,"SCAT_HAM: Parameter not found in param file\n"
	      "Necessary when HAM is set\n"
	      "(PROCESS ABORTED)\n\n"); exit(EXIT_FAILURE);}
	      
  /*Warnings*/
  if(!READ_LC && !strcmp(CODE_OPT, "LGC") &&
     strcmp(INFILE,"UNDEFINED")){
         fprintf(stderr,"INFILE is not used with the current set-up\n"
		 "(WARNING)\n\n");}
  if((READ_LC || !strcmp(CODE_OPT, "BOX")) &&
     strcmp(BASEFILE,"UNDEFINED")){
         fprintf(stderr,"BASEFILE is not used with the current set-up\n"
		 "(WARNING)\n\n");}
  if((READ_LC || !strcmp(CODE_OPT, "BOX")) &&
     strcmp(ENDFILE,"UNDEFINED")){
         fprintf(stderr,"ENDFILE is not used with the current set-up\n"
		 "(WARNING)\n\n");}
  if(!MASK && (ANG_DOWN || MIN_COMP!=FALSE)){
    fprintf(stderr,"ANG_DOWN and MIN_COMP are not used with the current set-up\n"
	    "(WARNING)\n\n");}
}
