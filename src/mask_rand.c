#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "hdf5.h"
#include "common.h"
#include "define.h"
#include "cosmo_tools.h"
#include <minimal_mangle.c>
#include <rng_gsl.c>

#define seed_ds 1357

void apply_mask_ply(lightcone lc_cat){
  MANGLE_PLY *ply;
  MANGLE_INT ipoly;
  int ii=0;
  struct drand48_data drand_buf;
  double weight = 1.0,ran1=0.0;

  ply = mply_read_file( MASKFILE );
#pragma omp parallel private(ipoly,weight,ran1,drand_buf) shared(ply) reduction(-:NSEL)
  {
    unsigned int seed=omp_get_thread_num()+(int)time(NULL);
    srand48_r(seed, &drand_buf);
    
#pragma omp for 
    for(ii=0; ii<lc_cat.np; ii++){
      if(lc_cat.lc[ii].sel){
	drand48_r(&drand_buf,&ran1);
	ipoly = mply_find_polyindex_radec(ply,lc_cat.lc[ii].ra,lc_cat.lc[ii].dec);
	if( ipoly < 0 ){
	  lc_cat.lc[ii].sel=FALSE;
	  NSEL--;
	  continue;}
	weight = mply_weight_from_index( ply, ipoly );
	if( weight < MIN_COMP){
	  lc_cat.lc[ii].sel=FALSE;
	  NSEL--;
	  continue;}
	if( ran1 > weight && ANG_DOWN){
	  lc_cat.lc[ii].sel=FALSE;
	  NSEL--;
	}
      }
    }
  }
}

void apply_mask_ran(ran_lc lc_cat){
  MANGLE_PLY *ply;
  MANGLE_INT ipoly;
  int ii=0;
  struct drand48_data drand_buf;
  double weight = 1.0,ran1=0.0;
  void *rng;

  ply = mply_read_file( MASKFILE );
#pragma omp parallel private(ipoly,weight,ran1,drand_buf) shared(ply) reduction(+:NS_RAN)
  {
    unsigned int seed=omp_get_thread_num()+(int)time(NULL);
    srand48_r(seed, &drand_buf);
    
#pragma omp for 
    for(ii=0; ii<lc_cat.np; ii++){
      ran1 = rng_uniform( rng );
      ipoly = mply_find_polyindex_radec(ply,lc_cat.lc[ii].ra,lc_cat.lc[ii].dec);
      if( ipoly < 0 )
	continue;
      weight = mply_weight_from_index( ply, ipoly );
      if( weight < MIN_COMP)
	continue;
      if( ran1 > weight && ANG_DOWN)
	continue;
      lc_cat.lc[ii].sel=TRUE;
      NS_RAN++;
    }
  }
}
  
void angular_random_inmask(ran_lc lc_ran){
  long ii=0;
  MANGLE_PLY *ply;
  MANGLE_INT ipoly;
  struct drand48_data drand_buf;
  double weight = 1.0,ran=0.0,ran1=0.0,ran2=0.0;
  float RA_MIN=RA_LOS-RA_SIZE*0.5;
  float RA_MAX=RA_LOS+RA_SIZE*0.5;
  float CDEC_MIN=cos(M_PI*0.5+(DEC_LOS-DEC_SIZE*0.5)*D2R);
  float CDEC_MAX=cos(M_PI*0.5+(DEC_LOS+DEC_SIZE*0.5)*D2R);

  ply = mply_read_file( MASKFILE );
  
#pragma omp parallel private(ipoly,weight,ran,ran1,ran2,drand_buf) shared(ply)
  {
    unsigned int seed=omp_get_thread_num()+(int)time(NULL);
    float ra,dec;
    srand48_r(seed, &drand_buf);
    
#pragma omp for nowait schedule(static)
    for(ii=0; ii<lc_ran.np; ii++){
      int FLAG=FALSE;
      while(!FLAG){
	drand48_r(&drand_buf,&ran1);
	drand48_r(&drand_buf,&ran2);
	ra=ran1*(RA_MAX-RA_MIN)+RA_MIN;
	if(ra<0.0)ra=ra+360.0;
	dec=R2D*(acos(ran2*(CDEC_MAX-CDEC_MIN)+CDEC_MIN)-M_PI*0.5);
	drand48_r(&drand_buf,&ran);
	ipoly = mply_find_polyindex_radec(ply,ra,dec);
	if( ipoly < 0 )
	  continue;
	weight = mply_weight_from_index( ply, ipoly );
	if( weight < MIN_COMP)
	  continue;
	if( ran > weight && ANG_DOWN)
	  continue;
	FLAG=TRUE;
      }
      lc_ran.lc[ii].ra=ra;
      lc_ran.lc[ii].dec=dec;
      lc_ran.lc[ii].sel=TRUE;
    }
  }
}

void angular_random(ran_lc lc_ran){
  long ii=0;
  float ra=0.0,dec=0.0;
  double ran1=0.0,ran2=0.0;
  struct drand48_data drand_buf;
  float RA_MIN=RA_LOS-RA_SIZE*0.5;
  float RA_MAX=RA_LOS+RA_SIZE*0.5;
  float CDEC_MIN=cos(M_PI*0.5+(DEC_LOS-DEC_SIZE*0.5)*D2R);
  float CDEC_MAX=cos(M_PI*0.5+(DEC_LOS+DEC_SIZE*0.5)*D2R);
  
#pragma omp parallel private(ra,dec,ran1,ran2,drand_buf)
  {
    unsigned int seed=omp_get_thread_num()+(int)time(NULL);
    srand48_r(seed, &drand_buf);
    
#pragma omp for 
    for(ii=0; ii<lc_ran.np; ii++){
      drand48_r(&drand_buf,&ran1);
      drand48_r(&drand_buf,&ran2);
      ra=ran1*(RA_MAX-RA_MIN)+RA_MIN;
      if(ra<0.0)ra=ra+360.0;
      dec=R2D*(acos(ran2*(CDEC_MAX-CDEC_MIN)+CDEC_MIN)-M_PI*0.5);
      lc_ran.lc[ii].ra=ra;
      lc_ran.lc[ii].dec=dec;
      lc_ran.lc[ii].sel=TRUE;
    }
  }
}


double generate_redshift(int ii){

  struct drand48_data drand_buf;
  unsigned int seed=(int)time(NULL)+ii;
  srand48_r(seed, &drand_buf);
  double rds=0;
  double xran=0.0,yran=0.0,zran=0.0;
  double X=0.0,Y=0.0,Z=0.0,R=-1.0;
  while(R<RC_MIN || R>=RC_MAX){
    drand48_r(&drand_buf,&xran);
    drand48_r(&drand_buf,&yran);
    drand48_r(&drand_buf,&zran);
    X=xran*RC_MAX;
    Y=yran*RC_MAX;
    Z=zran*RC_MAX;
    R=sqrt(X*X+Y*Y+Z*Z);
  }
  rds=dist2redshift(R);
  return rds;
}

void assign_redshift_randomcat(lightcone lc_cat, ran_lc r_cat){

  if(!SHUFFLE){
    struct drand48_data drand_buf;
    long ii=0;
    RC_MAX=HUBBLE_PAR*com_distance(MAXZ);
    RC_MIN=HUBBLE_PAR*com_distance(MINZ);
    read_radial_nbar(FALSE,NULL);

#pragma omp parallel
    {
      unsigned int seed=omp_get_thread_num()+(int)time(NULL);
      srand48_r(seed, &drand_buf);
      
#pragma omp for 
      for(ii=0;ii<r_cat.np;ii++){
	double ran=0.0;
	drand48_r(&drand_buf,&ran);
	float z=generate_redshift((int)(ran*1e6));
	while(ran > number_density(z)/MAX_NBAR){
	  drand48_r(&drand_buf,&ran);
	  z=generate_redshift((int)(ran*1e6));
	}
	r_cat.lc[ii].rds = z;
	r_cat.lc[ii].wfkp = 1.0/(1.0+P0*number_density(z));
      }
    }
    free_gsl_nbar();
  }
  if(SHUFFLE){
    long ii=0,jj=0;
    int kk=0;

    if(NO_SEL){
#pragma omp parallel for private(ii,jj)
      for(kk=0; kk<FRAN; kk++){
	ii=0;
	for(jj=0; jj<lc_cat.np;jj++){
	  r_cat.lc[ii+kk*NSEL].rds=lc_cat.lc[jj].rds;
	  r_cat.lc[ii+kk*NSEL].zreal=lc_cat.lc[jj].zreal;
	  r_cat.lc[ii+kk*NSEL].wfkp=lc_cat.lc[jj].wfkp;
	  ii++;
	}
      }
    }
    else{
#pragma omp parallel for private(ii,jj)
      for(kk=0; kk<FRAN; kk++){
	ii=0;
	for(jj=0; jj<lc_cat.np;jj++)
	  if(lc_cat.lc[jj].sel){
	    r_cat.lc[ii+kk*NSEL].rds=lc_cat.lc[jj].rds;
	    r_cat.lc[ii+kk*NSEL].zreal=lc_cat.lc[jj].zreal;
	    r_cat.lc[ii+kk*NSEL].wfkp=lc_cat.lc[jj].wfkp;
	    ii++;
	  }
      }
    }
  }
}

void write_random_ascii(ran_lc dat){
  FILE *fd;
  int  ii=0;
  char fname[CHAR_SIZE];

  sprintf(fname,"%s%s%.2f%s%.2f%s",OUTFILE,".minz.",MINZ,".maxz.",MAXZ,".ran");
  fd = fopen(fname,"w");
  for(ii=0; ii<dat.np; ii++)
    /*fprintf(fd,"%.6f %.6f %.6f %.6f\n", dat.lc[ii].ra,dat.lc[ii].dec,
      dat.lc[ii].rds,dat.lc[ii].wfkp);*/
    fprintf(fd,"%.6f %.6f %.6f %.6f %.6f\n", dat.lc[ii].ra,dat.lc[ii].dec,
	    dat.lc[ii].rds,dat.lc[ii].wfkp,dat.lc[ii].zreal);
  fclose(fd);  
}

void fill_random(ran_lc aux_cat,ran_lc lc_cat){
  long ii=0,jj=0;
  for(ii=0;ii<aux_cat.np;ii++){
    if(aux_cat.lc[ii].sel){
      lc_cat.lc[jj].ra  = aux_cat.lc[ii].ra;
      lc_cat.lc[jj].dec = aux_cat.lc[ii].dec;
      jj++;
    }
    if(jj==lc_cat.np)
      break;
  }
}

void construct_random(ran_lc *lc_ran,lightcone lc_cat){
  int NRANDOM=0;
  if(!strcmp(RANFILE, "UNDEFINED")){
    fprintf(stderr,"Building angular and radial random ...");
    NRANDOM=FRAN*NSEL;
    init_lightcone_rand(lc_ran);
    allocate_mem_random(NRANDOM,lc_ran);
    if(MASK)
      angular_random_inmask((*lc_ran));
    if(!MASK)
      angular_random((*lc_ran));
  }
  else{
    long nran=n_halos_ascii(RANFILE);
    fprintf(stderr,"Reading angular randoms ... ");
    if(MASK){
      ran_lc aux_ran;
      init_lightcone_rand(&aux_ran);
      allocate_mem_random(nran,&aux_ran);
      read_angular_random(aux_ran);
      apply_mask_ran(aux_ran);
      if(NS_RAN < FRAN*NSEL){
	FRAN=(int)(NS_RAN/NSEL);
	fprintf(stderr,"Small number of random in file %s\n"
		"FRAN modified to %d\n"
		"(WARNING)\n\n",RANFILE,FRAN);
      }
      NRANDOM=FRAN*NSEL;
      init_lightcone_rand(lc_ran);
      allocate_mem_random(NRANDOM,lc_ran);
      fill_random(aux_ran,(*lc_ran));
    }
    if(!MASK){
      NRANDOM=FRAN*NSEL;

      if(nran < FRAN*NSEL){
	FRAN=(int)(nran/NSEL);
	fprintf(stderr,"Small number of random in file %s\n"
		"FRAN modified to %d\n"
		"(WARNING)\n\n",RANFILE,FRAN);
      }
      init_lightcone_rand(lc_ran);
      allocate_mem_random(NRANDOM,lc_ran);
      read_angular_random((*lc_ran));
    }
    timer(1);
    fprintf(stderr,"Assigning redshifts to randoms ... ");
  }

  assign_redshift_randomcat(lc_cat, (*lc_ran));
  timer(1);
  fprintf(stderr,"Writing random catalogue...");
  write_random_ascii((*lc_ran));
  timer(1);
}	

