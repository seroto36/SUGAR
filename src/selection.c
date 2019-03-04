#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <time.h>  
#include "cosmo_tools.h"
#include "common.h"
#include "define.h"

void init_sheels(sheels data[],lightcone lc_cat){
  int ii=0;
  float bin_size=(MAXZ-MINZ)/(float)NPART;
    
  for(ii=0; ii<NPART; ii++){
    data[ii].np=0;
    data[ii].z_min=MINZ+(float)ii*bin_size;
    data[ii].z_max=MINZ+(float)(ii+1)*bin_size;
    data[ii].r_min=com_distance(data[ii].z_min);
    data[ii].r_max=com_distance(data[ii].z_max);
    data[ii].vol=Volume(data[ii].r_min,data[ii].r_max);
  }
  
  for(ii=0; ii<lc_cat.np; ii++){
  if(lc_cat.lc[ii].rds >= MAXZ || lc_cat.lc[ii].rds < MINZ) continue;
    else{
      int itemp=(int)((lc_cat.lc[ii].rds-MINZ)/bin_size);
      data[itemp].np++;
    }
  }

  for(ii=0; ii<NPART; ii++)
    data[ii].nbar=(float)data[ii].np/data[ii].vol;
}
  
void down_sampling_selection(lightcone lc_cat){
  int ii=0;
  double ran_nd;
  struct drand48_data drand_buf;
  unsigned int seed=(int)time(NULL);
  srand48_r(seed, &drand_buf);

  for(ii=0; ii<lc_cat.np; ii++){
    if(lc_cat.lc[ii].rds >= MAXZ || lc_cat.lc[ii].rds < MINZ)
      continue;
    float Nnbar=norm_number_density(lc_cat.lc[ii].rds);
    drand48_r(&drand_buf,&ran_nd);
    if(ran_nd <= Nnbar){
      float zobs = lc_cat.lc[ii].rds;
      lc_cat.lc[ii].sel = TRUE;
      lc_cat.lc[ii].wfkp=1.0/(1.0+P0*number_density(zobs));
      NSEL++;
    }
  }
}

void downsampling_box(box box_cat){
  int ii=0;
  float Vol=(float)(LBOX*LBOX*LBOX);
  float BOXNBAR=(float)box_cat.np/Vol;
  float H = 100.0*(sqrt((1.0+ZBOX)*(1.0+ZBOX)*(1.0+ZBOX)*OMEGA_M+OMEGA_L));
  float DeltaZ = (1.0+ZBOX)/H;
  struct drand48_data drand_buf;

#pragma omp parallel
  {
    unsigned int seed=omp_get_thread_num()+(int)time(NULL);
    double ran_nd=0.0;
    srand48_r(seed, &drand_buf);
    
#pragma omp for nowait schedule(static)
    for(ii=0; ii<box_cat.np; ii++){
      drand48_r(&drand_buf,&ran_nd);
      if(ran_nd > NBAR/BOXNBAR)
	continue;
      box_cat.bx[ii].sel  = TRUE;
      box_cat.bx[ii].zrds = box_cat.bx[ii].z+box_cat.bx[ii].vz*DeltaZ;
      if(box_cat.bx[ii].zrds >= LBOX)
	box_cat.bx[ii].zrds-=LBOX;
      if(box_cat.bx[ii].zrds < 0.0)
	box_cat.bx[ii].zrds+=LBOX;
    }
  }
}


void downsampling(lightcone lc_cat){
  sheels parts[NPART];
  init_sheels(parts,lc_cat);
  read_radial_nbar(TRUE,parts);
  down_sampling_selection(lc_cat);
  free_gsl_NORMnbar();
}

void select_ham_obj(hl_proxy *ham_px, lightcone lc_cat){
  long ii=0,jj=0,kk=0;
  if(FSAT==DEFAULT)
    for(ii=0; ii<NPART; ii++){
      for(jj=0; jj<ham_px[ii].nsel;jj++){
	kk=ham_px[ii].prx[jj].pst;
	float zobs = lc_cat.lc[kk].rds;
	lc_cat.lc[kk].sel = TRUE;
	lc_cat.lc[kk].wfkp=1.0/(1.0+P0*number_density(zobs));
	NSEL++;
      }}

  else
    for(ii=0; ii<NPART; ii++){
      int NFRAC=(int)(FSAT*ham_px[ii].nsel);
      for(jj=0; jj<ham_px[ii].np;jj++){
	kk=ham_px[ii].prx[jj].pst;
	if(lc_cat.lc[kk].str_id!=-1){
	  float zobs = lc_cat.lc[kk].rds;
	  lc_cat.lc[kk].sel = TRUE;
	  lc_cat.lc[kk].wfkp=1.0/(1.0+P0*number_density(zobs));
	  NSEL++;
	}
	if(NSEL==NFRAC)
	  break;
      }

      for(jj=0; jj<ham_px[ii].np;jj++){
	kk=ham_px[ii].prx[jj].pst;
	if(lc_cat.lc[kk].str_id==-1){
	  float zobs = lc_cat.lc[ham_px[ii].prx[jj].pst].rds;
	  lc_cat.lc[kk].sel = TRUE;
	  lc_cat.lc[kk].wfkp=1.0/(1.0+P0*number_density(zobs));
	  NSEL++;
	}
	if(NSEL==ham_px[ii].nsel)
	  break;
      }      
    }
}

void select_ham_box(proxy *ham_px, box box_cat){
  long jj=0,ii=0;
  long NOBJ=(long)(NBAR*LBOX*LBOX*LBOX);
  float H = 100.0*(sqrt((1.0+ZBOX)*(1.0+ZBOX)*(1.0+ZBOX)*OMEGA_M+OMEGA_L));
  float DeltaZ = (1.0+ZBOX)/H;

  if(FSAT==DEFAULT)
    for(jj=0; jj<NOBJ;jj++){
      ii=ham_px[jj].pst;
      box_cat.bx[ii].zrds = box_cat.bx[ii].z+box_cat.bx[ii].vz*DeltaZ;
      box_cat.bx[ii].sel = TRUE;
      if(box_cat.bx[ii].zrds >= LBOX)
	box_cat.bx[ii].zrds-=LBOX;
      if(box_cat.bx[ii].zrds < 0.0)
	box_cat.bx[ii].zrds+=LBOX;
      NSEL++;
    }

  else{
    int NFRAC=(int)(FSAT*NOBJ);
    for(jj=0; jj<box_cat.np;jj++){
      ii=ham_px[jj].pst;
      if(box_cat.bx[ii].str_id!=-1){
	box_cat.bx[ii].zrds = box_cat.bx[ii].z+box_cat.bx[ii].vz*DeltaZ;
	box_cat.bx[ii].sel = TRUE;
	if(box_cat.bx[ii].zrds >= LBOX)
	  box_cat.bx[ii].zrds-=LBOX;
	if(box_cat.bx[ii].zrds < 0.0)
	  box_cat.bx[ii].zrds+=LBOX;
	NSEL++;}
      if(NSEL==NFRAC)
	break;
    }
    for(jj=0; jj<box_cat.np;jj++){
      ii=ham_px[jj].pst;
      if(box_cat.bx[ii].str_id==-1){      
	box_cat.bx[ii].zrds = box_cat.bx[ii].z+box_cat.bx[ii].vz*DeltaZ;
	box_cat.bx[ii].sel = TRUE;
	if(box_cat.bx[ii].zrds >= LBOX)
	  box_cat.bx[ii].zrds-=LBOX;
	if(box_cat.bx[ii].zrds < 0.0)
	  box_cat.bx[ii].zrds+=LBOX;
	NSEL++;}
      if(NSEL==NOBJ)
	break;
    }
  }
}

float stellar_mass(float volume,int count){
  double steMASS=0.0;
  double h3=HUBBLE_PAR*HUBBLE_PAR*HUBBLE_PAR;
  double cou_per_vol=(double)count*h3/(double)volume;

  gsl_spline_eval_e(spline_CMF,cou_per_vol,Acc_CMF,&steMASS);

  return steMASS;
}

void select_ham_obj_sm(hl_proxy *ham_px, lightcone lc_cat){
  long ii=0,jj=0;
  for(ii=0; ii<NPART; ii++){
    for(jj=0; jj<ham_px[ii].nsel;jj++){
      float zobs = lc_cat.lc[ham_px[ii].prx[jj].pst].rds;
      float SM = stellar_mass(ham_px[ii].vol,jj+1);
      lc_cat.lc[ham_px[ii].prx[jj].pst].sel = TRUE;
      lc_cat.lc[ham_px[ii].prx[jj].pst].Mstar = SM;
      lc_cat.lc[ham_px[ii].prx[jj].pst].wfkp=1.0/(1.0+P0*number_density(zobs));

      NSEL++;
    }
  }
}


void ProbDistr_SM(hl_proxy *ham_px, lightcone lc_cat, mass_sel *smbins){

  long ii=0,jj=0,ll=0;
  float INV_MBIN_SIZE = (float)MASS_NBIN/(MAX_MASS-MIN_MASS);
  
  for(ii=0; ii<NPART; ii++){
    for(ll=0; ll<ham_px[ii].np;ll++){
      jj=ham_px[ii].prx[ll].pst;
      float SM = stellar_mass(ham_px[ii].vol,ll+1);
      int smp=(int)((SM-MIN_MASS)*INV_MBIN_SIZE);
      if(SM>=MAX_MASS || SM<MIN_MASS) continue;
      smbins[ii].i_obj[smp]++;
      lc_cat.lc[jj].Mstar = SM;
    }    
    
    float rest_obj = 0.0;
    for(ll=MASS_NBIN-1; ll>=0; ll--){      
      if(smbins[ii].f_obj[ll] > (float)smbins[ii].i_obj[ll]){
	rest_obj = rest_obj+smbins[ii].f_obj[ll] - (float)smbins[ii].i_obj[ll];
	smbins[ii].f_obj[ll] = (float)smbins[ii].i_obj[ll]; 
      }      
      else{
	smbins[ii].f_obj[ll] = smbins[ii].f_obj[ll] + rest_obj;
	rest_obj = 0.0;
	if(smbins[ii].f_obj[ll] > (float)smbins[ii].i_obj[ll]){
	  rest_obj = smbins[ii].f_obj[ll] - (float)smbins[ii].i_obj[ll];
	  smbins[ii].f_obj[ll] = (float)smbins[ii].i_obj[ll]; 
	}
      }
      if(smbins[ii].i_obj[ll]>0)
	smbins[ii].prob[ll]=smbins[ii].f_obj[ll]/(float)smbins[ii].i_obj[ll];
      else
	smbins[ii].prob[ll]=0;
    }
  }
}

void select_mass(hl_proxy *ham_px, lightcone lc_cat, mass_sel *smbins){

  int ii=0,ll=0,jj=0;
  double rand=0.0;
  float INV_MBIN_SIZE = (float)MASS_NBIN/(MAX_MASS-MIN_MASS);
  struct drand48_data drand_buf;
  unsigned int seed=(int)time(NULL);
  srand48_r(seed, &drand_buf);

  for(ii=0; ii<NPART; ii++){
    for(ll=0; ll<ham_px[ii].np; ll++){
      jj=ham_px[ii].prx[ll].pst;
      int smp=(int)((lc_cat.lc[jj].Mstar-MIN_MASS)*INV_MBIN_SIZE);
      if(smbins[ii].prob[smp]==0.0)
	continue;
      drand48_r(&drand_buf,&rand);
      if(rand>smbins[ii].prob[smp])
	continue;
      float zobs = lc_cat.lc[jj].rds;
      lc_cat.lc[jj].sel=TRUE;
      lc_cat.lc[jj].wfkp=1.0/(1.0+P0*number_density(zobs));
      NSEL++;
    }	
  } 
}

void halo_abundance(lightcone lc_cat){
  hl_proxy *ham_px;
  init_sheels_ham(&ham_px);
  allocate_sheels_ham(&ham_px,lc_cat);
  fill_sheels(ham_px,lc_cat);
  if(!CSMF)
    select_ham_obj(ham_px,lc_cat);
  if(CSMF){
    CSMF_calculator();
    if(!strcmp(SMFCAT, "UNDEFINED"))
      select_ham_obj_sm(ham_px,lc_cat);
    if(strcmp(SMFCAT, "UNDEFINED")){
      mass_sel *smbins;
      init_mass_sel(&smbins);
      ReadandClass_SM(smbins,ham_px);
      ProbDistr_SM(ham_px,lc_cat,smbins);
      select_mass(ham_px,lc_cat,smbins);
    }
    free_gsl_CMF();
  }
}

void halo_abundance_box(box box_cat){
  proxy *ham_px;
  long size = box_cat.np;
  ham_px  = (proxy *)malloc(size*sizeof(proxy));

  fill_proxy_box(ham_px,box_cat);

  if(!CSMF)
    select_ham_box(ham_px,box_cat);
  if(CSMF){
    fprintf(stderr,"Under construction!\n");
    exit(EXIT_FAILURE);
  }
}

void lightcone_selection(lightcone lc_cat){
  //Print lightcone with all halos within the geometry
  if(NO_SEL){
    fprintf(stderr,"Writting output ...");
    NSEL=lc_cat.np;
    if(WRITE_HDF5)
      write_hdf5(lc_cat);
    else
      write_ascii(lc_cat);
    timer(1);
  }
  else{
    /*Construct lightcone making a downsampling to reproduce
      the radial number density in nbar*/
    if(DOWNSAMPLE){
      fprintf(stderr,"Downsampling catalogue ...");
      downsampling(lc_cat);
      timer(1);
    }
    /*Halo abundance matching following Rodriguez-Torres et al.
      2016*/
    if(HAM){
      fprintf(stderr,"(Sub)Halo Abundance Matching ...");
      halo_abundance(lc_cat);
      timer(1);
    }
    free_gsl_nbar();

    if(MASK){
      fprintf(stderr,"Applying mask light-cone ...");
      apply_mask_ply(lc_cat);
      timer(1);}
    fprintf(stderr,"Writting output ...");
    if(WRITE_HDF5)
      write_hdf5_sel(lc_cat);
    else
      write_ascii_sel(lc_cat);
    timer(1);
  }  
}

void box_selection(box box_cat){

  /*Downsample the input catalogue to NBAR number density*/
  if(DOWNSAMPLE){
    fprintf(stderr,"Downsampling catalogue ...");
    downsampling_box(box_cat);
    timer(1);
  }

  /*Halo abundance matching following Rodriguez-Torres et al.
    2016*/
  if(HAM){
    fprintf(stderr,"(Sub)Halo Abundance Matching ...");
    halo_abundance_box(box_cat);
    timer(1);
  }
  fprintf(stderr,"Writting output ...");
  if(WRITE_HDF5)
    write_hdf5_sel_box(box_cat);
  else
    write_ascii_sel_box(box_cat);
  timer(1);
}
