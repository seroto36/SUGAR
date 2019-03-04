#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <time.h>  
#include "cosmo_tools.h"
#include "common.h"
#include "define.h"

void box_rotation(box box_cat,float alpha,float beta,float x0,float y0,float z0){

  long ii=0;
  double rr = alpha*D2R;
  double dr = beta*D2R;
  
#pragma omp parallel shared(dr,rr)
  {
#pragma omp for nowait schedule(static)
    for(ii=0; ii<box_cat.np; ii++){
      double x = box_cat.bx[ii].x-x0;
      double y = box_cat.bx[ii].y-y0;
      double z = box_cat.bx[ii].z-z0;
      double vx = box_cat.bx[ii].vx;
      double vy = box_cat.bx[ii].vy;
      double vz = box_cat.bx[ii].vz;
      
      box_cat.bx[ii].x = x*cos(dr)*cos(rr) - y*sin(rr) - z*sin(dr)*cos(rr);
      box_cat.bx[ii].y = x*cos(dr)*sin(rr) + y*cos(rr) - z*sin(dr)*sin(rr);
      box_cat.bx[ii].z = x*sin(dr) + z*cos(dr);

      box_cat.bx[ii].vx = vx*cos(dr)*cos(rr) - vy*sin(rr) - vz*sin(dr)*cos(rr);
      box_cat.bx[ii].vy = vx*cos(dr)*sin(rr) + vy*cos(rr) - vz*sin(dr)*sin(rr);
      box_cat.bx[ii].vz = vx*sin(dr) + vz*cos(dr);
    }
  }
}

void box_traslation(box box_cat,float x0,float y0,float z0){

   long ii=0;
#pragma omp parallel
   {
#pragma omp for nowait schedule(static)
     for(ii=0;ii<box_cat.np;ii++){
       box_cat.bx[ii].x = box_cat.bx[ii].x-x0;
       box_cat.bx[ii].y = box_cat.bx[ii].y-y0;
       box_cat.bx[ii].z = box_cat.bx[ii].z-z0;
     }
   }   
}

void remapping(box box_cat){

  long ii=0;
#pragma omp parallel for  
  for(ii=0; ii<box_cat.np; ii++){
    if(box_cat.bx[ii].x<0.0 && box_cat.bx[ii].y<0.0){
      box_cat.bx[ii].x = box_cat.bx[ii].x+sqrt(2)*LBOX*0.5;
      box_cat.bx[ii].y = box_cat.bx[ii].y+sqrt(2)*LBOX*0.5;
    }
    if(box_cat.bx[ii].x<0.0 && box_cat.bx[ii].y>=0.0){
      box_cat.bx[ii].x = box_cat.bx[ii].x+sqrt(2)*LBOX*0.5;
      box_cat.bx[ii].y = box_cat.bx[ii].y-sqrt(2)*LBOX*0.5;
    }
  }
}

void set_coordinates_box(int snap_num,box *box_cat){
  long ii=0;
  if(!strcmp(CONF_BOX, "NO_REMAP") ){
    if(RA_LOS != 0.0 || DEC_LOS != 0.0)
      for(ii = 0; ii < snap_num; ii++)
	box_rotation(box_cat[ii],RA_LOS,DEC_LOS,X_OBS,Y_OBS+LBOX*0.5,Z_OBS+LBOX*0.5);
    else
      for(ii = 0; ii < snap_num; ii++)
	box_traslation(box_cat[ii],X_OBS,Y_OBS+LBOX*0.5,Z_OBS+LBOX*0.5);
  }
  
  if(!strcmp(CONF_BOX, "SDSS_REMAP")){
    for(ii = 0; ii < snap_num; ii++){
      box_rotation(box_cat[ii],-45.0,0.0,LBOX*0.5,LBOX*0.5,LBOX*0.5);
      remapping(box_cat[ii]);
      box_traslation(box_cat[ii],X_OBS,Y_OBS,Z_OBS);
    }
    if(RA_LOS != 0.0 || DEC_LOS != 0.0)
      for(ii = 0; ii < snap_num; ii++)
	box_rotation(box_cat[ii],RA_LOS,DEC_LOS,0.0,0.0,0.0);
  }   
}

long count_and_sel_lc_halos(box box_cat,float min_dis,float max_dis){

  long ii=0,n_lc=0;
  float rc=0.0,ra=0.0,dec=0.0;   
  
#pragma omp parallel for private(rc,ra,dec)  reduction(+:n_lc)
  for(ii=0; ii<box_cat.np; ii++){
    float X2 = box_cat.bx[ii].x*box_cat.bx[ii].x;
    float Y2 = box_cat.bx[ii].y*box_cat.bx[ii].y;
    float Z2 = box_cat.bx[ii].z*box_cat.bx[ii].z;
    rc  = sqrt(X2+Y2+Z2);
    dec = asin(box_cat.bx[ii].z/rc)*R2D;
    ra=atan2(box_cat.bx[ii].y,box_cat.bx[ii].x)*R2D;

    if((fabs(ra-RA_LOS) <=RA_SIZE*0.5 || fabs(ra-RA_LOS-360.0) <=RA_SIZE*0.5 \
	|| fabs(ra-RA_LOS+360.0) <=RA_SIZE*0.5)
       && fabs(dec-DEC_LOS)<=DEC_SIZE*0.5){
      
      float LOS_x  = box_cat.bx[ii].x/rc;
      float LOS_y  = box_cat.bx[ii].y/rc;
      float LOS_z  = box_cat.bx[ii].z/rc;   
      float LOS_vx = box_cat.bx[ii].vx*LOS_x;
      float LOS_vy = box_cat.bx[ii].vy*LOS_y;
      float LOS_vz = box_cat.bx[ii].vz*LOS_z;
      float z      = dist2redshift(rc);
      float vr     = LOS_vx + LOS_vy + LOS_vz;
      float H      = 100.0*(sqrt((1+z)*(1+z)*(1+z)*OMEGA_M+OMEGA_L));
      float x_rds  = box_cat.bx[ii].x+(1.0+z)*vr*LOS_x/H;
      float y_rds  = box_cat.bx[ii].y+(1.0+z)*vr*LOS_y/H;
      float z_rds  = box_cat.bx[ii].z+(1.0+z)*vr*LOS_z/H;
      float rc0    = sqrt(x_rds*x_rds+y_rds*y_rds+z_rds*z_rds);
      
      if(rc0 < max_dis && rc0 >= min_dis){
	//if(rc < max_dis && rc >= min_dis){
	float zobs   = dist2redshift(rc0);
	if(ra<0.0)ra=ra+360.0;
	/*Here, the codes modifies the arrays of boxes, these data cannot
	  be used again to extract positions and velocities of boxes*/ 
	box_cat.bx[ii].aux = box_cat.bx[ii].vx;
	/*box_cat.bx[ii].M200  = box_cat.bx[ii].x;
	box_cat.bx[ii].Rv    = box_cat.bx[ii].y;
	box_cat.bx[ii].Vpeak = box_cat.bx[ii].z;*/
	box_cat.bx[ii].x   = ra;
	box_cat.bx[ii].y   = dec;
	box_cat.bx[ii].z   = zobs;
	box_cat.bx[ii].vx  = z;
	box_cat.bx[ii].sel = TRUE;
	n_lc++;
      }    
    }
  }
  fprintf(stderr,"TOTAL NUMBER OF OBJECTS %ld\n",n_lc);
  return n_lc;
}

void fill_lccat(long index,int snapnum, box box_cat,lightcone lc_cat){
  long ii=0;
  long jj=index;
  
  for(ii=0; ii<box_cat.np; ii++){
    if(box_cat.bx[ii].sel){
      lc_cat.lc[jj].ra      = box_cat.bx[ii].x;
      lc_cat.lc[jj].dec     = box_cat.bx[ii].y;
      lc_cat.lc[jj].rds     = box_cat.bx[ii].z;
      lc_cat.lc[jj].zreal   = box_cat.bx[ii].vx;
      /*lc_cat.lc[jj].vx      = box_cat.bx[ii].aux;
      lc_cat.lc[jj].vy      = box_cat.bx[ii].vy;
      lc_cat.lc[jj].vz      = box_cat.bx[ii].vz;*/
      lc_cat.lc[jj].Mv      = box_cat.bx[ii].Mv;
      lc_cat.lc[jj].Rv      = box_cat.bx[ii].Rv;
      lc_cat.lc[jj].M200    = box_cat.bx[ii].M200;
      lc_cat.lc[jj].Vpeak   = box_cat.bx[ii].Vpeak;
      lc_cat.lc[jj].Vmax    = box_cat.bx[ii].Vmax;
      lc_cat.lc[jj].str_id  = box_cat.bx[ii].str_id;
      lc_cat.lc[jj].id      = box_cat.bx[ii].id;
      lc_cat.lc[jj].lc_id   = jj;
      lc_cat.lc[jj].boxnum  = ii;
      lc_cat.lc[jj].snapnum = snapnum;
      lc_cat.lc[jj].sel     = FALSE;
      jj++;
    }
  }
}

void construct_lightcone(inter interv, box *box_cat, lightcone *lc_cat){

  fprintf(stderr,"Building light-cone... ");
  
  long N_OBJ[interv.ns_in],N_CUM[interv.ns_in];
  float h=HUBBLE_PAR;
  int ii=0;
  long index=0;
  set_coordinates_box(interv.ns_in,box_cat);

  for(ii = 0; ii < interv.ns_in; ii++){
    N_OBJ[ii]=count_and_sel_lc_halos(box_cat[ii],h*interv.lim_dis[ii+1],
				     h*interv.lim_dis[ii]);
    if(ii==0) N_CUM[ii] = N_OBJ[ii];
    else N_CUM[ii] = N_CUM[ii-1] + N_OBJ[ii];
  }
  
  allocate_mem_lightcone(N_CUM[interv.ns_in-1],lc_cat);

#pragma omp parallel for private(index)
  for(ii = 0; ii < interv.ns_in; ii++){
    if(ii==0) index=0;
    else index=N_CUM[ii-1];
    int snapnum = interv.snap_num[interv.first_sn+ii];    
    fill_lccat(index,snapnum,box_cat[ii],(*lc_cat));
  }

  timer(1);
}



