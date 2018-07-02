#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include "define.h"
#include "cosmo_tools.h"

double f_dis(double z, void * p){
  float h=HUBBLE_PAR;
  float Om=OMEGA_M;
  float Ol=OMEGA_L;
  return((299792.458/(100.0*h*sqrt((1.0+z)*(1.0+z)*(1.0+z)*Om+Ol))));
}

float com_distance(float z){
   
  if(z < 0.0){
    fprintf(stderr, "negative redshift");
    exit(1);
  }
  else{  
    if(z==0)return 0; 
    else{
      double error = 0.0, sum = 0.0;
      size_t neval;      
      gsl_function F;
      
      F.function = &f_dis;
      F.params = NULL;
      
      gsl_integration_qng(&F,0,z,0,1e-3,&sum,&error,&neval);
      
      return sum;
    }
  }
}

float dist2redshift(float r){
   
   float rc=r/HUBBLE_PAR;
   double z=0.0;

   gsl_spline_eval_e(spline_rc,rc,Acc_rc,&z);

   return z;
}

void r_to_z(){
  int i=0,steps=10000;
  double binZ=(maxZitp-minZitp)/(double)steps;
  double red[steps];
  double X[steps];

  for(i=0;i< steps;i++){
    red[i]=minZitp+i*binZ;
    double z=red[i];
    double r=com_distance(z);
    X[i]=r;
  }

  Acc_rc=gsl_interp_accel_alloc();
  spline_rc=gsl_spline_alloc(gsl_interp_cspline,steps);
  gsl_spline_init(spline_rc,X,red,steps);

}

void end_r_to_z(){
  gsl_interp_accel_free(Acc_rc);
  gsl_spline_free(spline_rc);
}

float Volume(float rmin,float rmax){
 
  float rs=RA_SIZE*D2R;
  float ds=DEC_SIZE*D2R;
  float ld=DEC_LOS*D2R;
  float radial=(rmax*rmax*rmax-rmin*rmin*rmin);
  float theta=sin(ld+ds*0.5)-sin(ld-ds*0.5);
  float V;
  if(READ_LC)
    V=AREA*D2R*D2R*radial/3.0;
  if(!READ_LC)
    V=rs*radial*theta/3.0;
  float h3=HUBBLE_PAR*HUBBLE_PAR*HUBBLE_PAR;
  return V*h3;
}

float number_density(double z){
  double nbar=0.0;
  gsl_spline_eval_e(spline_nbar,z,Acc_nbar,&nbar);
  return nbar;
}

float norm_number_density(double z){
  double nbar=0.0;
  gsl_spline_eval_e(spline_NORMnbar,z,Acc_NORMnbar,&nbar);
  return nbar;
}
