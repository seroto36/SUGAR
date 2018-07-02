#include <gsl/gsl_spline.h>
#include "define.h"

#define minZitp 0.0
#define maxZitp 20.0
#define P0 20000.0

/* gsl vectors */
gsl_interp_accel *Acc_NORMnbar;
gsl_interp_accel *Acc_nbar;
gsl_interp_accel *Acc_CMF;
gsl_interp_accel *Acc_rc;
gsl_spline *spline_NORMnbar;
gsl_spline *spline_nbar;
gsl_spline *spline_CMF;
gsl_spline *spline_rc;

float com_distance(float z);
double f_dis(double x, void * p);
float dist2redshift(float r);
void r_to_z();
void end_r_to_z();
float Volume(float rmin,float rmax);
float number_density(double z);
float norm_number_density(double z);
