const char *Xh5,*Yh5,*Zh5,*Vxh5,*Vyh5,*Vzh5;

//Structures
typedef struct {
  float ra,dec,rds,zreal,wfkp;
  float Mv,M200,Rv,Vpeak,Vmax,Mstar;
  long boxnum,str_id,id,lc_id;
  int snapnum,sel;
  float vx,vy,vz;
} lc_data;

typedef struct {
  float ra,dec,rds,wfkp,zreal;
  int sel;
} ran_data;

typedef struct{
  int sel;
  long id,str_id,box_id;
  float zrds,Mstar;
  float x,y,z,vx,vy,vz;
  float aux;
  float M200,Mv,Rv,Vmax,Vpeak;
} box_data;

typedef struct{
  long np;
  lc_data *lc;
} lightcone;

typedef struct{
  long np;
  ran_data *lc;
} ran_lc;

typedef struct{
  long np,ntot;
  box_data *bx;
}box;

typedef struct{
  long pst;
  float var;
} proxy;

typedef struct{
  long nsel;
  long np;
  float vol;
  proxy *prx;
}hl_proxy;

typedef struct{
  float ra,dec,rds,zreal,Mstar,wfkp;
  float Rv,Vmax,Mv,M200,Vpeak;
  long str_id,id,boxnum,lc_id;
  int snapnum; 
} hdf5_data;

typedef struct{
  float ra,dec,rds,zreal;
  float Rv,Vmax,Mv,M200,Vpeak;
  long str_id,id,boxnum,lc_id;
  int snapnum; 
} hdf5_data_r;

typedef struct{
  float x,y,z,vx,vy,vz,M200;
  float Rv,Vmax,Vpeak,Mv;
  long str_id,id;
} hdf5_box;

typedef struct{
  float x,y,z,vx,vy,vz,M200;
  float zrds,Mstar;
  float Rv,Vmax,Vpeak,Mv;
  long str_id,id;
} hdf5_box_w;

typedef struct{
  int ns_in; //number of snapshots within the rds. range
  int ns_tot; //total number of snapshots in the redshift file
  int n_lim; //total number of redshift divisions in the lightcone 
  int first_sn, last_sn; // positions first and lst snapshots used
  int *snap_num; //number of each snapshot
  float *snap_rds; //redshift of each snapshots
  float *lim_rds; //limits in redshift for each shell
  float *lim_dis; //limits in distance for each shell
} inter;

typedef struct{
  long np;
  float z_min,z_max;
  float r_min,r_max;
  float vol,nbar;
}sheels;

typedef struct{
  int   *tot_obj;
  int   *i_obj;
  float *f_obj;
  float *prob;
} mass_sel;

//initialisation functions
void timer(int ii);
long n_halos_ascii(char *file);
void ini_conf_values();
void ini_var_conf_lc();
void ini_var_conf_box();
void ini_var_param();
void check_conf_lc();
void check_conf_bx();
void check_params();
void init_lightcone_struct(lightcone *light_cat);
void init_lightcone_rand(ran_lc *light_cat);
void init_box_struct(box *box_cat);
void init_boxes_struct(inter interv,box **box_cat);
void init_inter_struct(inter *interv);
void interval_distances(inter *interv);
void init_sheels_ham(hl_proxy **ham_px);
void init_mass_sel(mass_sel **sm_sel);

//Allocation functions
void load_lightcone(lightcone *light_cat);
void load_box(box *box_cat);
void load_boxes(inter interv,box **box_cat);
void free_lightcone(lightcone cat);
void free_random(ran_lc cat);
void free_box(box cat);
void free_boxes(int num,box *cat);
void free_inter(inter interv);
void allocate_mem_random(long size,ran_lc *light_cat);
void allocate_mem_lightcone(long size,lightcone *light_cat);
void allocate_mem_box(long size,box *box_cat);
void allocate_sheels_ham(hl_proxy **ham_px,lightcone lc_cat);
void free_gsl_nbar();
void free_gsl_NORMnbar();
void free_gsl_CMF();

//reading fuctions
void read_conf_file();
void read_param_file();
void read_argc(int argc, char **argv);
void read_hdf5_lightcone(lightcone *data);
void read_lightcone_ascii(lightcone lc);
void read_box_ascii(char *fname,box box_cat);
void read_redshift(inter interv);
void read_hdf5_box(char *fname,box *box_cat);
void read_hdf5_boxes(char *fname,box *box_cat);
void read_radial_nbar(int NORM,sheels parts[]);
void read_angular_random(ran_lc lc_ran);
void ReadandClass_SM(mass_sel *sm_sel,hl_proxy *ham_px);
float fast_atof (char *p);
long n_halos_ascii(char *file);

//writting functions
void write_hdf5(lightcone lc_cat);
void write_ascii(lightcone dat);
void write_hdf5_sel(lightcone lc_cat);
void write_ascii_sel(lightcone dat);
void write_output_info(inter interv,int tot);
void write_ascii_sel_box(box dat);

//selection functions
void lightcone_selection(lightcone lc_cat);
void construct_random(ran_lc *lc_ran,lightcone lc_cat);
void apply_mask_ply(lightcone lc_cat);
void fill_sheels(hl_proxy *ham_px,lightcone lc_cat);
void fill_proxy_box(proxy *ham_px,box bx_cat);
void box_selection(box box_cat);

//error message functions
void error_line_numvar(int nline, char *file);
void error_open_file(char *fname);
void error_line_cfile(int nline);
void error_input_param(int nline,char *name,char *fname);
void error_missing_param(char *fname, char *param);
void error_repeat_param(int nline, char *fconf, char *name);
void error_malloc_allocation(char *var);

//set characteristics of light-cone
void set_redshift_intervals(inter *intervals);
void construct_lightcone(inter interv, box *box_cat, lightcone *lc_cat);
long count_lightcone_halos(inter interv, box box_cat,float min_dis,float max_dis);
void CSMF_calculator();
