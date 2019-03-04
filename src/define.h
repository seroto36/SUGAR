#define FALSE 0
#define TRUE !FALSE
#define DEFAULT -1
#define DEF_COOR 1e7
#define CHAR_SIZE 1024
#define SMCUT 11.0 //stellar mass cut SMF
#define MIN_MASS 9.5 //Minimum stellar mass
#define MAX_MASS 13.0 //Maximum stellar mass 
#define MASS_NBIN 50 //Number of bins for mass distribution
#define N_VAR_LC 13 //Number of options in lgc configure file
#define N_VAR_BX 9 //Number of options in box configure file
#define N_VAR 35 //Number of posible variables in the parameters file
#define DATASET "DM-lc"
#define DATABOX "DM-HALOS"
#define N_COL_LC_ASCII 11
#define N_COL_BOX_ASCII 13
//#define N_COL_BOX_ASCII 6
#define R2D 57.2957795131
#define D2R 0.017453292519943295

extern char *fconf; //Config file
extern char *fparam; //Param file
extern char *CODE_OPT; //option box or lightcone
extern int PAR_VAL_LC[N_VAR_LC];
extern char *PAR_NAME_LC[];
extern int PAR_VAL_BX[N_VAR_BX];
extern char *PAR_NAME_BX[];
extern char *VAR_VALUE[];
extern char *VAR_NAME[];
extern float MAX_NBAR;
extern long NSEL;
extern long NS_RAN;
extern float RC_MAX;
extern float RC_MIN;

//Config variables
extern int MASK;
extern int RANDOM;
extern int READ_HDF5;
extern int WRITE_HDF5;
extern int READ_LC;
extern int PART_CAT;
extern int HAM;
extern int HOD;
extern int MODIFIED_HAM;
extern int DOWNSAMPLE;
extern int COR_FUN;
extern int FIND_PARM;
extern int VAR_CUT;
extern int NO_SEL;

//Param variables
extern int SHUFFLE; //Make randoms using shuffle method
extern int ANG_DOWN; //Apply angular downsampling using .ply mask
extern int NPART;  //Number of divisions in redshift in case of n(z)
extern int FRAN; //Fraction between N randoms and data FRAN=NDAT/NRAN
extern int CSMF; //use cumulative stellar masss func. to include masses
extern int ROTBOX; //initial rotation of boxes
extern float ZBOX; //Redshift of box in BOX mode
extern float NBAR; //Number density of the final catalogue in BOX mode
extern float AREA; //Angular area of the input lightcones (deg^2)
extern float OMEGA_M; //From simulation
extern float OMEGA_L;  //From simulation
extern float HUBBLE_PAR;  //From simulation
extern float MINZ; //Minimum redshift of lightcone
extern float MAXZ; //Maximum redshift of lightcone
extern float RA_LOS; //RA size of lc
extern float DEC_LOS; //DEC size of lc
extern float RA_SIZE; //Line of sight of the lightcone's center
extern float DEC_SIZE; //Line of sight of the lightcone's center
extern float LBOX; //Line of sight of the lightcone's center
extern float MIN_COMP; //minimum completeness for selecting a polygon
extern float SCAT_HAM; //scatter halo abundance matching
extern float FSAT; //Fraction of subhalos in the final catalogue

/*Position of the observer will depends on the CONF_REMAP option. In the case
of SDSS_REMAP the center of the coordinates will be in the center of one of the
larger faces. Remember, that this code contruct the light-cone by putting the 
LOS along the x-axe. In the case of NO_REMAP the coordinate origin is in the 
middle of the Y-Z plane form by one of the faces of the box.*/
extern float X_OBS;
extern float Y_OBS;
extern float Z_OBS;

extern float DEC_SIZE; //Line of sight of the lightcone's center
extern float DEC_SIZE; //Line of sight of the lightcone's center
extern char CONF_BOX[CHAR_SIZE]; //SDSS_REMAP NO_REMAP
extern char INFILE[CHAR_SIZE]; //input file for single file reading
extern char BASEFILE[CHAR_SIZE]; //base file for snapshots reading
extern char ENDFILE[CHAR_SIZE]; //end extension file for snapshots reading
extern char NBARFILE[CHAR_SIZE]; //n(z) file
extern char OUTFILE[CHAR_SIZE]; //output file name
extern char RANFILE[CHAR_SIZE]; //random file name
extern char MASKFILE[CHAR_SIZE]; //mask file name
extern char F_RDS_SIM[CHAR_SIZE]; //redshift list of simulation
extern char SMFCAT[CHAR_SIZE]; //stellar mass catalogue
