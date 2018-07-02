#include "define.h"

char *CODE_OPT; /*BOX or LGC*/

//other variables
float MAX_NBAR=DEFAULT;
long NSEL=DEFAULT;
long NS_RAN=0;
float RC_MAX=DEFAULT;
float RC_MIN=DEFAULT;

//Config variables
int MASK=DEFAULT;
int RANDOM=DEFAULT;
int WRITE_HDF5=DEFAULT;
int READ_HDF5=DEFAULT;
int READ_LC=DEFAULT;
int PART_CAT=DEFAULT;
int HAM=DEFAULT;
int HOD=DEFAULT;
int MODIFIED_HAM=DEFAULT;
int DOWNSAMPLE=DEFAULT;
int COR_FUN=DEFAULT;
int FIND_PARM=DEFAULT;
int VAR_CUT=DEFAULT;
int NO_SEL=FALSE;

//param variables
int SHUFFLE=FALSE;
int ANG_DOWN=FALSE;
int CSMF=FALSE;
int NPART=DEFAULT;
int FRAN=1;
int ROTBOX=1;
float NBAR=DEFAULT;
float ZBOX=DEFAULT;
float AREA=DEFAULT;
float OMEGA_M=DEFAULT;
float OMEGA_L=DEFAULT;
float HUBBLE_PAR=DEFAULT;
float MINZ=DEFAULT;
float MAXZ=DEFAULT;
float RA_LOS=DEF_COOR;
float DEC_LOS=DEF_COOR;
float RA_SIZE=FALSE;
float DEC_SIZE=FALSE;
float X_OBS=DEF_COOR;
float Y_OBS=DEF_COOR;
float Z_OBS=DEF_COOR;
float LBOX=DEFAULT;
float MIN_COMP=FALSE;
float SCAT_HAM=DEFAULT;
float FSAT=DEFAULT;
char CONF_BOX[CHAR_SIZE]="UNDEFINED";
char F_RDS_SIM[CHAR_SIZE]="UNDEFINED";
char INFILE[CHAR_SIZE]="UNDEFINED";
char BASEFILE[CHAR_SIZE]="UNDEFINED";
char ENDFILE[CHAR_SIZE]="UNDEFINED";
char OUTFILE[CHAR_SIZE]="UNDEFINED";
char RANFILE[CHAR_SIZE]="UNDEFINED";
char MASKFILE[CHAR_SIZE]="UNDEFINED";
char NBARFILE[CHAR_SIZE]="UNDEFINED";
char SMFCAT[CHAR_SIZE]="UNDEFINED";

/*IMPORTANT MESSAGE: set N_VAR_LC before include new objects in PARAM_NAME
also include the global variables and modify the ini_var_conf_lc*/

char *fconf; //Config file
char *fparam; //Config file
int PAR_VAL_LC[N_VAR_LC]; 
char *PAR_NAME_LC[]={"MASK",
		     "RANDOM",
		     "READ_HDF5",
		     "WRITE_HDF5",
		     "READ_LC",
		     "PART_CAT",
		     "HAM",
		     "HOD",
		     "MODIFIED_HAM",
		     "DOWNSAMPLE",
		     "COR_FUN",
		     "VAR_CUT",
		     "FIND_PARM"};

/*IMPORTANT MESSAGE: set N_VAR_BX before include new objects in PARAM_NAME
also include the global variables and modify the ini_var_conf_box*/
int PAR_VAL_BX[N_VAR_BX];
char *PAR_NAME_BX[]={"COR_FUN",
		     "READ_HDF5",
		     "WRITE_HDF5",
		     "HAM",
		     "HOD",
		     "MODIFIED_HAM",
		     "DOWNSAMPLE",
		     "VAR_CUT",
		     "FIND_PARM"};

/*Varibles fix in param file, this variables will depends on the selecctions 
  in the configure file*/

char *VAR_VALUE[N_VAR]={0};
char *VAR_NAME[]={"OMEGA_M",
		  "OMEGA_L",
		  "AREA",
		  "HUBBLE_PAR",
		  "MINZ",
		  "MAXZ",
		  "RA_LOS",
		  "DEC_LOS",
		  "RA_SIZE",
		  "DEC_SIZE",
		  "FRAN",
		  "SHUFFLE",
		  "INFILE",
		  "ROTBOX",
		  "BASEFILE",
		  "ENDFILE",
		  "OUTFILE",
		  "RANFILE",
		  "SMFCAT",
		  "CSMF",
		  "MASKFILE",
		  "F_RDS_SIM",
		  "NPART",
		  "CONF_BOX",
		  "SCAT_HAM",
		  "X_OBS",
		  "Y_OBS",
		  "Z_OBS",
		  "LBOX",
		  "ZBOX",
		  "NBAR",
		  "MIN_COMP",
		  "ANG_DOWN",
		  "NBARFILE",
		  "FSAT"};
