/*
 * transit.c - Models the modulation produced by a Planet crossing in
 *             front of the star. Main component of the Transit program.
 *
 * Copyright (C) 2003 Patricio Rojo (pato@astro.cornell.edu)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */
//\omitfh

/* TD: Command line parameters parsing */
/* TD: calloc checks */
/* TD: atmosphere info  retrieval */

#include <transit.h>
#include <math.h>
/*
#include "transitstd.c"
#include "readlineinfo.c"
#include "../util/voigt.c"
#include "../util/sampling.c"
*/

/* Version history:
   0.3: First light, it only calculates Kappa values of it and of
        lorentz widening seem rather high. 102703. Patricio Rojo
   0.5: Lorentz and doppler width work all right now. 102803. PR
 */
static int version=1;
static int revision=-1;
extern int verblevel;              /* verbose level, greater than 10 
				      is only for debuging */

static inline PREC_ZREC interpolateTR(PREC_ZREC *x,
				      PREC_ZREC *y,
				      int n,
				      PREC_ATM val);

static inline int makegauTR(int *nsmr,
			    double **smr,
			    double width,
			    double dl,
			    float mgau);

static inline double integrateTR(PREC_RES *spectra,
				 int osw);

enum param {
  CLA_DUMMY=128;
  CLA_ATMOSPHERE,
  CLA_LINEDB
};


struct optdocs const var_docs[]={
  {NULL,HELPTITLE,0,
   NULL,"GENERAL ARGUMENTS"},
  {"version",no_argument,'v',
   NULL,"Prints version number and exit"},
  {"help",no_argument,'h',
   NULL,"Prints list of possible parameters"},
  {"defaults",no_argument,'d',
   NULL,"Prints default values of the different variable"},

  {NULL,HELPTITLE,0,
   NULL,"INPUT/OUTPUT"},
  {"output",required_argument,'o',
   "outfile","Change output file name, a dash (-) "
   "directs to standard output"},
  {"atmosphere",required_argument,CLA_ATMOSPHERE,
   "atmfile","File containing atmospheric info (Radius, "
   "pressure, temperature). A dash (-) indicates alternative "
   "input."},
  {"linedb",required_argument,CLA_LINEDB,
   "linedb","File containing line information (TWII format, "
   "as given by 'lineread'"},

  {NULL,HELPTITLE,0,
   NULL,"RADIUS OPTIONS (all in planetary radii units)"},
  {"radius",no_argument,'r',
   NULL,"Interactively input radius parameters"},
}

/* \fcnfh
   Output a syntax help message.

   @returns Never It either exit with the syntax help or with a error
                  message if number of parameters are wrong
*/
void synhelp_transit(const char *unknown,const char par,
		     const struct transithint *th)
{
  //Following is the help text  
  char message[]="Syntax:\n\ttransit [..options..]\n\n"
    " Where, sorted by effect, the available options are:\n"
    "\n GENERAL OPTIONS:\n"
    "  -h              Show this help\n"
    "  -v [+..][-..]   Increase or decrease verbose level by one per\n"
    "                  each + or -. 0 is the quietest, %i is the\n"
    "                  noisiest (%i)\n"
    "  -V              Show program's version and exit.\n"
    "\n FILE OPTIONS\n"
    "  -f a<atm_file>  Name of the atmospheric data file (\"%s\")\n"
    "  -f l<line_file> Name of the line info file produced by\n"
    "                  lineread (\"%s\")\n"
    "  -f o<out_file>  Name of the output file (\"%s\")\n"
    "\n RADIUS OPTIONS (all in planetary radius units)\n"
    "  -r i<low_rad>   Lower radius. 0 if you want to use atmospheric\n"
    "                  data minimum (%g).\n"
    "  -r f<high_rad>  Upper radius. 0 if you want to use atmospheric\n"
    "                  data maximum (%g).\n"
    "  -r d<delta_rad> Radius spacing. 0 if you want to use atmospheric\n"
    "                  data spacing (%g).\n"
    "\n WAVELENGTH OPTIONS (all in nanometers)\n"
    "  -w i<low_wav>   Lower wavelength. 0 if you want to use line\n"
    "                  data minimum (%g).\n"
    "  -w f<high_wav>  Upper wavelength. 0 if you want to use line\n"
    "                  data maximum (%g).\n"
    "  -w d<delta_wav> Wavelength spacing. 0 if you want to use line\n"
    "                  data spacing (%g).\n"
    "  -w o<delta_wav> Wavelength oversampling (%i).\n"
    "  -m <margin>     Not trustable range in microns at boundary\n"
    "                  of line databases. Also transitions this\n"
    "                  much away from the requested range will be\n"
    "                  considered (%g)\n"
    "\n WAVENUMBER OPTIONS (all in cm-1)\n"
    "  -n i<low_wn>    Lower wavenumber. 0 if you want to use\n"
    "                  equivalent of the wavelength maximum (%g).\n"
    "  -n f<high_wn>   Upper wavenumber. 0 if you want to use\n"
    "                  equivalent of the wavelength minimum (%g).\n"
    "  -n d<delta_wn>  Wavenumber spacing. 0 if you want to have\n"
    "                  the same number of points as in the\n"
    "                  wavelength sampling in an equispaced\n"
    "                  grid. (%g)\n"
    "  -n o<delta_wn>  Wavenumber oversampling. 0 if you want\n"
    "                  the same value as for the wavelengths (%i).\n"
    "\n OPACITY CALCULATION OPTIONS:\n"
    "  -s <subbinning> Number of fine-bins to calculate the Voigt\n"
    "                  function (%i)\n"
    "  -a <timesalpha> Number of the max-alpha (the greater of Voigt\n"
    "                  or Doppler widths) that need to be contained\n"
    "                  in a calculated Voigt profile (%g)\n"
    "  -d <maxratio>   Ratio of maximum allowed doppler change\n"
    "                  before recalculating profile (%g)\n"
    "\n OBSERVATIONAL OPTIONS:\n"
    "  -t <tel_res>    Telescope resolution in nm. (%g)\n"
    "\n";

  //Check if I'm here because of a bad command. If so complain and
  //suggest help
  if(unknown){
    fprintf(stderr,
	    "Option '-%c %s' not available. Try 'transit -h' for help\n"
	    ,par,unknown);
  }
  //Else output syntax help
  else{
    fprintf(stderr,message, 
	    th->verbnoise, verblevel, 
	    th->f_atm,  th->f_line, th->f_out, 
	    th->rads.i, th->rads.f, th->rads.d,
	    th->wavs.i, th->wavs.f, th->wavs.d, th->wavs.o, th->m,
	    th->wns.i,  th->wns.f,  th->wns.d,  th->wns.o,
	    th->voigtfine, th->timesalpha, th->maxratio_doppler,
	    th->t
	    );

  }

  //stop execution of program no matter what.
  exit(EXIT_FAILURE);
}



//\fcnfh
int main (int argc,		/* Number of variables */
	  char **argv)		/* Variables*/
{

  //Initialization of data structure's pointers. Note that s\_lt is not
  //assigned because the array that is going to point has not been
  //initialized yet.
  struct transit transit;
  struct transithint trh;
  memset(&transit, 0, sizeof(struct transit));
  memset(&trh, 0, sizeof(struct transithint));
  transit.ds.th=&trh;

  //Initialization of wavelength sampling in nanometers.
  //Fields from 'trh.wavs'.
  //'i' and 'f' fields are for initial and final wavelength, if 0 then
  //the most extreme common minimum or maximum of the line databases is
  //chosen. 
  //'d' field is the wavelength spacing.
  //'n' is the number of sampled elements: it should always be set to 0
  //whenever 'd' is changed, so that 'makewavsample()' can update
  //correctly.
  //'o' is the oversampling.
  //'v' is the value's array.
  prop_samp *samp=&trh.wavs;
  samp->i=0;
  samp->f=0;
  samp->d=0.188;
  samp->n=0;
  samp->o=100;
  samp->v=NULL;
  trh.na|=TRH_WI|TRH_WF|TRH_WD|TRH_WO;

  //Initialization of radius sampling in planetary radius.
  //Fields from 'trh.rads'.
  //'.i', '.f', '.d', '.n', '.v' are equivalent to the wavelegth
  //sampling above but for the radius. If zeroed then '.i', '.f' or '.d'
  //are taken from the atmosphere datafile.
  //'o' equal to 0 is necessary to dissable oversampling in radius.
  samp=&trh.rads;
  samp->i=0;
  samp->f=0;
  samp->d=0.1;
  samp->n=0;
  samp->v=NULL;
  samp->o=0;
  trh.na|=TRH_RI|TRH_RF|TRH_RD;

  //Initialization of wavenumber sampling in cm-1. Note that the
  //extincitons are calculated in wavenumber and then optionally
  //interpolated to wavelength.
  //Fields from 'trh.wns'.
  //'.i', '.f', '.d', '.n', '.v' are equivalent to the wavelegth
  //sampling above but for the radius. If zeroed then '.i', '.f' are
  //taken from the wavelength transformation. '.d' is calculated so that
  //the same number of points as requested for wavelength are outputted.
  //'o' equal to 0 is necessary to dissable oversampling in radius.
  samp=&trh.wns;
  samp->i=0;
  samp->f=0;
  samp->i=5716.5;
  samp->f=5718;
  samp->d=0;
  samp->n=0;
  samp->v=NULL;
  samp->o=0;
  trh.na|=TRH_WNI|TRH_WNF|TRH_WND|TRH_WNO;

  //Initialization of general variables.
  //'rc' and 'rn' are the general auxiliary variables.
  //'file\_out' is the name of the output file, a '-' indicates standard
  //output.
  //'trh.verbnoise' is the noisiest possible verbose level
  //controlled by 'verbose' when in a non-debugging compilation, a value
  //of 10 or more is only used for debugging.
  //'defile\_out' is the default output filename
  int rn;
  char rc;
  trh.verbnoise=4;
  char defile_out[]="-";
  trh.f_out=(char *)calloc(strlen(defile_out)+1,sizeof(char));
  strcpy(trh.f_out,defile_out);
  trh.na|=TRH_FO;
#ifdef NODEBUG_TRANSIT
  verblevel=2;
#else
  verblevel=20;
#endif /* NODEBUG_TRANSIT */

  //Initialization of line database variables. 
  //'file\_line' is the name of the info file that is finally used
  //'trh.m' is the amount of microns that cannot be trusted at the
  //boundaries of databases range, and also how much extra out of
  //requested range it has to look for transitions.
  //'defile\_line' is the default name of the line info file.
  //(i.e. modifiable by the user)
  trh.m=0.001;
  trh.na|=TRH_WM;
  char defile_line[]="./res/lineread.inf";
  trh.f_line=(char *)calloc(strlen(defile_line)+1,sizeof(char));
  strcpy(trh.f_line,defile_line);
  trh.na|=TRH_FL;


  //Initialization of atmospheric parameters.
  //'.f\_atm' is the name of the file with the atmospheric parameters,
  //a '-' indicates that a one-point solution is desired.
  char defile_atm[]="-";
  trh.f_atm=(char *)calloc(strlen(defile_atm)+1,sizeof(char));
  strcpy(trh.f_atm,defile_atm);
  trh.na|=TRH_FA;
  trh.fl|=TRU_ATMHARDC1P|TRU_SAMPLIN;

  //Initialization of extinction parameters.
  //'.voigtfine' is the fine binning of voigt
  trh.voigtfine=5;
  trh.timesalpha=50;
  trh.maxratio_doppler=0.001;
  trh.na|=TRH_VF|TRH_TA|TRH_DR;
  trh.fl|=TRU_EXTINPERISO;

  //Command line parameters' processing
  opterr=0;
  while(1){
    rn=getopt(argc,argv,"f:Vhv:m:r:w:n:a:s:d:");
    if (rn==-1)
      break;

    transitDEBUG(20,verblevel,
		 "Processing option '%c', argum: %s\n"
		 ,rn,optarg);

    switch(rn){
    case 'w':			//Change wavelength sampling
      samp=&trh.wavs;
    case 'n':			//Change wavenumber sampling
      if(rn=='n') samp=&trh.wns;
    case 'r':			//Change radius sampling
      if(rn=='r') samp=&trh.rads;
      switch(*optarg++){
      case 'i':			//initial value
	samp->i=atof(optarg);
	break;
      case 'f':			//final value
	samp->f=atof(optarg);
	break;
      case 'd':			//spacing
	samp->d=atof(optarg);
	samp->n=0;
	samp->v=NULL;
	break;
      case 'o':			//oversampling (only for wavelength)
	if(rn=='w'){
	  samp->o=atof(optarg);
	  break;
	}
      default:
	synhelp_transit(optarg-1,rn,&trh);
	break;
      }
      break;

    case 'f':			//Change filenames.
      rc=*optarg;
      while(*++optarg==' ');	//Get rid of extra blanks.
      switch(rc){
      case 'a':			//Atmosphere file
	trh.f_atm=(char *)realloc(trh.f_atm,strlen(optarg)+1);
	strcpy(trh.f_atm,optarg);
	break;
      case 'l':			//Lineinfo file
	trh.f_line=(char *)realloc(trh.f_line,strlen(optarg)+1);
	strcpy(trh.f_line,optarg);
	break;
      case 'o':			//Output file
	trh.f_out=(char *)realloc(trh.f_out,strlen(optarg)+1);
	strcpy(trh.f_out,optarg);
	break;
      default:
	*--optarg=rc;
	synhelp_transit(optarg,rn,&trh);
	break;
      }
      break;

    case 't':			//Telescope resolution
      trh.t=atof(optarg);
      break;

    case 'm':			//Change 'trh.m' margin
      trh.m=atof(optarg);
      break;

    case 'd':			//Change Doppler's maximum accepted
				//ratio before recalculating
      trh.maxratio_doppler=atof(optarg);
      break;

    case 's':			//Change voigt fine-binning
      trh.voigtfine=atoi(optarg);
      break;

    case 'a':			//Change times of alphas in profile
      trh.timesalpha=atof(optarg);
      break;

    case 'v':			//Increase/decrease verbose level
      optarg--;
      while(*++optarg){
	if(*optarg=='+')
	  verblevel++;
	else if(*optarg=='-')
	  verblevel--;
      }
      if(verblevel<0)
	verblevel=0;
      break;

    case 'V':			//Print version number and exit
      printf("This is 'transit' version %i.%i\n\n",version,revision);
      exit(EXIT_SUCCESS);
      break;

    case '?':
      rn=optopt;
    default:			//Ask for syntax help
      synhelp_transit(optarg,rn,&trh);
      break;
    case 'h':
      synhelp_transit(NULL,0,&trh);
      break;
    }
  }

  //No program warnings if verblevel is 0 or 1
  if(verblevel<2)
    transit_nowarn=1;

  //Presentation
  transitprint(1,verblevel,
	       "                TRANSIT v%i.%i\n"
	       "-----------------------------------------------\n"
	       ,version,revision);

  //Read line info
  if((rn=readlineinfo(&transit))!=0)
    transiterror(TERR_SERIOUS,
		 "readlineinfo() returned error code %i\n"
		 ,rn);

  transitDEBUG(20,verblevel,
	       "Final limits: %g %g\n"
	       ,transit.wavs.i,transit.wavs.f);

  //Read Atmosphere information
  if((rn=getatm(&transit))!=0)  
    transiterror(TERR_SERIOUS,
		 "getatm() returned error code %i\n"
		 ,rn);


  //Make wavelength binning
  if((rn=makewavsample(&transit))<0)
    transiterror(TERR_SERIOUS,
		 "makewavsample() returned error code %i\n"
		 ,rn);
  if(rn>0)
    transitprint(1,verblevel,
		 "makewavsample() modified some parameters according\n"
		 "to returned flag: 0x%x\n"
		 ,rn);

  //Make wavenumber binning
  if((rn=makewnsample(&transit))<0)
    transiterror(TERR_SERIOUS,
		 "makewnsample() returned error code %i\n"
		 ,rn);
  if(rn>0)
    transitprint(1,verblevel,
		 "makewnsample() modified some parameters according\n"
		 "to returned flag: 0x%x\n"
		 ,rn);

  //Make radius binning and interpolate data to new value
  if((rn=makeradsample(&transit))<0)
    transiterror(TERR_SERIOUS,
		 "makeradsample() returned error code %i\n"
		 ,rn);
  if(rn>0)
    transitprint(1,verblevel,
		 "makeradsample() modified some parameters according\n"
		 "to returned flag: 0x%x\n"
		 ,rn);


  //Calculates extinction coefficient
  if((rn=extwn(&transit))!=0)
    transiterror(TERR_SERIOUS,
		 "extwn() returned error code %i\n"
		 ,rn);


  printf("#wavenumber[cm-1]\twavelength[nm]\textinction[cm-1]\tcross-section[cm2]\n");
  for(rn=0;rn<transit.wns.n;rn++)
    printf("%10.4f%10.4f%15.5g%15.5g\n"
	   ,transit.wns.v[rn],WNU_O_WLU/transit.wns.v[rn],
	   transit.ds.ex->k[0][0][rn],
	   transit.ds.ex->k[0][0][rn]*transit.isof[0].m/transit.isov[0].d[0]);




  return EXIT_SUCCESS;
}

/* TD: atmosphere info  retrieval */
/*\fcnfh
   getatm: Get atmospheric parameters from a given file. It has to get
   all the extra isotopes data that was not in line database.

   @returns 0 On success
            -1 If no atmospheric info file was specified, and no
               defaults are allowed
	    -2 Default handling mode inexistent
	    -3 Something really bad happened!
*/
int getatm(struct transit *tr) /* Containing filename of atmosphere
				       info file */
{
  //'nmb' auxiliar variable for a number.
  int nmb,i;
  //'at' and 'hints' are the atmospheric and hint's structure,
  //respectively
  struct transithint *th=tr->ds.th;
  static struct atm_data st_at;
  memset(&st_at,0,sizeof(struct atm_data));
  tr->ds.at=&st_at;
  /*  struct atm_data *at=tr->ds.at;*/

  //Hard coded values.
  //'hc\_n\_e' number of extra isotopes (eiso).
  //'hc\_name' names.
  //'hc\_t' temperature.
  //'hc\_mass' masses.
  //'hc\_cs' cross section.
  //'hc\_dens' density.
#define HC_N_E 1
#define MAX_NAME 30
  int hc_n_e=0;
  char hc_name[HC_N_E][MAX_NAME]={""};
  PREC_ZREC hc_t=1350;
  PREC_ZREC hc_mass[HC_N_E]={0};
  PREC_CS hc_cs[HC_N_E]={0};
  PREC_ATM hc_dens[hc_n_e + tr->n_i];
  PREC_ATM hc_pres=1.0e3;
  //write defult density for the standard isotopes
    /* TD: Density units */
  hc_dens[0]=6.5e-4*hc_pres*(tr->isof[0].m)/(KB)/(hc_t);
  for(i=1;i<(tr->n_i);i++)
    hc_dens[i]=6e-10*hc_pres*(tr->isof[i].m)/(KB)/(hc_t);
  //and now for the extra isotopes if any
  if(hc_n_e){
    hc_dens[tr->n_i+1]=0;
  }

#undef HC_N_E
#undef MAX_NAME

  //Pass atmospheric flags into info struct
  transitacceptflag(tr->fl,th->fl,TRU_ATMBITS);

  //If filename is not given or is "-", then do not use atmosphere file
  //and use defaults instead. 
  if(th->f_atm==NULL||strcmp(th->f_atm,"-")==0){
    transitaccepthint(tr->f_atm,th->f_atm,th->na,TRH_FA);
    st_at.rads.n=1;
    st_at.rads.v=(PREC_ATM *)calloc(st_at.rads.n,sizeof(PREC_ATM));
    st_at.atm.t=(PREC_ATM *)calloc(st_at.rads.n,sizeof(PREC_ATM));
    st_at.rads.v[0]=1.0;
    //See which way does the user wants the default to be handled
    switch(tr->fl&TRU_ATM1PBITS){
    case TRU_ATMHARDC1P:	/* wants to use hard-coded values */
      st_at.atm.t[0]=hc_t;
      st_at.n_niso=hc_n_e;
      nmb=tr->n_e=tr->n_i+st_at.n_niso;
      tr->isof=(prop_isof *)realloc(tr->isof,nmb*sizeof(prop_isof));
      tr->isov=(prop_isov *)realloc(tr->isov,nmb*sizeof(prop_isov));
  transitDEBUG(20,verblevel,
	       "First isotope: %s\n"
	       ,tr->isof[0].n);
      st_at.isov=(prop_isov *)calloc(nmb,sizeof(prop_isov));
      for(i=0;i<nmb;i++){
	st_at.isov[i].d=(PREC_ATM *)calloc(1,sizeof(PREC_ATM));
	st_at.isov[i].d[0]=hc_dens[i];
      }
      nmb=tr->n_i;
      for(i=0;i<st_at.n_niso;i++){
	tr->isof[nmb+i].n=(char *)calloc(strlen(hc_name[i]),
					 sizeof(char));
	strcpy(tr->isof[nmb+i].n,hc_name[i]);
	tr->isof[nmb+i].m=hc_mass[i];
	st_at.isov[nmb+i].c=(PREC_CS *)calloc(1,sizeof(PREC_CS));
	st_at.isov[nmb+i].c[0]=hc_cs[i];
      }
      
      break;
    case TRU_ATMNODEF:		/* Wants an error message */
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "getatm():: No atmospheric file specified");
      return -1;
      break;
    case TRU_ATMGIVEN1P:	/* Wants to get it from somewhere */
    case TRU_ATMASK1P:		/* Wants to give it to standard input */
      transiterror(TERR_CRITICAL,
		   "getatm():: Defaults handling mode (0x%x) not\n"
		   "yet implemented\n"
		   ,tr->fl&TRU_ATM1PBITS);
      break;
    default:			/* Oops */
      transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
		   "getatm():: Unexistent default handling mode (0x%x)\n"
		   "requested\n"
		   ,tr->fl&TRU_ATM1PBITS);
      return -2;
      break;
    }

    //Let user know that defaults are being used
    transitprint(1,verblevel,
		 "You are using one point atmospheric conditions:\n"
		 " Temperature: %g K\n"
		 ,st_at.atm.t[0]);
    //Densities for all isotopes
    for(i=0;i<tr->n_i;i++)
      transitprint(1,verblevel,
		   " %s: density %8g g/cm3\n"
		   ,tr->isof[i].n,st_at.isov[i].d[0]);
    //Density and cross section for extra isotopes
    for(;i<tr->n_e;i++)
      transitprint(1,verblevel,
		   " %s: density %.8g g/cm3\tcross section %8g\n"
		   ,tr->isof[i].n,st_at.isov[i].d[0],
		   st_at.isov[i].c[0]);

  }
  //Otherwise if info is to be taken from file.
  else{
  /* TD: atmosphere info retrieval */
    //If a filename was given check its existence and that it can be
    //opened.
    verbfileopen(th->f_atm,&tr->fp_atm,"Atmospheric info ");

    transiterror(TERR_CRITICAL,
		 "Read atmospheric data file not yet implemented!\n");
  }
  transitASSERT(tr->n_e>tr->n_i,
		"Uyuyuyy!, number of isotopes after extension (%i)\n"
		"is smaller than number of isotopes before that (%i)\n"
		,tr->n_e,tr->n_i);

  //Return succes and set progress indicator
  tr->pi|=TRPI_GETATM;
  return 0;
}


/* \fcnfh
   Creates the sample points from hinted values

   @returns TRH\_S?<<bitshift for modified input
             0 if nothing was changed but there is a sampled array
	    -1 if hinted initial is bigger than maximum allowed.
	    -2 if hinted final is smaller than minimum allowed.
	    -3 if accepted initial value is greater or equal to final one.
	    -4 if both spacing and number of elements were hinted.
	    -5 if none or both of spacing or number of elements were in
               the referenced.
	    -6 Not valid oversampling was given by ref when requested.
*/
int makesample(prop_samp *samp,	/* Resulting sampled data */
	       prop_samp *hint,	/* Proposed sampling */
	       prop_samp *ref,	/* Reference values */
	       const long fl,
	       const int bitsshift,
	       const float margin)
{
  //'res' is the returned value
  //'n' and 'v' are auxiliary variables to produced sampling array
  int res=0,n;
  PREC_RES *v;
  double osd,si;

  //check initial value
  if(!(fl&(TRH_SI<<bitsshift))||hint->i<=0||hint->i<ref->i+margin){
    samp->i=ref->i+margin;
    res|=TRH_SI;
  }
  else if(hint->i>ref->f-margin){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Hinted initial value for %s sampling, is bigger than\n"
		 " maximum allowed final value %.8g. Consider margin %.8g\n"
		 ,TRH_SFTNAME(bitsshift),ref->f-margin,margin);
    return -1;
  }
  else
    transitaccepthint(samp->i,hint->i,fl,TRH_SI<<bitsshift);
  si=samp->i;

  //check final value
  if(!(fl&(TRH_SF<<bitsshift))||hint->f<=0||hint->f>ref->f-margin){
    samp->f=ref->f-margin;
    res|=TRH_SF;
  }
  else if(hint->f<ref->i+margin){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Hinted final value for %s sampling is smaller than\n"
		 "minimum allowed initial value %.8g.\n"
		 "Consider margin %.8g\n"
		 ,TRH_SFTNAME(bitsshift),ref->i+margin,margin);
    return -2;
  }
  else
    transitaccepthint(samp->f,hint->f,fl,TRH_SF<<bitsshift);

  //check that resultant range makes sense
  if(samp->f<=si){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Initial accepted sampling value (%g) is greater or\n"
		 "equal than final accepted sample value(%g).\n"
		 "%s was being hinted\n"
		 ,si,samp->f,TRH_SFTNAME(bitsshift));
    return -3;
  }

  transitprint(20,verblevel,
	       "Flags: 0x%lx    hint.d:%g   hint.n:%li\n"
	       ,fl,hint->d,hint->n);
  //check that only one of spacing or number of elements field have been
  //hinted
  if((fl&(TRH_SD<<bitsshift))&&(fl&(TRH_SN<<bitsshift))){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Both spacing(%g) and number of elements(%i) has\n"
		 "been hinted. That doesn't makes sense!. %s was being sampled.\n"
		 ,hint->d,hint->n,TRH_SFTNAME(bitsshift));
    return -4;
  }

  //if none has been hinted then use ref's
  if((!(fl&(TRH_SD<<bitsshift))||hint->d<=0)&&
     (!(fl&(TRH_SN<<bitsshift))||hint->n<=0)
     ){
    //If none or both of ref's exist then error
    if((ref->d<=0&&ref->n<=0)||(ref->d>0&&ref->n>0)){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "Spacing and number of elements were either both(%i) or\n"
		   "none(%i) in the reference for %s sampling. "
		   "And yes, none was hinted.\n"
		   ,ref->d>0&&ref->n>0,ref->d<=0&&ref->n<=0
		   ,TRH_SFTNAME(bitsshift));
      return -5;
    }
    //if spacing exists
    if(ref->d>0){
      samp->d=ref->d;
      samp->n=-1;
    }
    //otherwise use set array
    else{
      //If initial or final value were modified, then warn that the
      //array might be wrong!.
      if(res){
	/* TD: Check if they really changed */
	transiterror(TERR_WARNING,
		     "Array of length %i was given as reference\n"
		     "for %s sampling, but either (or both) the\n"
		     "initial(%g -> %g) and final(%g -> %g)\n"
		     "values might have been modified.\n"
		     ,ref->n,TRH_SFTNAME(bitsshift)
		     ,ref->i,si,ref->f,samp->f);
      }
      samp->n=ref->n;
      samp->d=-1;
      samp->v=(double *)calloc(ref->n,sizeof(double));
      memcpy(samp->v,ref->v,ref->n*sizeof(double));
      if(ref->o!=0)
	transiterror(TERR_WARNING,
		     "Fixed sampling array of length %i was referenced\n"
		     "But also oversampling was given (%g). Ignoring in\n"
		     "%s sampling\n"
		     ,samp->n,ref->o,TRH_SFTNAME(bitsshift));

      //return any possible modification
      return res;
    }
  }
  //else if spacing was hinted, then it has to be positive at this point
  else if(fl&(TRH_SD<<bitsshift)){
    transitASSERT(hint->d<=0,
		  "OOPS!, logic test 1 failed in %s's makesample()!!\n"
		  ,TRH_SFTNAME(bitsshift));
    transitaccepthint(samp->d,hint->d,fl,TRH_SD<<bitsshift);
  }
  //otherwise we have a positive hinted n
  else{
    transitASSERT(hint->n<=0,
		  "OOPS!, logic test 2 failed in %s's makesample()!!\n"
		  ,TRH_SFTNAME(bitsshift));
    transitaccepthint(samp->n,hint->n,*fl,TRH_SN<<bitsshift);
    samp->d=-1;
    samp->v=(double *)calloc(hint->n,sizeof(double));
    memcpy(samp->v,hint->v,hint->n*sizeof(double));
    if(hint->o!=0)
      transiterror(TERR_WARNING,
		   "Fixed sampling array of length %i was hinted\n"
		   "But also oversampling was given (%g). Ignoring in\n"
		   "%s sampling\n"
		   ,samp->n,hint->o,TRH_SFTNAME(bitsshift));
    return res;
  }

  //At this points make the sampling if it was not given as an
  //array.
  n=samp->n=(samp->f - si)/samp->d+1;

  //if there is an oversampling, check whether a value is not hinted
  if(!(fl&TRH_SO<<bitsshift)||hint->o<=0){
    //if so, check if we have a valid ref or error
    if(ref->o<=0){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "Not valid oversampling in the reference for\n"
		   "%s sampling.\n"
		   ,TRH_SFTNAME(bitsshift));
      return -6;
    }
    samp->o=ref->o;
  }
  else
    transitaccepthint(samp->o,hint->o,*fl,TRH_SO<<bitsshift);

  n=samp->n=(samp->n-1)*samp->o+1;
  osd=samp->d/(double)samp->o;

  //allocate and fill sampling array
  v=samp->v=(PREC_RES *)calloc(n,sizeof(PREC_RES));
  *v=si;
  v+=--n;
  while(n)
    *v--=si+n--*osd;
  //check the final point
  if(samp->v[samp->n-1]!=samp->f)
    transiterror(TERR_WARNING,
		 "Final sampled value (%g) of the\n"
		 "%li points doesn't coincide exactly with required\n"
		 "value (%g). %s sampling with pre-oversampling\n"
		 "spacing of %g.\n"
		 ,samp->v[samp->n-1],samp->n,samp->f
		 ,TRH_SFTNAME(bitsshift),samp->d);

  //return the flags of accepted values.
  return res;
}


/* \fcnfh
 Calls makesample with the appropiate parameters and set the flags

 @returns makesample() output
          -10 neither spacing nor number was hinted.
*/
int makewavsample(struct transit *tr)
{
  //'res' will be the result status
  int res;
  int fl=tr->ds.th->na;
  prop_samp *samp=&tr->ds.th->wavs;

  transitcheckcalled(tr->pi,"makewavsample",1,
		     "chkrange",TRPI_CHKRNG);

  if((!(fl&TRH_WD)||samp->d<=0)&&
     (!(fl&TRH_WN)||samp->n<=0)){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Spacing or number must be hinted for wavelength,\n"
		 "cannot just guess them.\n"
		 );
    return -10;
  }


  //make the sampling
  res=makesample(&tr->wavs,samp,&tr->ds.li->wavs,
		 fl,TRH_WAVSFT,tr->m);

  //set progress indicator if sampling was successful and return status
  if(res>=0)
    tr->pi|=TRPI_MAKEWAV;
  return res;
}


/* \fcnfh
 Calls makesample with the appropiate parameters and set the flags

 @returns makesample() output
*/
int makewnsample(struct transit *tr)
{
  //'res' will be the result status
  int res;
  prop_samp fromwav;
  memset(&fromwav,0,sizeof(prop_samp));
  struct transithint *trh=tr->ds.th;
  prop_samp *nsamp=&trh->wns;
  prop_samp *wsamp=&tr->wavs;

  transitcheckcalled(tr->pi,"makewnsample",1,
		     "makewavsample",TRPI_MAKEWAV);

  //convert from wavelength maximum
  fromwav.i=WNU_O_WLU/wsamp->f;

  //convert from wavelength minimum
  fromwav.f=WNU_O_WLU/wsamp->i;

  //set spacing such that the grid has the same number of points as the
  //wavelength one.
  fromwav.d=(fromwav.f-fromwav.i)*wsamp->o/wsamp->n;
  transitprint(20,verblevel,
	       "wavenumber spacing: %g "
	       "  comes from wavelength's number: %li\n"
	       "  in wavenumber range: %g to %g\n"
	       ,fromwav.d,wsamp->n,fromwav.i,fromwav.f);

  //use wavelength's oversampling
  fromwav.o=wsamp->o;

  //don't give a fixed array.
  fromwav.n=-1;

  //make the sampling
  res=makesample(&tr->wns,nsamp,&fromwav,trh->na,
		 TRH_WNSFT,tr->m*fromwav.f*fromwav.f/WNU_O_WLU);

  //set progress indicator if sampling was successful and return status
  if(res>=0)
    tr->pi|=TRPI_MAKEWN;
  return res;
}


/* \fcnfh
 Calls makesample with the appropiate parameters and set the flags

 @returns makesample() output
          1 Only one point value was requested
*/
int makeradsample(struct transit *tr)
{
  //'res' will be the result status
  //'nrad', 'niso' and 'neiso' will be the number of radius, regular and
  //extra isotopes points, respectively.
  //'limit' are limiting values that the sampling can take.
  //'rad', 'isovs', 'atms' and 'in' are auxiliary pointers to the radius
  //sampling, variable isotope info, atmosphere content, and variable
  //isotope structure pre-info, respectively.
  //'isov' and 'atmt' are structure pointers to where the info is going
  //to be stored after resampling.
  int res,i,j,iso1db;
  int nrad,niso=tr->n_i,neiso=tr->n_e,ndb=tr->n_db;

  prop_isov *isovs;

  struct atm_data *atms=tr->ds.at;
  prop_samp *rsamp=&atms->rads;

  struct iso_noext *in=tr->ds.in;
  prop_isov *isovt=tr->isov;
  prop_atm *atmt=&tr->atm;

  prop_samp *rad=&tr->rads;

  //getatm() and readinfo\_twii() must have been called before
  transitcheckcalled(tr->pi,"makeradsample",2,
		     "getatm",TRPI_GETATM,
		     "readinfo_twii",TRPI_READINFO);
  transitASSERT(atms->rads.n<1||!ndb||!neiso||!niso,
		"makeradsample():: called but essential variables are\n"
		"missing!\n");

  //We need to set-up limit so that the hinted values are compatible
  //with the atmosphere.
  //If there is only one atmospheric point then no sense in doing a
  //radius sampling
  if(rsamp->n==1){
    rad->f=rad->i=rsamp->v[0];
    nrad=rad->n=1;
    rad->d=-1;
    rad->v=(PREC_RES *)calloc(1,sizeof(PREC_RES));
    rad->v[0]=rsamp->v[0];
    //return succes as it would have makesample()
    res=0;
    /* TD: warn that hinted values are going to be useless */
  }
  //If there is more than one atmospheric point
  else{
    //If no initial value is hinted, we take atmosphere minimum.
    /* TD: more than one atmospheric point: completion of 'limit', don't
       forget to set nrad, set 'rsamp' from 'tr->ds.th.rads' if there
       are zeroes */
    
    nrad=rad->n=0;
    //do the sampling
    res=makesample(rad,&tr->ds.th->rads,rsamp,
		   tr->ds.th->na,TRH_RADSFT,0);
  }

  //Allocate arrays that will receive the interpolated data
  isovt->z=(PREC_ZREC *)calloc(nrad*neiso,sizeof(PREC_ZREC));
  isovt->d=(PREC_ATM *)calloc(nrad*neiso,sizeof(PREC_ATM));
  isovt->c=(PREC_CS *)calloc(nrad*neiso,sizeof(PREC_CS));
  for(i=0;i<neiso;i++){
    isovt[i].z=isovt->z+i*nrad;
    isovt[i].d=isovt->d+i*nrad;
    isovt[i].c=isovt->c+i*nrad;
  }
  atmt->t=(PREC_ATM *)calloc(nrad,sizeof(PREC_ATM));

  /* TD: interpolation */
  //interpolate temperature values according to radius
  resamplex(tr->fl,nrad,rad->v,rsamp->n,rsamp->v);
  resampley(tr->fl,1,atms->atm.t,atmt->t);

  //Now for the isotope. 
  //First, for the isotopes that were added by getatm() (extended). We
  //can use same x axis as for atmospheric sample, because they were set
  //in the same file!
  //Density for all the isotopes
  /* TD-BUG: Find out why it fails if I take the brackets away. */
  for(i=0;i<neiso;i++){
    resampley(tr->fl,1,atms->isov[i].d,isovt[i].d);
  }
  //and cross section of only the extended isotopes.
  for(i=niso;i<neiso;i++){
    resampley(tr->fl,1,atms->isov[i].c,isovt[i].c);
  }

  //Second, non-extended isotopes:
  //We have to go to each database separately
  for(i=0;i<ndb;i++){
    //position in the first isotope of the database
    iso1db=tr->db[i].s;
    isovs=in->isov+iso1db;

    //interpolate variable isotope info respect to temperature
    resamplex(tr->fl,in->db[i].t,in->db[i].T,nrad,atmt->t);
    for(j=0;j<tr->db[i].i;j++){
      transitASSERT(iso1db+j>neiso-1,
		    "trying to reference an isotope (%i) outside\n"
		    "the extended limit (%i)\n"
		    ,iso1db+j,neiso-1);
      resampley(tr->fl,2,isovs[j].z,isovt[iso1db+j].z,
		isovs[j].c,isovt[iso1db+j].c);
    }
  }

  if(res>=0)
    tr->pi|=TRPI_MAKERAD;
  return res;
}


/*\fcnfh
  extwn: Scattering parameters should be added at some point here

  @returns 0 on success
           -2 if no radii has been specified.
           -3 if only 2 wavelengths points
           -5 no isotopes selected!
	   -6 fine binning less than 1
	   -7 timesalpha less than 1
	   -8 maxratio less than 0
*/
int extwn (struct transit *tr)
{
  prop_samp *rad=&tr->rads;
  static struct extinction st_ex;
  tr->ds.ex=&st_ex;
  struct extinction *ex=&st_ex;
  struct line_transition *line=tr->lt;
  PREC_RES *k,**kiso,*wn,dwn,wavn,iniwn,wnmar,wni,wnf;
  PREC_NSAMP nrad,nwn;
  int neiso,niso;
  int r,i,ln;
  int w,*wa,subw;
  int j,maxj,minj,*nwnh;
  PREC_VOIGT ***profile, *profwn;
  PREC_VOIGTP *alphal,*alphad;
  int *cp,*wrc;
  double propto_adop,propto_alor, propto_k;
  PREC_ATM temp, *densiso;
  PREC_CS *csiso;
  PREC_ZREC *ziso;
  PREC_ZREC *mass;
  _Bool extinctperiso;

  transitcheckcalled(tr->pi,"kapwl",5,
		     "getatm",TRPI_GETATM,
		     "readinfo_twii",TRPI_READINFO,
		     "readdatarng",TRPI_READDATA,
		     "makewnsample",TRPI_MAKEWN,
		     "makeradsample",TRPI_MAKERAD
		     );
  transitacceptflag(tr->fl,tr->ds.th->fl,TRU_EXTBITS);

  //Check hinted values.
  //'.voigtfine' is the number of fine-bins of the Voigt function
  if(tr->ds.th->voigtfine<1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Fine binning of Voigt function  has to be\n"
		 "positive: %i\n"
		 ,tr->ds.th->voigtfine);
    return -6;
  }
  transitaccepthint(ex->vf,tr->ds.th->voigtfine,
		    tr->ds.th->na,TRH_VF);

  //'.timesalpha' is the number of alphas from the maximum of either
  //doppler or lorenz that the profile calculation have to consider.
  if(tr->ds.th->timesalpha<1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Times of maximum width has to be greater than\n"
		 "one: %i\n"
		 ,tr->ds.th->voigtfine);
    return -7;
  }
  transitaccepthint(ex->ta,tr->ds.th->timesalpha,
		    tr->ds.th->na,TRH_TA);

  //'.maxratio' is the maximum allowed ratio change before recalculating
  //profile array.
  if(tr->ds.th->maxratio_doppler<0){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Maximum allowed doppler width ratio change has to\n"
		 "be 0 or positive (%g)\n"
		 ,tr->ds.th->maxratio_doppler);
    return -8;
  }
  transitaccepthint(ex->maxratio,tr->ds.th->maxratio_doppler,
		    tr->ds.th->na,TRH_DR);
     

  iniwn=tr->wns.i;
  wn=tr->wns.v;
  wni=wn[0]-wnmar;
  wnf=wn[nwn-1]+wnmar;
  wnmar=tr->wnm=tr->m*wn[0]*wn[0];
  nwn=tr->wns.n;
  dwn=tr->wns.d/tr->wns.o;
  nrad=rad->n;
  neiso=tr->n_e;
  niso=tr->n_i;
  if(nrad<1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "There are no atmospheric parameters specified\n"
		 "I need at least one atmospheric point to calculate\n"
		 "a spectra");
    return -2;
  }
  if(nwn<2){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "I need at least 2 wavenumber points to compute\n"
		 "anything; I need resolution\n"
		 );
    return -3;
  }
  if(neiso<1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "You are requiring a spectra of zero isotopes!\n"
		 );
    return -5;
  }

  //following to store isotope dens
  densiso=(PREC_ATM *) calloc(niso,sizeof(PREC_ATM ));
  csiso=(PREC_CS *) calloc(niso,sizeof(PREC_CS));
  ziso=(PREC_ZREC *) calloc(niso,sizeof(PREC_ZREC));
  mass=(PREC_ZREC *) calloc(niso,sizeof(PREC_ZREC));
  nwnh=(int *) calloc(niso,sizeof(int));

  //allocate array for the voigt profile
  wa=(int *)calloc(niso,sizeof(int));
  profile=(PREC_VOIGT ***)calloc(niso,sizeof(PREC_VOIGT **));
  *profile=(PREC_VOIGT **)calloc(niso*ex->vf,sizeof(PREC_VOIGT *));

  //allocate indicative arrays. '.lw' stores info about line's width,
  //and '.cp' indicate when to change arrays.
  //'.recalc' index from which a recalculation of the voigt profile is
  //required.
  ex->lw=(int *)calloc(nwn,sizeof(int));
  wrc=(int *)calloc(niso,sizeof(int));
  ex->recalc=(int **)calloc(niso,sizeof(int *));
  ex->recalc[0]=(int *)calloc(nwn*niso,sizeof(int));
  for(i=1;i<niso;i++){
    ex->recalc[i]=ex->recalc[0]+i*nwn;
    profile[i]=*profile+i*ex->vf;
  }

  extinctperiso=(tr->fl&TRU_EXTINPERISO);
  //allocate array for extinctions and widths, initialization of second
  //index will be done in the loop (\lin{kini})
  ex->k=(PREC_RES ***)calloc(nrad,sizeof(PREC_RES **));
  ex->al=(PREC_VOIGTP **)calloc(nrad,sizeof(PREC_VOIGTP *));
  ex->ad=(PREC_VOIGTP **)calloc(nrad,sizeof(PREC_VOIGTP *));
  kiso=*ex->k=(PREC_RES **)calloc(niso*nrad,sizeof(PREC_RES *));
  i=extinctperiso?niso:1;
  ex->k[0][0]=(PREC_RES *)calloc(nrad*i*nwn,sizeof(PREC_RES));
  *ex->al=(PREC_VOIGTP *)calloc(nrad*niso,sizeof(PREC_VOIGTP));
  *ex->ad=(PREC_VOIGTP *)calloc(nrad*niso,sizeof(PREC_VOIGTP));

  //For each radius (index 'r')
  for(r=0;r<nrad;r++){

    transitprint(2,verblevel,"Radius %i: %g[planetary radius]\n",r,rad->v[r]);

    //Initialization of 2nd dimension of extinction array.
    //\linelabel{kini}
    i=extinctperiso?niso:1;
    kiso=ex->k[r]=*ex->k+r*i*nwn;
    alphal=ex->al[r]=*ex->al+r*niso;
    alphad=ex->ad[r]=*ex->ad+r*niso;

    //set some auxiliary variables.
    temp=tr->atm.t[r];

    //'propto\_adop' is proportional to the doppler width, which in its
    //total splendor is
    //\[
    //\alpha_D=\frac{\Wn}{\sqrt{m}}\underbrace{\frac{\sqrt{2k_BT}}{c}}
    //_{\mathrm{propto\_adop}}
    //\label{dopbr}
    //\]
    propto_adop=sqrt(2*KB*temp)/LS;

    //'propto\_alor' is proportional to the Lorenz width, which in its
    //total splendor is
    //\[
    //\alpha_L=\frac{\sigma_c}{\pi c}
    //\sqrt{\frac{2\kb T}{\pi}} \sum_{\mathrm{collisioners}}n_i
    //\sqrt{\left(\frac{1}{m_r}
    //+\frac{1}{m_i}\right)}
    //\label{lorwidth}
    //\]
    propto_alor=sqrt(temp*2*KB/PI)/LS/PI;

    //Initialize a voigt profile for every isotope as well for the
    //mass, ziso, densiso and csiso arrays
    for(i=0;i<niso;i++){
      kiso[i]=kiso[0];
      if(extinctperiso)
	kiso[i]+=nwn*i;
      mass[i]=tr->isof[i].m;
      ziso[i]=tr->isov[i].z[r];
      densiso[i]=tr->isov[i].d[r];
      csiso[i]=tr->isov[i].c[r];

      //Calculate lorentz
      alphal[i]=0;
      for(j=0;j<niso;j++)
	alphal[i]+=densiso[i]*csiso[i]/mass[i]
	  *sqrt(1/mass[i] + 1/tr->isof[j].m);
      alphal[i]*=propto_alor;

      //Following calculates doppler divided by central wavenumber.
      alphad[i]=propto_adop/sqrt(mass[i]);

      if((nwnh[i]=newprofile(profile[i],ex->vf,&ex->lw[0],dwn,
			     wn[0]*alphad[i],alphal[i],ex->ta)
	  )<1)
	transiterror(TERR_CRITICAL,
		     "newprofile() returned error code %i on its\n"
		     "first try for isotope %i.\n"
		     ,nwnh[i],i);
      w=nwn-1;
      cp=ex->recalc[i];
      cp[w]=(int)(ex->maxratio*wn[w]/dwn+0.5);
      if(!cp[w])
	cp[w]=1;
      wrc[i]=w-cp[w];
      if(wrc[i]>=0)
	cp[wrc[i]]=1;
      wa[i]=w;
    }

    for(ln=0;ln<tr->n_l;ln++){
      /*
      if(ln!=10000&&ln!=10702&&ln!=10402)
	continue;
      if(ln<9000||ln>11000)
	continue;
      */

      wavn=WNU_O_WLU/line[ln].wl;
      /* 
       if(wavn<wni||wavn>wnf)
       continue;
      */
      //if it is beyond lower limit
      if(wavn<iniwn)
	continue;
      //      w=-(long)((iniwn-wavn)/dwn+1)
      else
	w=(wavn-iniwn)/dwn;
      transitDEBUG(20,verblevel,
		   "wavn:%g lgf:%g\n"
		   ,wavn,line[ln].lgf);
      //If it is beyond the last then just
      //skip that line
      if(w>=nwn)
	continue;

      subw=ex->vf*(wavn-w*dwn-iniwn)/dwn;
      i=line[ln].isoid;
      k=kiso[i];

      cp=ex->recalc[i];

      transitASSERT(wa[i]!=-1&&wa[i]<w,
		    "Database is not ordered!, previous wavenumber was\n"
		    "at index %i, new one at %i (it should have been smaller)\n"
		    ,wa,w);
      if(w<=wrc[i]){
	//Find number of wavenumbers until the next recalculation
	cp[wrc[i]]=0;
	cp[w]=(int)(ex->maxratio*wn[w]/dwn+0.5);
	if(!cp[w])
	  cp[w]=1;
	wrc[i]=w-cp[w];
	if(wrc[i]>=0)
	  cp[wrc[i]]=1;
	transitDEBUG(22,verblevel,
		     "Recalculating voigt for isotope %i... current\n"
		     "wavenumber %i, next wavenumber %i/%i\n"
		     ,i,w,wrc[i],nwn);

	free(profile[i]);

	if((nwnh[i]=newprofile(profile[i],ex->vf,&ex->lw[w],dwn,
			       wn[w]*alphad[i],alphal[i],ex->ta)
	    )<1)
	  transiterror(TERR_CRITICAL,
		       "newprofile() returned error code %i for\n"
		       "isotope %i\n"
		       ,nwnh[i],i);
      }

      /* CAVEATS: _mass_ densitty 
                  _log_ gf */
      propto_k=densiso[i]	         //mass density
	*SIGCTE			         //Constant in sigma
	*line[ln].lgf		         //Log(gf)
	*exp(-EXPCTE*line[ln].elow/temp) //Level population
	*(1-exp(-EXPCTE*wavn/temp))      //induced emission
	/mass[i]	        	 //mass
	/ziso[i];		         //Partition function

      transitDEBUG(20,verblevel,
		   "i=%i   temp=%g   Elow=%g\n"
		   "k= %10.3g  //densiso[i] \n"
		   "  *%10.3g  //SIGCTE\n"
		   "  *%10.3g  //line[ln].lgf\n"
		   "  *%10.3g  //exp(-EXPCTE*line[ln].elow/temp)\n"
		   "  *%10.3g  //(1-exp(-EXPCTE*wavn/temp))\n"
		   "  /%10.3g  //mass[i]\n"
		   "  /%10.3g  //ziso[i]\n"
		   " = %10.3g   //extinction\n"
		   ,i,temp,line[ln].elow
		   ,densiso[i]
		   ,SIGCTE
		   ,line[ln].lgf
		   ,exp(-EXPCTE*line[ln].elow/temp)
		   ,(1-exp(-EXPCTE*wavn/temp))
		   ,mass[i]
		   ,ziso[i]
		   ,propto_k);

      //set 'profwn' such that the index mimic wavenumber's array
      profwn=profile[i][subw]+nwnh[i]-w;

      //set upper and lower limits.
      minj=w-nwnh[i];
      if(minj<0)
	minj=0;
      maxj=w+nwnh[i];
      if(maxj>=nwn)
	maxj=nwn-1;

      //distribute the oscillator strength according to the voigt
      //profile
      for(j=minj;j<maxj;j++)
	k[j]+=propto_k
	  *profwn[j];


      wa[i]=w;
    }
    for(i=0;i<niso;i++)
      free(profile[i][0]);
  }

  return 0;
}

/* 
   calculates a new voigt profile

   @returns number of points to center wavelength
*/
inline int newprofile(PREC_VOIGT **pr, /* output 2d profile */
		      int vf, /* 1st already allocated dimension */
		      int *lw,	/* where width of the profile (2nd
				   dimension) is to be stored. */
		      PREC_RES dwn, /* waavenumber spacing */
		      PREC_VOIGT dop, /* doppler width */
		      PREC_VOIGT lor, /* lorenz width */
		      float ta)	/* times of alpha */
{
  PREC_VOIGTP bigalpha;
  PREC_VOIGTP wvgt;
  int nvgt,j;

  //look for the biggest alpha (Uses Doppler width given by the
  //second index)
  bigalpha=dop;
  if(bigalpha<lor)
    bigalpha=lor;

  //Find the width of the profile such that at least 'op->ta' widths
  //are present on it.
  //Afterwards set the number of points in the array, an odd number
  //is preferred so that the central bin can have a maximum value.
  wvgt=bigalpha*ta;
  *lw=nvgt=2*wvgt/dwn+1;

  //Initialize array that will hold the profile.
  *pr=(PREC_VOIGT *)calloc(nvgt*vf,sizeof(PREC_VOIGT));
  for(j=1;j<vf;j++)
    pr[j]=pr[0]+j*nvgt;

  //calculate voigt
  if((j=voigtn(vf, nvgt, wvgt, lor,dop,pr, -1, 0))!=1)
    transiterror(TERR_CRITICAL,
		 "voigtn() returned error code %i\n"
		 ,j);

  return nvgt/2;
}


//\delfh
#if 0
/* TD: change 'will be' and implement methods */

/*\fcnfh
  telresconv: Convolve to telescope resolution, several 'mode'(int) will
  be available.

  @returns 0  on success or width == 0
           -1 if mode is an unaccepted value
	   -2 width array for TRC\_GIVEN is not the same length as 'nwl'
	   -3 wavelength array too small or not sorted.
	   -4 'osw' less than 1
*/
int telresconv(PREC_RES **spectra, //pointer to resulting spectra
	       int *nsp,           //Number of new points
	       PREC_RES *presp,    //Pre-convol spectra
	       PREC_RES *wl,       //wavelength points
	       int nwl,            //Number of data points for presp
	       int osw,            //Oversampling
	       int mode,           //what kind of convolution
	       ...)                //depends on the mode selected
{

  va_list ap;
  double w1,w2,*w;
  float mgau;
  int nw;
  int i,k;   //GO through presp, spectra
  double *smr,dl,width,tmp;
  int nsmr;
  int idxi,idxf,idx;

  /* TD: Give the user the optioni to change the following default */
  mgau=5;

  if(nwl<2||wl[1]<wl[0]){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "I need at least 2 wavelengths points (%i given) and\n"
		 "they have to be sorted. telresconv()\n"
		 ,nwl);
    return -3;
  }
  dl=wl[1]-wl[0];

  *spectra=(PREC_RES *)calloc(nwl,sizeof(PREC_RES));
  *nsp=(int)(tmp=nwl/osw);
  if(*nsp!=tmp)
    transiterror(TERR_WARNING,
		 "The total number of wavelength points(%i) before\n"
		 "telescope resolution convolution was not an integer\n"
		 "times the oversampling(%i).\n"
		 ,nwl,osw);

  va_start(ap,mode);
  switch(mode&TRC_MODEBITS){
  case TRC_FFT:
    transiterror(TERR_CRITICAL,
		 "Requested mode (%xi) has not yet being implemented\n"
		 "in telresconv()\n\n"
		 ,mode&TRC_MODEBITS);
    return -1000;
    /* TD: a clean exit with width ==0 */
    break;
  default:
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Wrong requested mode (%xi) for convolution of\n"
		 "Telescope resolution\n\n"
		 ,mode&TRC_MODEBITS);
    return -1;
    break;
  case TRC_LINEAR:
    w1=va_arg(ap,double);
    if(w1==0){
      free(*spectra);
      *spectra=presp;
      *nsp=nwl;
      return 0;
    }
    w2=va_arg(ap,double);
    nw=nwl;
    tmp=(w2-w1)/(nw-1);
    w=(double *)calloc(nw,sizeof(double));
    for(i=0;i<nw;i++)
      w[i]=w1+i*tmp;
    break;
  case TRC_GIVEN:
    w=va_arg(ap,double *);
    if(w==NULL){
      free(*spectra);
      *spectra=presp;
      *nsp=nwl;
      return 0;
    }
    nw=va_arg(ap,int);
    if (nw!=*nsp){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "Number of elements of the width array(%i)\n"
		   "for convolution in mode 'TRC_GIVEN' is not the\n"
		   "same as the number of wavelength elements(%i)\n\n"
		   ,nw,nwl);
      return -2;
      break;
    }
    w1=0;               //just to prevent 'uninitialized' warnings 
    break;
  case TRC_FIX:
    w1=va_arg(ap,double);
    if(w1==0){
      free(*spectra);
      *spectra=presp;
      *nsp=nwl;
      return 0;
    }
    nw=0;
    w=NULL;             //just to prevent 'uninitialized' warnings 
    break;
  }

  va_end(ap);
  width=0;
  smr=NULL;             //just to prevent 'uninitialized' warnings 
  nsmr=0;               //just to prevent 'uninitialized' warnings 

  for(i=0;i<nwl;i++){
    if((i%osw)==0){
      k=(int)(i/osw);
      do{
	if(nw){
	  width=w[k];
	  if(i)
	    free(smr);
	}
	else{
	  if(i)
	    break;
	  width=w1;
	}
	switch(mode&TRC_METHODBITS){
	default:
	  transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		       "Wrong requested method (%xi) for convolution of\n"
		       "telescope resolution\n\n"
		       ,mode&TRC_METHODBITS);
	  return -1;
	  break;
	case TRC_GAUSSIAN:
	  makegauTR(&nsmr, &smr,width,dl,mgau);
	  break;
	case TRC_BOX:
	  nsmr=(int)(width/dl);
	  nsmr+=nsmr&1?0:1;
	  smr=(double *)calloc(nsmr,sizeof(double));
	  for(i=0;i<nsmr;i++){
	    smr[i]=1.0/(double)nsmr;
	  }
	}
	
      }while(0);
    }

    idxi=i-nsmr/2;
    idx=0;
    if(idxi<0){
      idx=-idxi;
      idxi=0;
    }

    idxf=i+nsmr/2;
    if(idxf>nwl)
      idxf=nwl;

    for(;idxi<idxf;idxi++){
      (*spectra)[idxi]+=presp[i]*smr[idx];
      transitDEBUG(22,verblevel,
		   "for %i: sp %i: gau %i\t\tpre:%g post:%g\n"
		   ,i,idxi,idx,presp[i],(*spectra)[idxi]);
      idx++;
    }

  }

  if(osw<1){
    transiterror(TERR_SERIOUS,
		 "Oversampling has to be greater or equal to 1 (%i)\n\n"
		 ,osw);
    return -4;
  }
  else if (osw>1){
    for(i=0;i<*nsp;i++)
      (*spectra)[i]=integrateTR(*spectra+i*osw,osw)*dl;

    *spectra=(PREC_RES *)realloc(*spectra,*nsp*sizeof(PREC_RES));
  }

  return 0;
}

/*\fcnfh
  Private function that returns a gaussian shape.
 */
static inline int makegauTR(int *nsmr,
			    double **smr,
			    double width,
			    double dl,
			    float mgau)
{
  int med;
  double *mark,expo;
  int n;

  *nsmr=(int)(width*mgau/dl);
  *nsmr+=(*nsmr&1)==0;
  n=*nsmr;
  mark=*smr=(double *)calloc(*nsmr,sizeof(double));

  med=*nsmr/2+1;

  while(n){
    expo=(n-med)/width;
    *(mark++)=(ONEOSQRT2PI)/width*exp(-expo*expo*0.5);
    transitDEBUG(22,verblevel,"g: %3i:  %f\n",n,*(mark-1));
    n--;
  }

  return 0;
}

/*\fcnfh
  Private function that returns the quasi-integrated value of an array
  of given length, assuming that the spacing is unity. If it is not,
  just multiply result by bin width.
*/
static inline PREC_RES integrateTR(PREC_RES *array,
				   int length)
{

  return 0;
  
}

/*\fcnfh
  Given X and Y-arrays interpolate for x-position val.
 */
static inline PREC_ZREC interpolateTR(PREC_ZREC *x,
			       PREC_ZREC *y,
			       int n,
			       PREC_ATM val)
{
  int i,f,new;


  i=0;
  f=n-1;

  if(*x>val||x[f]<val)
    transiterror(TERR_SERIOUS,
		 "interpolateTR(), looked for value (%f) is not\n"
		 "within the array of paramaters (%f - %f) and I\n"
		 "cannot extrapolate\n"
		 ,val,*x,x[n-1]);

  do{
    new=(i+f)/2;
    if(val>x[new])
      i=new;
    else
      f=new;
  }while(f-i>1&&x[i+1]<val);

  transitDEBUG(22,verblevel,"Found value %f between (%f,%f) and (%f,%f)\n"
	       "Returned value is %f\n",val,x[i],y[i],x[i+1],y[i+1],
	       y[i]+(val-x[i])*(y[i+1]-y[i])/(x[i+1]-x[i]));

  return y[i]+(val-x[i])*(y[i+1]-y[i])/(x[i+1]-x[i]);

}
#endif
//\delfh

#if 0  /* Alternative obsolete in kapwl() */
    lwl=uwl=0;

    /* Fill for each wavelength */
    for(j=0;j<nwl;j++){

      transitdot(3,verblevel);
      transitDEBUG(22,verblevel,"NWL: %i %i %f\n",j,nwl,(wl)[4254]);
      wl2=wl[j]*wl[j];
      if(wl2/lastdwn2>deltaw&&j<nwl-1){
	deltaw=(wl[j+1]-wl[j]);
	deltawn=WNU_O_WLU*deltaw/wl[j]/wl[j];
	lastdwn2=wl[j]*wl[j];
	wnwidth=alor>adopprop?alor*prtimesw:adopprop*prtimesw;
	wlhwidth=lastdwn2*wnwidth/2.0/WNU_O_WLU;
	nwn=(int)(wnwidth/deltawn);
	vpro=(PREC_VOIGT *)realloc(vpro,nwn*sizeof(PREC_VOIGT));
	if((rc=voigtn(nwn, deltawn, alor,adopprop,vpro, eps))!=1){
	  transiterror(TERR_CRITICAL,
		       "voigtn() returned error code %i in kapwl().\n");
	}
      }	

      while(line[lwl].wl<wl[j]-wlhwidth)
	lwl++;
      while(uwl<nlines&&line[uwl].wl<wl[j]+wlhwidth)
	uwl++;

      if(uwl<=lwl){
	transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
		     "There are no lines from which to make the\n"
		     "spectra. Invalid range %i to %i was selected\n"
		     "out of maimum %i\n"
		     ,lwl,uwl,nwl);
	return -4;
      }

      kaptmp=0;
      pos=0;
      for(k=lwl;k<uwl;k++){
	pos=nwn/2-(int)((line[k].wl-wl[j])*WNU_O_WLU/wl2/deltawn+0.5);
	if(pos<0)
	  pos=0;
	if(pos>=nwn)
	  pos=nwn-1;
	transitDEBUG(24,verblevel,
		     "In wavelength %f, we found a distance of %i\n"
		     "units of wavenumber (%f) (as opposed to wavelength\n"
		     "%f) to the center of the line %f.\n"
		     ,wl[j],nwn/2-pos,deltawn,wl2*deltawn/WNU_O_WLU,
		     line[k].wl);
	rc=line[k].isoid;
	kaptmp+=
	  *densisorad[i]	/* Density */
	  *vpro[pos]		/* Line Profile */
	  *(1-exp(-EXPCTE*(wl[j]/Trad[i]))) /* Induced emission */
	  /* Following are factors within sigma */
	  *SIGCTE	        /* Constant */
	  *exp(line[k].lgf)	/* gf: degeneracy*osc.strength */
	  *exp(-line[k].elow/(KB)/Trad[i])  /* Level population */
	  /mmass[rc]            /* Isotope mass */
	  /Zisorad[rc][i];	/* Partition function */
      }
      (*kap)[i][j]=kaptmp;
      //%      transitDEBUG(20,verblevel,"NW: %f\n",wl[4254]);
      transitDEBUG(24,verblevel,
		   "%i/%i\t%f\t%f linmax:%i/%i line:%i/%i\n"
		   ,j,nwl,wl[j],(*kap)[i][j],pos,i,uwl,nlines);
      //%      break;
    }
    transitDEBUG(20,verblevel,"Next radius\n");
#endif /* obsolete */

#if 0  /* Integration & convolution at once in telresconv()... too
	  cumbersome  */
    idxt=nwl-i+nsmr/2;
    if(idxt>nsmr)
      idxt=nsmr;

    idxi=i-nsmr/2;
    idsi=0;
    if(idxi<0){
      idsi=-idxi;
      idxt+=idxi;
      idxi=0;
    }

    idib=idxi%osw;
    fdob=(idxi+idxt)/osw;
    for(idob=idxi/osw;idob<fdob;idob++){
      (*spectra)[idob]=convinteg(smr+i dsi,presp+idob*osw+idib,osw-idib);
      idsi+=osw-idib;
      idib=0;
    }
    if(fdob==idxi/osw)
      (*spectra)[idob]=convinteg(smr+idsi,presp+idob*osw+idib,idxt);
    else
      (*spectra)[idob]=convinteg(smr+idsi,presp+idob*osw+idib,
				 (idxi+idxt)%osw);
#endif /* obsolete */
//\deluh
