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
#include <util/procopt.h>
#include <math.h>

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
  //'f\_line' is the name of the twii data file 
  //'trh.m' is the amount of microns that cannot be trusted at the
  //boundaries of databases range, and also how much extra out of
  //requested range it has to look for transitions.
  //'defile\_line' is the default name of the line info file.
  //(i.e. modifiable by the user)
  trh.m=0.001;
  trh.na|=TRH_WM;
  char defile_line[]="./res/lineread.twii";
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
  if((rn=processparameters(argc,argv,&trh))!=0)
    transiterror(TERR_SERIOUS,
		 "processparameters() returned error code %i\n"
		 ,rn);

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
    /*    if(rn%transit.wns.o==0)*/
      printf("%10.4f%10.4f%15.5g%15.5g\n"
	     ,transit.wns.v[rn],WNU_O_WLU/transit.wns.v[rn],
	     transit.ds.ex->k[0][0][rn],
	     AMU*transit.ds.ex->k[0][0][rn]*transit.isof[0].m/transit.isov[0].d[0]);




  return EXIT_SUCCESS;
}


/* \fcnfh
   process command line options, saving them in the hint structure.

   @returns 0 on success
 */
int processparameters(int argc, /* number of command line arguments */
		      char **argv, /* command line arguments */
		      struct transithint *hints) /* structure to store
						    hinted parameters */
{
  //different non short options
  enum param {
    CLA_DUMMY=128,
    CLA_ATMOSPHERE,
    CLA_LINEDB,
    CLA_RADLOW,
    CLA_RADHIGH,
    CLA_RADDELT,
    CLA_WAVLOW,
    CLA_WAVHIGH,
    CLA_WAVDELT,
    CLA_WAVOSAMP,
    CLA_WAVMARGIN,
    CLA_WAVNLOW,
    CLA_WAVNHIGH,
    CLA_WAVNDELT,
    CLA_WAVNOSAMP,
    CLA_WAVNMARGIN,
    CLA_ONEPT,
    CLA_ONEABUND,
    CLA_ONEINT,
    CLA_ONEEXTRA,
  };

  //General help-option structure
  struct optdocs var_docs[]={
    {NULL,HELPTITLE,0,
     NULL,"GENERAL ARGUMENTS"},
    {"version",no_argument,'V',
     NULL,"Prints version number and exit"},
    {"help",no_argument,'h',
     NULL,"Prints list of possible parameters"},
    {"defaults",no_argument,'d',
     NULL,"Prints default values of the different variable"},
    {"verb",required_argument,'v',
     "[+..][-..]","Increase or decrease verbose level by one per\n"
     "each + or -. 0 is the quietest"},

    {NULL,HELPTITLE,0,
     NULL,"INPUT/OUTPUT"},
    {"output",required_argument,'o',
     "outfile","Change output file name, a dash (-)\n"
     "directs to standard output"},
    {"atmosp",required_argument,CLA_ATMOSPHERE,
     "atmfile","File containing atmospheric info (Radius,\n"
     "pressure, temperature). A dash (-) indicates alternative\n"
     "input"},
    {"linedb",required_argument,CLA_LINEDB,
     "linedb","File containing line information (TWII format,\n"
     "as given by 'lineread'"},

    {NULL,HELPTITLE,0,
     NULL,"RADIUS OPTIONS (planetary radii units, unless stated "
     "otherwise)"},
    {"radius",no_argument,'r',
     NULL,"Interactively input radius parameters"},
    {"rad-low",required_argument,CLA_RADLOW,
     "radius","Lower radius. 0 if you want to use atmospheric\n"
     "data minimum"},
    {"rad-high",required_argument,CLA_RADHIGH,
     "radius","Higher radius. 0 if you want to use atmospheric\n"
     "data maximum"},
    {"rad-delt",required_argument,CLA_RADDELT,
     "spacing","Radius spacing. 0 if you want to use atmospheric\n"
     "data spacing"},
    {"oneptn",required_argument,CLA_ONEPT,
     "press,temp,extra_iso","Don't calculate transit spectra, just\n"
     "obtain spectra for a given pressure and temperature. Unless\n"
     "oneabund is also specified and has the correct number of\n"
     "isotopes, the abundances will be asked interactively"},
    {"oneextra",required_argument,CLA_ONEEXTRA,
     "mass1name1,mass2name2,...","It only has effect with --onept,\n"
     "a list of the atomic mass and names for the hitherto specified\n"
     "extra isotopes. If it doesn't have the right amount of values,\n"
     "the program will ask interactively.\n"},
    {"oneabund",required_argument,CLA_ONEABUND,
     "q1,...","It also only has effect with --onept, a list of the\n"
     "abundances of the different isotopes. If it is omitted or\n"
     "doesn't have the right amount of values, the program will\n"
     "ask interactively. Note that the order of isotopes is the\n"
     "same given in the TWII data file"},
    {"onept-interactive",no_argument,CLA_ONEINT,
     NULL,"Wants to give abundances and pressure and temperature\n"
     "interactively through terminal input"},

    {NULL,HELPTITLE,0,
     NULL,"WAVELENGTH OPTIONS (all in nanometers)"},
    {"wavelength",no_argument,'w',
     NULL,"Interactively input wavelength parameters"},
    {"wl-low",required_argument,CLA_WAVLOW,
     "wavel","Lower wavelength. 0 if you want to use line\n"
     "data minimum"},
    {"wl-high",required_argument,CLA_WAVHIGH,
     "wavel","Upper wavelength. 0 if you want to use line\n"
     "data maximum"},
    {"wl-delt",required_argument,CLA_WAVDELT,
     "spacing","Wavelength spacing. 0 if you want to use line\n"
     "data spacing"},
    {"wl-osamp",required_argument,CLA_WAVOSAMP,
     "integer","Wavelength oversampling"},
    {"wl-marg",required_argument,CLA_WAVMARGIN,
     "boundary","Not trustable range in microns at boundary\n"
     "of line databases. Also transitions this much away from\n"
     "the requested range will be considered"},

    {NULL,HELPTITLE,0,
     NULL,"WAVENUMBER OPTIONS (all in cm-1)"},
    {"wavenumber",no_argument,'n',
     NULL,"Interactively input wavenumber parameters"},
    {"wn-low",required_argument,CLA_WAVNLOW,
     "waven","Lower wavenumber. 0 if you want to use equivalent\n"
     "of the wavelength maximum"},
    {"wn-high",required_argument,CLA_WAVNHIGH,
     "waven","Upper wavenumber. 0 if you want to use equivalent\n"
     "of the wavelength minimum"},
    {"wn-delt",required_argument,CLA_WAVNDELT,
     "spacing","Wavenumber spacing. 0 if you want to have the\n"
     "same number of points as in the wavelength sampling"},
    {"wn-osamp",required_argument,CLA_WAVNOSAMP,
     "integer","Wavenumber oversampling. 0 if you want the same\n"
     "value as for the wavelengths"},
    {"wn-marg",required_argument,CLA_WAVNMARGIN,
     "boundary","Not trustable range in cm-1 at boundaries.\n"
     "Transitions this much away from the requested range will\n"
     "be considered"},

    {NULL,HELPTITLE,0,
     NULL,"OPACITY CALCULATION OPTIONS:"},
    {"finebin",required_argument,'f',
     "integer","Number of fine-bins to calculate the Voigt\n"
     "function"},
    {"nwidth",required_argument,'a',
     "number","Number of the max-widths (the greater of Voigt\n"
     "or Doppler widths) that need to be contained in a calculated\n"
     "Voigt profile"},
    {"maxratio",required_argument,'u',
     "uncert","Maximum allowed uncertainty in doppler width before\n"
     "recalculating profile"},

    {NULL,HELPTITLE,0,
     NULL,"OBSERVATIONAL OPTIONS:"},
    {"telres",required_argument,'t',
     "width","Gaussian width of telescope resolution in nm"},

    {NULL,0,0,NULL,NULL}
  };

  struct optcfg var_cfg={NULL,NULL,NULL,
			 "Patricio Rojo <pato@astro.cornell.edu>"
			 ,NULL,NULL,0};
  int rn,i;
  //  int longidx;
  prop_samp *samp;
  char name[20],rc;
  char *sampv[]={"Initial","Final","Spacing","Oversampling integer for"};
  double rf;

  opterr=0;
  while(1){
    /* This is for old style
       rn=getopt(argc,argv,"f:Vhv:m:r:w:n:a:s:d:");*/
    rn=getprocopt(argc,argv,var_docs,&var_cfg,NULL);
    if (rn==-1)
      break;

    transitDEBUG(20,verblevel,
		 "Processing option '%c', argum: %s\n"
		 ,rn,optarg);

    switch(rn){
    case CLA_ATMOSPHERE:	//Atmosphere file
      hints->f_atm=(char *)realloc(hints->f_atm,strlen(optarg)+1);
      strcpy(hints->f_atm,optarg);
      break;
    case CLA_LINEDB:		//line database file
	hints->f_line=(char *)realloc(hints->f_line,strlen(optarg)+1);
	strcpy(hints->f_line,optarg);
	break;
    case 'o':			//output file
      hints->f_out=(char *)realloc(hints->f_out,strlen(optarg)+1);
      strcpy(hints->f_out,optarg);
      break;

    case 'r':
      samp=&hints->rads;
      fprintpad(stderr,1,"In units of planetary radius...\n");
      strcpy(name,"radius");
    case 'w':
      if(rn=='w'){
	samp=&hints->wavs;
	fprintpad(stderr,1,"In nanometers units...\n");
	strcpy(name,"wavelength");
      }
    case 'n':
      if(rn=='n'){
	samp=&hints->wns;
	fprintpad(stderr,1,"In cm-1 units...\n");
	strcpy(name,"wavenumber");
      }

      for(i=0;i<4;i++){
	if(i==3&&samp==&hints->rads)
	  break;
	while(rn){
	  fprintf(stderr,"- %s %s: ",sampv[i],name);
	  samp->i=readd(stdin,&rc);
	  switch(rc){
	  case 0:
	    break;
	    /*
	      case '?':
	      fprintpad(stderr,1,"%s\n",var_docs[longidx].doc);
	    */
	  default:
	    fprintf(stderr,"\nTry again\n");
	  }
	}
	rn=1;
	/*	longidx++;*/
      }
      break;

    case CLA_ONEPT:
      if((rn=getnd(3,',',optarg,&hints->onept.p,&hints->onept.t,
		   &rf))!=3){
	if(rn>0)
	  transiterror(TERR_SERIOUS,
		       "At least one of the values given for pressure (%g),\n"
		       "temperature (%g), or number of extra isotopes (%g),\n"
		       "was not a correct floating point value. Actually, the\n"
		       "latter should be an integer\n"
		       ,hints->onept.p,hints->onept.t,rf);
	else
	  transiterror(TERR_SERIOUS,
		       "There was %i comma-separated fields instead of 3 for\n"
		       "'--onept' option"
		       ,-rn);
      }
      hints->onept.ne=(int)rf;
      if(rf!=hints->onept.ne)
	transiterror(TERR_SERIOUS,
		     "A non integer(%g) number of extra isotopes was given with\n"
		     "the option --oneptn\n"
		     ,rf);
      hints->onept.one=1;
      hints->fl=(hints->fl&~TRU_ATM1PBITS)|TRU_ATMGIVEN1P;
      break;
    case CLA_ONEABUND:
      if((hints->onept.nq=getad(0,',',optarg,&hints->onept.q))<1)
	transiterror(TERR_SERIOUS,
		     "None of the given isotope abundances were accepted\n"
		     "%s\n"
		     ,optarg);
      transitprint(2,verblevel,
		   "%i abundance isotopes were correctly given: %s"
		   ,hints->onept.nq,optarg);
      break;

    case CLA_ONEEXTRA:
      /* TD from here: extra isotope name and mass as parameters */
      break;

    case CLA_ONEINT:
      hints->fl=(hints->fl&~TRU_ATM1PBITS)|TRU_ATMASK1P;
      break;
    case CLA_RADLOW:
      hints->rads.i=atof(optarg);
      break;
    case CLA_RADHIGH:
      hints->rads.f=atof(optarg);
      break;
    case CLA_RADDELT:
      hints->rads.d=atof(optarg);
      break;

    case CLA_WAVLOW:
      hints->wavs.i=atof(optarg);
      break;
    case CLA_WAVHIGH:
      hints->wavs.f=atof(optarg);
      break;
    case CLA_WAVDELT:
      hints->wavs.d=atof(optarg);
      hints->wavs.n=0;
      hints->wavs.v=NULL;
      break;
    case CLA_WAVOSAMP:
      hints->wavs.o=atof(optarg);
      break;
    case CLA_WAVMARGIN:
      hints->m=atof(optarg);
      break;
    case CLA_WAVNLOW:
      hints->wavs.i=atof(optarg);
      break;
    case CLA_WAVNHIGH:
      hints->wavs.f=atof(optarg);
      break;
    case CLA_WAVNDELT:
      hints->wavs.d=atof(optarg);
      hints->wavs.n=0;
      hints->wavs.v=NULL;
      break;
    case CLA_WAVNOSAMP:
      hints->wavs.o=atof(optarg);
      break;
    case CLA_WAVNMARGIN:
      hints->wnm=atof(optarg);
      break;
    case 't':			//Telescope resolution
      hints->t=atof(optarg);
      break;

    case 'u':			//Change Doppler's maximum accepted
				//ratio before recalculating
      hints->maxratio_doppler=atof(optarg);
      break;

    case 'f':			//Change voigt fine-binning
      hints->voigtfine=atoi(optarg);
      break;

    case 'a':			//Change times of alphas in profile
      hints->timesalpha=atof(optarg);
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

    case 'd':			//show defaults TD!
      break;
    case '?':
      rn=optopt;
    default:			//Ask for syntax help
      fprintf(stderr,
	      "Unknown or unsupported option of code %i(%c) passed\n"
	      "as argument\n"
	      ,rn,(char)rn);
      prochelp(EXIT_FAILURE,var_docs,&var_cfg);
      break;
    case 'h':
      prochelp(EXIT_SUCCESS,var_docs,&var_cfg);
      break;

    }
  }

  return 0;
}
