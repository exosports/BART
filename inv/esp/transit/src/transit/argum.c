/*
 * argum.c
 * argum.txc - Process command line parameters for the transit program.
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


#include <transit.h>
#include <transit_version.h>
#include <math.h>
#include <util/procopt.h>

const static transit_ray_solution *raysols[] = {
  &slantpath,
  NULL
};


#ifndef EXTRACFGFILES
#define PREPEXTRACFGFILES ""
#else
#define PREPEXTRACFGFILES ","EXTRACFGFILES
#endif



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
    CLA_NUMBERQ,
    CLA_ALLOWQ,
    CLA_EXTPERISO,
    CLA_NOEXTPERISO,
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
    {"paramf",ADDPARAMFILE,'p',
     "filename","Use filename to read parameters in addition to\n"
     "default file(s): './.transitrc"PREPEXTRACFGFILES"'"},

    {NULL,HELPTITLE,0,
     NULL,"INPUT/OUTPUT"},
    {"output",required_argument,'o',
     "outfile","Change output file name, a dash (-)\n"
     "directs to standard output"},
    {"atm",required_argument,CLA_ATMOSPHERE,
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


    {NULL,HELPTITLE,0,
     NULL,"ATMPOSPHERE OPTIONS"},
    {"number-abund",no_argument,CLA_NUMBERQ,
     NULL,"Indicates that given abundances are by number rather\n"
     "than by mass"},
    {"oneptn",required_argument,CLA_ONEPT,
     "press,temp,extra_iso","Don't calculate transit spectra, just\n"
     "obtain spectra for a given pressure and temperature. Unless\n"
     "oneabund is also specified and has the correct number of\n"
     "isotopes, the abundances will be asked interactively"},
    {"oneextra",required_argument,CLA_ONEEXTRA,
     "mass1name1,mass2name2,...","It only has effect with --onept,\n"
     "a list of the atomic mass and names for the hitherto specified\n"
     "extra isotopes. If it doesn't have the right amount of values,\n"
     "the program will ask interactively"},
    {"oneabund",required_argument,CLA_ONEABUND,
     "q1,...","It also only has effect with --onept, a list of the\n"
     "abundances of the different isotopes. If it is omitted or\n"
     "doesn't have the right amount of values, the program will\n"
     "ask interactively. Note that the order of isotopes is the\n"
     "same given in the TWII data file"},
    {"onept-interactive",no_argument,CLA_ONEINT,
     NULL,"Wants to give abundances and pressure and temperature\n"
     "interactively through terminal input"},
    {"allowq",required_argument,CLA_ALLOWQ,
     "value","How much less than one is accepted, so that no\n"
     "warning is issued if abundances don't ad up to that"},

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
     NULL,"EXTINCTION CALCULATION OPTIONS:"},
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
    {"per-iso",no_argument,CLA_EXTPERISO,
     NULL,"Calculates extinction per isotope, this allow displaying\n"
     "contribution from different isotopes, but also consumes more\n"
     "memory"},
    {"no-per-iso",no_argument,CLA_NOEXTPERISO,
     NULL,"Do not calculate extinction per isotope. Saves memory\n"
     "(this is the default)\n"},

    {NULL,HELPTITLE,0,
     NULL,"RESULTING RAY OPTIONS:"},
    {"solution",required_argument,'s',
     "sol_name","Name of the kind of output solution ('slant path'\n"
     "is currently the only availabale alternative)"},

    {NULL,HELPTITLE,0,
     NULL,"OBSERVATIONAL OPTIONS:"},
    {"telres",required_argument,'t',
     "width","Gaussian width of telescope resolution in nm"},

    {NULL,0,0,NULL,NULL}
  };

  struct optcfg var_cfg;
  memset(&var_cfg,0,sizeof(var_cfg));
  var_cfg.contact="Patricio Rojo <pato@astro.cornell.edu>";
  var_cfg.files="./.transitrc"PREPEXTRACFGFILES;

  int rn,i;
  prop_samp *samp;
  char name[20],rc,*lp;
  char *sampv[]={"Initial","Final","Spacing","Oversampling integer for"};
  double rf;

  procopt_debug=1;
  opterr=0;
  while(1){
    rn=getprocopt(argc,argv,var_docs,&var_cfg,NULL);
    if (rn==-1)
      break;

    transitDEBUG(21,verblevel,
		 "Processing option '%c', argum: %s\n"
		 ,rn,optarg);

    switch(rn){
    case 's':
      hints->solname=(char *)realloc(hints->solname,strlen(optarg)+1);
      strcpy(hints->solname,optarg);
      break;
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
	  if(!rc)
	    break;
	  switch(rc){
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
    case CLA_ALLOWQ:
      hints->allowrq=atoi(optarg);
      break;
    case CLA_NUMBERQ:
      hints->mass=0;
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
      break;
    case CLA_ONEABUND:
      if((hints->onept.nq=getad(0,',',optarg,&hints->onept.q))<1)
	transiterror(TERR_SERIOUS,
		     "None of the given isotope abundances were accepted\n"
		     "%s\n"
		     ,optarg);
      transitprint(2,verblevel,
		   "%i abundance isotopes were correctly given: %s\n"
		   ,hints->onept.nq,optarg);
      break;

    case CLA_ONEEXTRA:
      rn=ncharchg(optarg,',','\0')+1;
      if(hints->onept.n){
	free(hints->onept.n);
	free(hints->onept.n[0]);
	free(hints->onept.m);
      }
      hints->onept.n=(char **)calloc(rn,sizeof(char *));
      hints->onept.n[0]=(char *)calloc(rn*maxeisoname,sizeof(char));
      hints->onept.m=(PREC_ZREC *)calloc(rn,sizeof(PREC_ZREC));

      for(i=0;i<rn;i++){
	hints->onept.n[i]=hints->onept.n[0]+i*maxeisoname;
	hints->onept.m[i]=strtod(optarg,&lp);
	if(lp==optarg)
	  transiterror(TERR_SERIOUS,
		       "Bad format in the field #%i of --oneextra. It doesn't have\n"
		       " a valid value for mass. The field should be <mass1><name1>\n"
		       " with only an optional dash betweeen the mass and name:\n %s\n"
		       ,i+1,optarg);
	if(*lp=='-') lp++;
	strncpy(hints->onept.n[i],lp,maxeisoname);
	optarg=lp;
	hints->onept.n[i][maxeisoname-1]='\0';
	while(*optarg++);
	if(lp==optarg)
	  transiterror(TERR_SERIOUS,
		       "Bad format in the field #%i of --oneextra. It doesn't have\n"
		       " a valid isotope name The field should be <mass1><name1>\n"
		       " with only an optional dash betweeen the mass and name:\n %s\n"
		       ,i+1,optarg);
      }
      hints->onept.nm=rn;
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
      hints->wns.i=atof(optarg);
      break;
    case CLA_WAVNHIGH:
      hints->wns.f=atof(optarg);
      break;
    case CLA_WAVNDELT:
      hints->wns.d=atof(optarg);
      hints->wns.n=0;
      hints->wns.v=NULL;
      break;
    case CLA_WAVNOSAMP:
      hints->wns.o=atof(optarg);
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
      if(revision<0) snprintf(name,20,"pre%i",-revision);
      else snprintf(name,20,".%i",revision);
      printf("This is 'transit' version %i%s\n\n",version,name);
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
      prochelp(EXIT_FAILURE);
      break;
    case 'h':
      prochelp(EXIT_SUCCESS);
      break;
    case CLA_EXTPERISO:
      hints->fl|=TRU_EXTINPERISO;
      break;
    case CLA_NOEXTPERISO:
      hints->fl&=~TRU_EXTINPERISO;
      break;
    }
  }


  return 0;
}


/* \fcnfh
   Look for a match in the solution name

   @returns 0 on success
            -1 if no matching name was available
 */
int acceptsoltype(transit_ray_solution **sol,
		  char *hname)
{
  *sol=(transit_ray_solution *)raysols[0];
  int len;

  len=strlen(hname);

  while(*sol){
    if(strncasecmp(hname,(*sol)->name,len)==0)
      return 0;
    *sol++;
  }

  return -1;
}


/* \fcnfh
   All the hints that need to be accepted at the beginning are set
   here.

   @returns 0 on success
 */
int
acceptgenhints(struct transit *tr) /* transit structure */
{
  struct transithint *th=tr->ds.th;

  if(th->na&TRH_FO)
    transitaccepthint(tr->f_out,th->f_out,th->fl,TRH_FO);
  //Accept hinted solution type if it exists
  if(!(th->na&TRH_ST)||acceptsoltype(&tr->sol,th->solname)!=0){
    transit_ray_solution **sol=(transit_ray_solution **)raysols;
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Solution kind '%s' is invalid!\n"
		 "Currently Accepted are:\n"
		 ,th->solname);
    while(*sol)
      transiterror(TERR_SERIOUS|TERR_NOPREAMBLE|TERR_ALLOWCONT,
		   " %s\n",(*sol++)->name);
    exit(EXIT_FAILURE);
  }

  return 0;
}
