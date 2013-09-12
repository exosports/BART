/*
 * argum.c   - Process command line parameters for the transit program.
 *
 * Copyright (C) 2003-2006 Patricio Rojo (pato@astro.cornell.edu)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of version 2 of the GNU General 
 * Public License as published by the Free Software Foundation.
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


#include <time.h>
#include <transit.h>
#include <version_tr.h>
#include <math.h>
#include <pu/procopt.h>
#ifdef _USE_GSL
# include <gsl/gsl_version.h>
#endif /* _USE_GSL */
    

const static transit_ray_solution *raysols[] = {
  &slantpath,
  NULL
};

#ifndef EXTRACFGFILES
#  define PREPEXTRACFGFILES ""
#endif

#ifdef NODOTDEFAULT
#  define DOTCFGFILE ""
#  define DOTCFGFILENM "NO DEFAULT FILE"
#else
#  define DOTCFGFILE "./.transitrc"
#  define DOTCFGFILENM DOTCFGFILE
#  ifdef EXTRACFGFILES
#    define PREPEXTRACFGFILES ","EXTRACFGFILES
#  endif
#endif

/* \fcnfh
   process command line options, saving them in the hint structure.

   @returns 0 on success
 */
int processparameters(int argc, /* number of command line arguments */
		      char **argv, /* command line arguments */
		      struct transit *tr) /* structure to store
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
    CLA_GORBPAR,
    CLA_GORBPARFCT,
    CLA_GTIME,
    CLA_GTIMEFCT,
    CLA_GMASSRAD,
    CLA_GMASSRADFCT,
    CLA_OUTTAU,
    CLA_TOOMUCH,
    CLA_OUTTOOMUCH,
    CLA_RADFCT,
    CLA_WAVFCT,
    CLA_WNFCT,
    CLA_OUTSAMPLE,
    CLA_TAULEVEL,
    CLA_MODLEVEL,
    CLA_BLOWEX,
    CLA_TAUISO,
    CLA_MINELOW,
    CLA_CLOUDRAD,
    CLA_CLOUDFCT,
    CLA_CLOUDE,
    CLA_TRANSPARENT,
    CLA_DETEXT,
    CLA_DETCIA,
    CLA_DETTAU,
    CLA_CIAFILE,
    CLA_SAVEEXT,
    CLA_STARRAD,
  };

  //General help-option structure
  struct optdocs var_docs[]={
    {NULL,0,HELPTITLE,NULL,
     NULL,"GENERAL ARGUMENTS"},
    {"version",'V',no_argument,NULL,
     NULL,"Prints version number and exit"},
    {"help",'h',no_argument,NULL,
     NULL,"Prints list of possible parameters"},
    {"defaults",'d',no_argument,NULL,
     NULL,"Prints default values of the different variable"},
    {"verb",'v',no_argument,NULL,
     NULL,"Increase the verbose level by one\n"},
    {"quiet",'q',no_argument,NULL,
     NULL,"Decrease the verbose level to the minimum\n"},
    {"paramf",'p',ADDPARAMFILE,NULL,
     "filename","Use filename to read parameters in addition to\n"
     "default file(s): '" DOTCFGFILENM PREPEXTRACFGFILES"'"},

    {NULL,0,HELPTITLE,NULL,
     NULL,"INPUT/OUTPUT"},
    {"output",'o',required_argument,"-",
     "outfile","Change output file name, a dash (-)\n"
     "directs to standard output"},
    {"atm",CLA_ATMOSPHERE,required_argument,"-",
     "atmfile","File containing atmospheric info (Radius,\n"
     "pressure, temperature). A dash (-) indicates alternative\n"
     "input"},
    {"linedb",CLA_LINEDB,required_argument,"./res/lineread.tli",
     "linedb","File containing line information (TLI format,\n"
     "as given by 'lineread'"},
    {"outtoomuch",CLA_OUTTOOMUCH,required_argument,NULL,
     "filename","Ouputs depth where toomuch optical depth has been\n"
     "attained as a function of wavelength\n"},
    {"outsample",CLA_OUTSAMPLE,required_argument,NULL,
     "filename","Outputs sampling information. A dash (-) indicates\n"
     "standard input. By default there is no such output"},

    {NULL,0,HELPTITLE,NULL,
     NULL,"RADIUS OPTIONS (0-valued defaults would mean to use the"
     " values given by the atmosphere file)"},
    {"radius",'r',no_argument,NULL,
     NULL,"Interactively input radius parameters"},
    {"rad-low",CLA_RADLOW,required_argument,"0",
     "radius","Lower radius. 0 if you want to use atmospheric\n"
     "data minimum"},
    {"rad-high",CLA_RADHIGH,required_argument,"0",
     "radius","Higher radius. 0 if you want to use atmospheric\n"
     "data maximum"},
    {"rad-delt",CLA_RADDELT,required_argument,".5",
     "spacing","Radius spacing. 0 if you want to use atmospheric\n"
     "data spacing"},
    {"rad-fct",CLA_RADFCT,required_argument,"0",
     "factor","Radius factor. Multiplicating radius values by this\n"
     "gives centimeters.  If 0 then use atmosphere-file factor."},

    {NULL,0,HELPTITLE,NULL,
     NULL,"ATMPOSPHERE OPTIONS"},
    {"number-abund",CLA_NUMBERQ,no_argument,NULL,
     NULL,"Indicates that given abundances are by number rather\n"
     "than by mass"},
    {"onept",CLA_ONEPT,required_argument,NULL,
     "press,temp,extra_iso","Don't calculate transit spectra, just\n"
     "obtain spectra for a given pressure and temperature. Unless\n"
     "oneabund is also specified and has the correct number of\n"
     "isotopes, the abundances will be asked interactively"},
    {"oneextra",CLA_ONEEXTRA,required_argument,NULL,
     "mass1name1,mass2name2,...","It only has effect with --onept,\n"
     "a list of the atomic mass and names for the extra isotopes \n"
     "specified with --onept. If it doesn't have the right amount of\n"
     "values, the program will ask interactively"},
    {"oneabund",CLA_ONEABUND,required_argument,NULL,
     "q1,...","It also only has effect with --onept, a list of the\n"
     "abundances of the different isotopes. If it is omitted or\n"
     "doesn't have the right amount of values, the program will\n"
     "ask interactively. Note that the order of isotopes is the\n"
     "same given in the TLI data file"},
    {"onept-interactive",CLA_ONEINT,no_argument,NULL,
     NULL,"Wants to give abundances and pressure and temperature\n"
     "interactively through terminal input"},
    {"allowq",CLA_ALLOWQ,required_argument,"0.01",
     "value","How much less than one is accepted, so that no\n"
     "warning is issued if abundances don't ad up to that"},

    {NULL,0,HELPTITLE,NULL,
     NULL,"WAVELENGTH OPTIONS (all in fct units)"},
    {"wavelength",'w',no_argument,NULL,
     NULL,"Interactively input wavelength parameters"},
    {"wl-low",CLA_WAVLOW,required_argument,"0",
     "wavel","Lower wavelength. 0 if you want to use line\n"
     "data minimum"},
    {"wl-high",CLA_WAVHIGH,required_argument,"0",
     "wavel","Upper wavelength. 0 if you want to use line\n"
     "data maximum"},
    {"wl-delt",CLA_WAVDELT,required_argument,".00002",
     "spacing","Wavelength spacing. It cannot be 0 or less"},
    {"wl-osamp",CLA_WAVOSAMP,required_argument,"100",
     "integer","Wavelength oversampling. It cannot be 0 or less"},
    {"wl-fct",CLA_WAVFCT,required_argument,"0",
     "factor","Wavelength factor. Multiplicating wavelength values by\n"
     "this gives centimeters. If 0 or 1 then use centimeters"},
    {"wl-marg",CLA_WAVMARGIN,required_argument,"0.000001",
     "boundary","Not trustable range at boundary\n"
     "of line databases. Also transitions this much away from\n"
     "the requested range will be considered"},

    {NULL,0,HELPTITLE,NULL,
     NULL,"WAVENUMBER OPTIONS (all in cm-1)"},
    {"wavenumber",'n',no_argument,NULL,
     NULL,"Interactively input wavenumber parameters"},
    {"wn-low",CLA_WAVNLOW,required_argument,"0",
     "waven","Lower wavenumber. 0 if you want to use equivalent\n"
     "of the wavelength maximum"},
    {"wn-high",CLA_WAVNHIGH,required_argument,"0",
     "waven","Upper wavenumber. 0 if you want to use equivalent\n"
     "of the wavelength minimum"},
    {"wn-delt",CLA_WAVNDELT,required_argument,"0",
     "spacing","Wavenumber spacing. 0 if you want to have the\n"
     "same number of points as in the wavelength sampling"},
    {"wn-osamp",CLA_WAVNOSAMP,required_argument,"0",
     "integer","Wavenumber oversampling. 0 if you want the same\n"
     "value as for the wavelengths"},
    {"wn-fct",CLA_WNFCT,required_argument,"0",
     "factor","Output wavenumber factor. Multiplicating wavenumber\n"
     "values by this gives centimeters. If 0 then use wavelength's\n"
     "value. Note that this only applies to output, internally\n"
     "wavenumbers will always be in cm-1."},
    {"wn-marg",CLA_WAVNMARGIN,required_argument,"0",
     "boundary","Not trustable range in cm-1 at boundaries.\n"
     "Transitions this much away from the requested range will\n"
     "be considered. Use the maximum of the wavelength boundaries\n"
     "if this value is 0"},

    {NULL,0,HELPTITLE,NULL,
     NULL,"EXTINCTION CALCULATION OPTIONS:"},
    {"finebin",'f',required_argument,"5",
     "integer","Number of fine-bins to calculate the Voigt\n"
     "function"},
    {"nwidth",'a',required_argument,"50",
     "number","Number of the max-widths (the greater of Voigt\n"
     "or Doppler widths) that need to be contained in a calculated\n"
     "Voigt profile"},
    {"maxratio",'u',required_argument,"0.001",
     "uncert","Maximum allowed uncertainty in doppler width before\n"
     "recalculating profile"},
    {"per-iso",CLA_EXTPERISO,no_argument,NULL,
     NULL,"Calculates extinction per isotope, this allow displaying\n"
     "contribution from different isotopes, but also consumes more\n"
     "memory"},
    {"no-per-iso",CLA_NOEXTPERISO,no_argument,NULL,
     NULL,"Do not calculate extinction per isotope. Saves memory\n"
     "(this is the default)\n"},
    {"blowex",CLA_BLOWEX,required_argument,"1",
     "factor","Blow extinction by factor before computing tau. No\n"
     "physical significance of this variable, but only debugging"},
    {"minelow",CLA_MINELOW,required_argument,"0",
     "low-energy","Only use transitions with this minimum low energy\n"
     "(in cm-1)"},
    {"cloudrad",CLA_CLOUDRAD,required_argument,NULL,
     "radup,raddown","Make a cloud appear linearly from radup to raddown\n"
     "Units specified with '--cloudfct', or use radfct if there is none"},
    {"cloudfct",CLA_CLOUDFCT,required_argument,NULL,
     "factor","cloud radius values specified by '--cloudrad' will be\n"
     " multiplied by this to convert to cgs units\n"},
    {"cloudext",CLA_CLOUDE,required_argument,NULL,
     "extinction","Maximum extinction of the cloud, which opacity will\n"
     "linearly increase from 'radup' to 'raddown'\n"},
    {"detailext",CLA_DETEXT,required_argument,NULL,
     "filename:wn1,wn2,...","Save extinction at the particular\n"
     "wavenumbers in the specified filename"},
    {"detailcia",CLA_DETCIA,required_argument,NULL,
     "filename:wn1,wn2,...","Save extinction due to CIA at the particular\n"
     "wavenumbers in the specified filename"},
    {"cia", CLA_CIAFILE, required_argument, NULL,
     "filenames", "Use the indicated filenames for CIA opacities,\n"
     "it is a comma separated list"},
    {"saveext", CLA_SAVEEXT, required_argument, NULL,
     "filename", "Save extinction array in this file which won't need to\n"
     "be recomputed if only the radius scale (scale height) changes"},

    {NULL,0,HELPTITLE,NULL,
     NULL,"RESULTING RAY OPTIONS:"},
    {"solution",'s',required_argument,"Slant Path",
     "sol_name","Name of the kind of output solution ('slant path'\n"
     "is currently the only availabale alternative)"},
    {"toomuch",CLA_TOOMUCH,required_argument,"20",
     "optdepth","If optical depth for a particular path is larger\n"
     "than optdepth, then do not proceed to lower radius"},
    {"tauiso",CLA_TAUISO,required_argument,"0",
     "isoid","Compute tau only for isotope indexed in isoid (index which\n"
     "can actually be different from what you expect)"},
    {"outtau",CLA_OUTTAU,required_argument,"0",
     "#radius","Output is optical depth instead of modulation. It will be\n"
     "asked which radius to plot\n"},
    {"taulevel",CLA_TAULEVEL,required_argument,"1",
     "integer","Do a level integer integration for optical depth. 1\n"
     "is for constant index of refraction (better precision), use 2\n"
     "if it is variable."},
    {"modlevel",CLA_MODLEVEL,required_argument,"1",
     "integer","Do an integration of level <integer> to compute modulation.\n"
     "1 doesn't consider limb darkening. -1 doesn't consider limb darkening\n"
     "and additionally only returns the moduated radius at which extinction becomes\n"
     "one."},
    {"detailtau",CLA_DETTAU,required_argument,NULL,
     "filename:wn1,wn2,...","Save optical depth at the particular\n"
     "wavenumbers in the specified filename"},

    {NULL,0,HELPTITLE,NULL,
     NULL,"GEOMETRY PARAMETERS"},
    {"starrad",CLA_STARRAD, required_argument, "1.125",
     "radius_sun", "Stellar radius in solar radius"},
    {"g-orbpar",CLA_GORBPAR, required_argument, NULL,
     "smaxis,time,incl,ecc,long_node,arg_per","Orbital parameters, in"
     " the above order, to use the default of any of these (1,0,0,0,0,0),"
     " leave the corresponding field blank"},
    {"g-orbparfct",CLA_GORBPARFCT,required_argument,NULL,
     "unitsof:smaxis,time,incl,ecc,long_node,arg_per","Units of orbital"
     " parameters, in the above order, to use the default of any of these"
     " (AU,deg,hours,,deg,deg), leave the corresponding field blank"},
    {"transparent",CLA_TRANSPARENT,no_argument,NULL,
     NULL,"If selected the planet will have a maximum optical depth\n"
     " given by toomuch, it will never be totally opaque"},

    {NULL,0,0,NULL,NULL,NULL}
  };

  struct optcfg var_cfg;
  memset(&var_cfg,0,sizeof(var_cfg));
  var_cfg.contact="Patricio Rojo <pato@das.uchile.cl>";
  var_cfg.files=DOTCFGFILE PREPEXTRACFGFILES;
  var_cfg.columns=70;

  static struct transithint st_trh;
  memset(&st_trh, 0, sizeof(struct transithint));
  struct transithint *hints = tr->ds.th=&st_trh;

  st_trh.fl|=TRU_ATMASK1P|TRU_SAMPSPL|TRH_MASS;
  st_trh.verbnoise=4;
  st_trh.mass=1;
  /* TD: have this option user selectable */
  st_trh.tauiso=0;

  int rn,i;
  prop_samp *samp=NULL;
  char name[20],rc,*lp;
  char *sampv[]={"Initial","Final","Spacing","Oversampling integer for"};
  double rf;

  //preset values for detailed output
  struct detailfld *det=&hints->det.tau;
  strcpy(det->name,"optical depth");
  det=&hints->det.ext;
  strcpy(det->name,"extinction");
  det=&hints->det.cia;
  strcpy(det->name,"CIA extinction");
  

  procopt_debug=1;
  opterr=0;
  while(1){
    rn=procopt(argc,argv,var_docs,&var_cfg,NULL);
    if (rn==-1)
      break;

    transitDEBUG(21,verblevel,
		 "Processing option '%c', argum: %s\n"
		 ,rn,optarg);

    switch(rn){
    case CLA_CIAFILE:
      hints->ncia    = nchar(optarg,',') + 1;
      hints->ciafile = splitnzero_alloc(optarg,',');
      break;
    case CLA_SAVEEXT:
      hints->save.ext = strdup(optarg);
      break;
    case CLA_DETCIA:
      det=&hints->det.cia;
    case CLA_DETTAU:
      if(rn==CLA_DETTAU)
	det=&hints->det.tau;
    case CLA_DETEXT:
      if(rn==CLA_DETEXT)
	det=&hints->det.ext;

      //search for filename and save it
      if (det->n)
	free(det->ref);
      rn = sizeof(det->file);
      if(strlen(optarg)<rn) rn = strlen(optarg);
      i=0;
      while(i<rn)
	if(optarg[i++]==':') break;
      optarg[i-1]='\0';
      strcpy(det->file,optarg);
      //get wavenumbers to save
      det->n=getad(0,',',optarg+i,&det->ref);
      if(det->n<1||i==rn-1)
	transiterror(TERR_SERIOUS,
		     "Bad format for detailed %s parameter, no valid"
		     " wavenumbers\n"
		     ,det->name);
      break;
    case CLA_MINELOW:
      hints->minelow=atof(optarg);
      break;
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
    case CLA_OUTSAMPLE:
      if(hints->f_outsample) free_null(hints->f_outsample);
      hints->f_outsample=(char *)calloc(strlen(optarg)+1,sizeof(char));
      strcpy(hints->f_outsample,optarg);
      break;
    case CLA_OUTTOOMUCH:
      if(hints->f_toomuch) free_null(hints->f_toomuch);
      if(*optarg!='\0'){
	hints->f_toomuch=(char *)calloc(strlen(optarg)+1,sizeof(char));
	strcpy(hints->f_toomuch,optarg);
      }
      break;
    case CLA_OUTTAU:
      if(atoi(optarg))
	hints->fl|=TRU_OUTTAU;
      hints->ot=atoi(optarg)-1;
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
		       " temperature (%g), or number of extra isotopes (%g),\n"
		       " was not a correct floating point value. Actually,\n"
		       " the latter should be an integer\n"
		       ,hints->onept.p,hints->onept.t,rf);
	else
	  transiterror(TERR_SERIOUS,
		       "There was %i comma-separated fields instead of 3 for\n"
		       " '--onept' option"
		       ,-rn);
      }
      hints->onept.ne=(int)rf;
      if(rf!=hints->onept.ne)
	transiterror(TERR_SERIOUS,
		     "A non integer(%g) number of extra isotopes was given\n"
		     " with the option --onept\n"
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
		       "Bad format in the field #%i of --oneextra. It\n"
		       " doesn't have a valid value for mass. The field\n"
		       " should be <mass1><name1> with only an optional\n"
		       " dash betweeen the mass and name:\n %s\n"
		       ,i+1,optarg);
	if(*lp=='-') lp++;
	strncpy(hints->onept.n[i],lp,maxeisoname);
	optarg=lp;
	hints->onept.n[i][maxeisoname-1]='\0';
	while(*optarg++);
	if(lp==optarg)
	  transiterror(TERR_SERIOUS,
		       "Bad format in the field #%i of --oneextra. It\n"
		       " doesn't have a valid isotope name The field should\n"
		       " be <mass1><name1> with only an optional dash\n"
		       " betweeen the mass and name:\n %s\n"
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
    case CLA_RADFCT:
      hints->rads.fct=atof(optarg);
      break;

    case CLA_WAVLOW:
      hints->wavs.i=atof(optarg);
      break;
    case CLA_WAVHIGH:
      hints->wavs.f=atof(optarg);
      break;
    case CLA_WAVDELT:
      hints->wavs.d=atof(optarg);
      if(hints->wavs.d<=0)
	transiterror(TERR_SERIOUS,
		     "Wavelength spacing has to strictly positive,\n"
		     "greater than 0 instead of %g.\n"
		     ,hints->wavs.d);
      hints->wavs.n=0;
      hints->wavs.v=NULL;
      break;
    case CLA_WAVFCT:
      hints->wavs.fct=atof(optarg);
      break;
    case CLA_WAVOSAMP:
      hints->wavs.o=atof(optarg);
      break;
    case CLA_WAVMARGIN:
      hints->margin=atof(optarg);
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
    case CLA_WNFCT:
      hints->wns.fct=atof(optarg);
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

    case 'v':			//Increase verbose level
      verblevel++;
      break;

    case 'q':			//decrease verbose level
      verblevel=0;
      break;

    case 'V':			//Print version number and exit
      if(version_rc>0) snprintf(name,20,"-rc%i",version_rc);
      else name[0]='\0';
      printf("This is 'transit' version %i.%i%s\n\n",version,revision,name);
      exit(EXIT_SUCCESS);
      break;

    case 'd':			//show defaults TD!
      break;
    case '?':
      rn=optopt;
      transiterror(TERR_SERIOUS,
		   "Unknown, unsupported, or missing parameter to option "
		   "of code %i(%c) passed\n"
		   "as argument, use '-h' to see accepted options.\n"
		   ,rn,(char)rn);
      break;
    default:			//Ask for syntax help
      transiterror(TERR_CRITICAL,
		   "Even though option of code %i(%c) had a valid structure\n"
		   "element, it had no switch control statement. Code\n"
		   "need to be revised.\n"
		   ,rn,(char)rn);
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
    case CLA_BLOWEX:
      hints->blowex = atof(optarg);
      break;

    case CLA_STARRAD:
      hints->sg.starrad = atof(optarg);
      break;
    case CLA_GORBPAR:
      getnd(6, ',', optarg,
	    &hints->sg.smaxis,&hints->sg.time,&hints->sg.incl,
	    &hints->sg.ecc,&hints->sg.lnode,&hints->sg.aper);
      break;
    case CLA_GORBPARFCT:
      getnd(6,',',optarg,
	    &hints->sg.smaxisfct,&hints->sg.timefct,&hints->sg.inclfct,
	    &hints->sg.eccfct,&hints->sg.lnodefct,&hints->sg.aperfct);
      break;
    case CLA_TRANSPARENT:
      hints->sg.transpplanet=1;
      break;

    case CLA_TOOMUCH:
      hints->toomuch=atof(optarg);
      break;
    case CLA_TAUISO:
      hints->tauiso=atoi(optarg);
      break;
    case CLA_TAULEVEL:
      hints->taulevel=atoi(optarg);
      break;
    case CLA_MODLEVEL:
      hints->modlevel=atoi(optarg);
      break;
    case CLA_CLOUDRAD:
      hints->cl.rini=strtod(optarg,&optarg);
      if(*optarg!=','||optarg[1]=='\0')
	transiterror(TERR_SERIOUS,
		     "Syntax error in option '--cloudrad', parameters need\n"
		     " to be radup,raddown");
      hints->cl.rfin=strtod(optarg+1,NULL);
      if(hints->cl.rini<hints->cl.rfin ||
	 (hints->cl.rfin<=0&&hints->cl.rini!=0))
	transiterror(TERR_SERIOUS,
		     "Syntax error in option '--cloudrad', radup(%g) needs\n"
		     " to be bigger than raddown (%g) and both greater than\n"
		     " zero\n"
		     ,hints->cl.rini,hints->cl.rfin);
      break;
    case CLA_CLOUDFCT:
      hints->cl.rfct=atof(optarg);
      break;
    case CLA_CLOUDE:
      hints->cl.maxe=atof(optarg);
      break;
    }

    
  }

  procopt_free();

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
    sol++;
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


  //Accept output file
  if(th->f_out)
    tr->f_out=th->f_out;
  else{
    tr->f_out=(char *)calloc(2,sizeof(char));
    tr->f_out[0]='-';
  }

  //Accept toomuch output file
  tr->f_toomuch=th->f_toomuch;
  tr->f_outsample=th->f_outsample;

  //Accept hinted solution type if it exists
  if(acceptsoltype(&tr->sol,th->solname)!=0){
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

  //set hinted geometry hints
  setgeomhint(tr);

  //accept hints for detailed output
  static struct detailout det;
  memcpy(&det,&th->det,sizeof(struct detailout));
  tr->ds.det=&det;

  //return success
  return 0;
}


/*\fcnfh
  Saves hints structure
*/
void
savehint(FILE *out,
	 struct transithint *hints)
{
  //save main structure
  fwrite(hints,sizeof(struct transithint),1,out);

  //save strings
  savestr(out,hints->f_atm);
  savestr(out,hints->f_line);
  savestr(out,hints->f_out);
  savestr(out,hints->f_toomuch);
  savestr(out,hints->f_outsample);
  savestr(out,hints->solname);
  for(int i=0 ; i<hints->ncia ; i++){
    savestr(out, hints->ciafile[i]);
  }

  //save sub-structures
  savesample_arr(out,&hints->rads);
  savesample_arr(out,&hints->wavs);
  savesample_arr(out,&hints->wns);
  savesample_arr(out,&hints->ips);

  saveonept_arr(out,&hints->onept);
}

/* \fcnfh
   Restore hints structure, the structure needs to have been allocated
   before

   @returns 0 on success
            -1 if not all the expected information is read
	    -2 if info read is wrong
	    -3 if cannot allocate memory
	    1 if information read was suspicious
*/
int 
resthint(FILE *in,
	 struct transithint *hint)
{
  int rn,res=0;

  //Restore main structure
  rn=fread(hint,sizeof(struct transithint),1,in);
  if(rn<0) return rn; else res+=rn;

  //restore strings
  rn=reststr(in,&hint->f_atm);
  if(rn<0) return rn; else res+=rn;
  rn=reststr(in,&hint->f_line);
  if(rn<0) return rn; else res+=rn;
  rn=reststr(in,&hint->f_out);
  if(rn<0) return rn; else res+=rn;
  rn=reststr(in,&hint->f_toomuch);
  if(rn<0) return rn; else res+=rn;
  rn=reststr(in,&hint->f_outsample);
  if(rn<0) return rn; else res+=rn;
  rn=reststr(in,&hint->solname);
  if(rn<0) return rn; else res+=rn;
  for(int i=0 ; i<hint->ncia ; i++){
    rn=reststr(in,hint->ciafile+i);
    if(rn<0) return rn; else res+=rn;
  }

  //restore sub-structures
  restsample_arr(in,&hint->rads);
  if(rn<0) return rn; else res+=rn;
  restsample_arr(in,&hint->wavs);
  if(rn<0) return rn; else res+=rn;
  restsample_arr(in,&hint->wns);
  if(rn<0) return rn; else res+=rn;
  restsample_arr(in,&hint->ips);
  if(rn<0) return rn; else res+=rn;

  restonept_arr(in,&hint->onept);
  if(rn<0) return rn; else res+=rn;

  return res;
}


void
printintro()
{
  char rcname[20];
  if(version_rc>0) snprintf(rcname,20,"-rc%i",version_rc);
  else rcname[0]='\0';
  transitprint(1,verblevel,
	       "-----------------------------------------------\n"
	       "                TRANSIT v%i.%i%s\n"
	       "-----------------------------------------------\n"
	       ,version,revision,rcname);
  time_t tim = time(NULL);
  transitprint(2,verblevel,
	       "Started on %s\n", ctime(&tim));
}


/* \fcnfh
   Frees hints structure
*/
void
freemem_hints(struct transithint *h)
{
  //free strings which are copied into transit
  free(h->f_atm);
  free(h->f_line);
  free(h->f_out);
  free(h->f_toomuch);
  free(h->f_outsample);

  //free other strings
  free(h->solname);
  if (h->ncia){
    free(h->ciafile[0]);
    free(h->ciafile);
  }

  //free sub-structures
  freemem_onept(&h->onept);
  freemem_samp(&h->rads);
  freemem_samp(&h->wavs);
  freemem_samp(&h->wns);
  freemem_samp(&h->ips);
  /* TD: Free sve if it is ever enabled 
     freesaves(&h->save); */

  freemem_cloud(&h->cl);
  freemem_detailout(&h->det);


}

void
freemem_cloud(struct extcloud *c)
{
}

void
freemem_detailout(struct detailout *d)
{
  freemem_detailfld(&d->ext);
  freemem_detailfld(&d->tau);
  freemem_detailfld(&d->cia);
}

void
freemem_detailfld(struct detailfld *f)
{
  if (f->n)
    free(f->ref);
}
