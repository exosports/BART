
/* \fcnfh
   Computes the extinction in cm-1 given the temperature and the wavenumber
   sampling

   @returns 0 on succes
*/
int
computecia(double *res,		/* Result is stored in this previously
				   allocated and zeroed array of size
				   wn->n */ 
	   float amagat,
	   double temp,
	   prop_samp *wn,
	   struct cia *ciast)
{
  struct ciapair *cia=ciast->p;
  int np=ciast->n,nt,nwn;
  PREC_CIA *e;
  long i,f,l;
  PREC_RES *wnv=wn->v;
  PREC_CIA *t,*a;

  //process each of the databases, skiping those which don't have the
  //right temperature range.
  while(np--){
    e=alloca(cia->nwn*sizeof(PREC_CIA));
    nt=cia->nt;
    t=cia->t;
    nwn=cia->nwn;
    if(temp<t[0]||temp>t[nt-1]||wnv[0]>cia->wn[nwn-1]||
       wnv[wn->n-1]<cia->wn[0])
      continue;

    //search for first element
    f=0;
    while(cia->wn[f++]<wnv[0]);
    f--;
    if(f) f--;

    for(i=f;i<cia->nwn;i++){
      a=cia->a[i];
#ifdef _USE_GSL
      gsl_interp_accel acc={0,0,0};
      gsl_interp *spl=gsl_interp_alloc(gsl_interp_cspline,nt);
      gsl_interp_init(spl,t,a,nt);
      e[i]=gsl_interp_eval(spl,t,a,temp,&acc);
      gsl_interp_free(spl);
#else
#error We cannot spline interpolate without GSL to obtain CIA opacities
#endif

      //skip from the last element;
      if(cia->wn[i]>wnv[nwn-1]){
	i++;
	break;
      }
    }
    l=i;

    //interpolate for the right wavenumbers skipping of limits values
#ifdef _USE_GSL
    gsl_interp_accel acc={0,0,0};
    gsl_interp *spl=gsl_interp_alloc(gsl_interp_cspline,l-f);
    gsl_interp_init(spl,cia->wn+f,e+f,l-f);
    for(i=0;i<wn->n;i++){
      if(wnv[i]<cia->wn[f]||wnv[i]>cia->wn[l])
	continue;
      res[i]+=gsl_interp_eval(spl,cia->wn+f,e+f,wnv[i],&acc)*amagat*amagat;
    }
    gsl_interp_free(spl);
#else
#error We cannot spline interpolate without GSL to obtain CIA opacities
#endif

    cia++;
  }

  return 0;
}


#if 0				/* Obsolete */

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

#endif /* obsolete */


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


  //\delfh
#if 0 				/* READLININFO.TXC */
/* tokenize a comma separated list of files and open them */
  //Find the maximum number of line info files, and check which of those
  //are available.
  char *file;
  int of=0,cf=Pnchr(th->f_line,',');
  transit->fp_line=(FILE **)calloc(cf+1,sizeof(FILE *));
  transit->f_line=(char *)calloc(strlen(th->f_line)+1,sizeof(char));
  file=strtok(th->f_line,",");
  th->f_line=NULL;
  //For each suggested file,
  while(file!=NULL){
    //if they can't be opened keep them in 'th'
    if((rn=fileexistopen(file,&(transit->fp_line[of])))!=1){
      transiterror(TERR_WARNING,
		   "Line info file '%s' is not available.\n"
		   " fileexistopen() error code %i.\n"
		   ,file,rn);
      strcat(th->f_line,file);
      strcat(th->f_line,",");
    }
    //otherwise move the name to 'transit'
    else{
      strcat(transit->f_line,file);
      strcat(transit->f_line,",");
      of++;
      th->na&=~TRH_FL;
    }
    file=strtok(NULL,",");
  }
  //If nothing was succesfully opened.
  if(!of)
    transiterror(TERR_SERIOUS,
		 "No line info files were available from:\n%s\n"
		 ,th->f_line);
  //Clean a bit, in 'transit' and 'th' if there was something left,
  //erase the final `,', else free the array and nullify the
  //pointer. Also fix the size of FILE pointers to what was used.
  th->f_line=(char *)realloc(th->f_line,rn=strlen(th->f_line));
  if(rn) th->f_line[rn-1]='\0';
  else th->f_line=NULL;
  transit->f_line=(char *)realloc(transit->f_line,rn=strlen(transit->f_line));
  if(rn) transit->f_line[rn-1]='\0';
  else transit->f_line=NULL;
  transit->fp_line=(FILE **)realloc(transit->fp_line,of*sizeof(FILE *));
  //Update flags to indicate non-accepted functions
  cf-=of;
  if(cf>TRH_FL_SOME) cf=TRH_FL_SOME;
  th->na|=cf;
#endif /* not used tokenization */

#if 0 				/* obsolte */
/*
  readdatastart: Initialize data file to be read from iniw from
  lineread. This function doesn't check for correct boundaries of
  data. You should run checkrange() before running this

  @returns 1  on success
           -1 on unexpected EOF
*/
int readdatastart(char *datafile,

		  double iniw,

		  PREC_LNDATA dbini,
		  float dwmark,
		  PREC_NREC *wmark,

		  FILE **fp)
{
  PREC_NREC i,alloc,mb;
  PREC_LNDATA wltmp;

  alloc=8;

  if((*fp=fopen(datafile,"r"))==NULL){
    transiterror(TERR_SERIOUS,
		 "Cannot open data file '%s', please check and try again.\n"
		 ,datafile);
  }    

  /* Finding starting point in datafile */

  i=(int)((iniw-dbini)/dwmark);
  datafileBS(*fp, wmark[i], wmark[i+1], iniw, &mb, sizeof(struct line_transition));
  transitDEBUG(20,verblevel,"Beginning found at position %li ",mb);
  if (mb){
    do{
      fseek(*fp,--mb*sizeof(struct line_transition),SEEK_SET);
      fread(&wltmp,sizeof(PREC_LNDATA),1,*fp);
    }while(wltmp>=iniw);
    mb++;
  }
  transitDEBUG(20,verblevel,"and then slide to %li\n",mb);
  fseek(*fp,mb*sizeof(struct line_transition),SEEK_SET);

  return 1;
}

/*
  readdatanext: Read next line information from datafile previously
  opened by readdatanext

  @returns 1 on success
           -1 on EOF
*/
inline int readdatanext(FILE *fp,
			struct line_transition *res,

			short *dbiso,    //Read from which database/isotopes?
			short ndbiso)      // that many
{
  int j;

  do{
    if(fread(res,sizeof(struct line_transition),1,fp)==0){
      transiterror(TERR_WARNING,
		   "End-of-file while reading single record. "
		   "Last wavelength read was %f\n\n"
		   ,(*res).wl);
      return -1;
    }
    transitprint(21,verblevel,"Wavelength:%f iso:%i\n",(*res).wl,
		 (*res).isoid);

    /* Following blocks checks for wanted isotopes */
    for(j=0;j<ndbiso;j++){
      if(dbiso[j]==(*res).isoid)
	break;
    }
    if(j>0&&j==ndbiso)    //Skipping this record because it is an
                          //unwanted isotope.
      continue;
  }while(0);

  return 1;
}
#endif /*obsolete*/
//\deluh

}

  //Integrate towards a constant value of extinction if we are in the
  //outmost layer
  else{
    if(ex[1]==ex[0])
      res= ex[0] * r0 * ( sqrt( rad[1] * rad[1] / r0 / r0 - 1) );
    else{
      PREC_RES alpha = ( ex[1] - ex[0] ) / dr;
      PREC_RES rm    = rad[1];
      if(alpha<0)
	res= - alpha * ( (rm + 2*ex[1]/alpha) * sqrt( rm * rm - r0 * r0)
			 - r0 * r0 * log( sqrt( rm * rm / r0 / r0 - 1) 
					  + rm / r0 )
			 ) / 2.0;
      else
	res=   alpha * ( (rm - 2*r0 - 2*ex[0]/alpha) * sqrt( rm*rm - r0*r0 )
			 - 2 * r0 * r0 * log( sqrt( rm * rm / r0 / r0 - 1)
					      + rm / r0 ) 
			 - r0*r0 * log( 2*r0 )
			 ) / 2.0;
    }
  }

