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
