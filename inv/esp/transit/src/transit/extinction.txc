/*
 * extinction.c 
 * extinction.txc - Computes extinction coefficient. Component of the
 *                  Transit program.
 *
 * Copyright (C) 2004 Patricio Rojo (pato@astro.cornell.edu)
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

/*\fcnfh
  This funciton fills up the extinction information in tr->ds.ex. It
  uses 
  extwn: Scattering parameters should be added at some point here.

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
  struct line_transition *lt=&(tr->ds.li->lt);
  PREC_LNDATA *ltgf=lt->gf;
  PREC_LNDATA *ltelow=lt->elow;
  PREC_LNDATA *ltwl=lt->wl;
  short *ltisoid=lt->isoid;
  double efct=lt->efct;
  double wfct=lt->wfct;
  PREC_RES *k,**kiso,*wn,dwn,wavn,iniwn,wni,wnf;
  PREC_NSAMP nrad,nwn;
  int neiso,niso,nisoalloc;
  int r,i,ln;
  int w,*wa,subw;
  int j,maxj,minj,*nwnh;
  PREC_VOIGT ***profile, *profwn;
  PREC_VOIGTP *alphal,*alphad;
  int *wrc;
  double propto_adop,propto_alor, propto_k;
  PREC_ATM temp, *densiso;
  PREC_CS *csiso;
  PREC_ZREC *ziso;
  PREC_ZREC *mass;
  _Bool extinctperiso;
  struct isotopes *iso=tr->ds.iso;
  long nlines=tr->ds.li->n_l;

  transitcheckcalled(tr->pi,"extwn",4,
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
  nwn=tr->wns.n;
  wni=wn[0]-tr->wnmi;
  wnf=wn[nwn-1]+tr->wnmf;
  dwn=tr->wns.d/tr->wns.o;
  nrad=rad->n;
  neiso=iso->n_e;
  nisoalloc=niso=iso->n_i;
  /*TD: enable nisoalloc, currently alloc is used */
  //Do not allocate memory in multidimensional arrays if we are ignoring
  //that particular isotope. We are not using this in unidimensional
  //arrays because of the extra hassle.
  for(i=0;i<niso;i++)
    if(iso->isodo[i]==ignore)
      nisoalloc--;
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

  //following to store isotope dens, these are auxiliary variables.
  densiso=(PREC_ATM *) alloca(niso*sizeof(PREC_ATM ));
  csiso=  (PREC_CS *)  alloca(niso*sizeof(PREC_CS  ));
  ziso=   (PREC_ZREC *)alloca(niso*sizeof(PREC_ZREC));
  mass=   (PREC_ZREC *)alloca(niso*sizeof(PREC_ZREC));
  nwnh=   (int *)      alloca(niso*sizeof(int));

  //allocate array for the voigt profile, also auxiliary.
  wa=(int *) alloca(niso*sizeof(int));
  wrc=(int *)alloca(niso*sizeof(int));
  profile=(PREC_VOIGT ***)calloc(niso,sizeof(PREC_VOIGT **));
  *profile=(PREC_VOIGT **)calloc(niso*ex->vf,sizeof(PREC_VOIGT *));
  for(i=1;i<niso;i++)
    profile[i]=*profile+i*ex->vf;

  extinctperiso=(tr->fl&TRU_EXTINPERISO);
  //allocate array for extinctions and widths, initialization of second
  //index will be done in the loop (\lin{kini})
  ex->e=(PREC_RES ***)calloc(nrad,sizeof(PREC_RES **));
  *ex->e=(PREC_RES **)calloc(niso*nrad,sizeof(PREC_RES *));
  i=extinctperiso?niso:1;
  if((ex->e[0][0]=(PREC_RES *)calloc(nrad*i*nwn,sizeof(PREC_RES)))==NULL)
    transiterror(TERR_CRITICAL|TERR_ALLOC,
		 "Unable to allocate %li = %li*%li*%li to calculate\n"
		 "extinction for every radii, %stry to shorten the wavenumber\n"
		 "range\n"
		 ,nrad*i*nwn,nrad,i,nwn,extinctperiso?"try disabling exctinction per\n"
		 "isotope (option --no-per-iso), or ":"");

  alphal=(PREC_VOIGTP *)calloc(niso,sizeof(PREC_VOIGTP));
  alphad=(PREC_VOIGTP *)calloc(niso,sizeof(PREC_VOIGTP));

  //For each radius (index 'r')
  for(r=0;r<nrad;r++){

    transitprint(2,verblevel,"Radius %i: %g[cm]\n",r,rad->fct*rad->v[r]);

    //Initialization of 2nd dimension of extinction array.
    //\linelabel{kini}
    kiso=ex->e[r]=*ex->e+r*niso;
    i=extinctperiso?niso:1;
    ex->e[r][0]=ex->e[0][0]+r*i*nwn;

    //set some auxiliary variables.
    temp=tr->atm.t[r];

    //'propto\_adop' is proportional to the Doppler width, which in its
    //total splendor is
    //\[
    //\alpha_D = \frac{\Wn}{\sqrt{m}}\underbrace{
    //\frac{\sqrt{2k_B T \ln 2}}{c}}_{\mathrm{propto\_adop}}
    //\label{dopbr}
    //\]
    propto_adop=sqrt(2*KB*temp/AMU)*SQRTLN2/LS;

    //'propto\_alor' is proportional to the Lorenz width, which in its
    //total splendor is
    //\[
    //\alpha_L = \underbrace{\Gamma_{nat}}_{\mathrm{ignored}}
    //+\underbrace{\frac{1}{\pi c}\sqrt{\frac{2\kb T}{\pi}}}
    //_{\mathrm{propto\_alor}}
    //\sigma_c\sum_{\mathrm{collisioners}}n_i \sqrt{\left(\frac{1}{m_r}
    //+\frac{1}{m_i}\right)}.
    //\label{lorwidth}
    //\]
    propto_alor=sqrt(temp*2*KB/PI/AMU)/AMU/LS/PI;

    //Initialize a voigt profile for every isotope as well for the
    //mass, ziso, densiso and csiso arrays
    for(i=0;i<niso;i++){
      //If this isotope is marked as ignore (no density info) continue
      //with the next one.
      if(iso->isodo[i]==ignore)
	continue;

      kiso[i]=kiso[0];
      if(extinctperiso)
	kiso[i]+=nwn*i;
      mass[i]=iso->isof[i].m;
      ziso[i]=iso->isov[i].z[r];
      densiso[i]=iso->isov[i].d[r];
      csiso[i]=iso->isov[i].c[r];

      //Calculate lorentz with every isotope except those ignored
      alphal[i]=0;
      for(j=0;j<neiso;j++)
	if(iso->isodo[j]!=ignore)
	  alphal[i]+=iso->isov[j].d[r]/iso->isof[j].m
	    *sqrt(1/mass[i] + 1/iso->isof[j].m);

      alphal[i]*=csiso[i]*propto_alor;

      //Following calculates doppler divided by central wavenumber.
      alphad[i]=propto_adop/sqrt(mass[i]);

      //Now get a new profile, 'profile[i]' has dimensions
      //'ex->vf'x(2*'nwnh[i]'+1). The latter is calculated by newprofile,
      //though.
      if((nwnh[i]=newprofile(profile[i],ex->vf,dwn,
			     wn[0]*alphad[i],alphal[i],ex->ta)
	  )<1)
	transiterror(TERR_CRITICAL,
		     "newprofile() returned error code %i on its\n"
		     "first try for isotope %i.\n"
		     ,nwnh[i],i);
      //set 'w' to the last wavenumber index (because the lines are
      //sorted by wavelength.
      //And set 'j' to how many other wavenumbers bins to
      //recalculate the Voigt profile, it has to be at least 1. Then set
      //'wrc[i]' to the wavelength it has to be recalculated.
      //verbline{voigtrec}
      w=nwn-1;
      j=(int)(ex->maxratio*wn[w]/dwn+0.5);
      if(!j) j=1;
      wrc[i]=w-j;

      //initialize 'wa[i]' to last wavenumber
      wa[i]=w;
    }

    //Compute the spectra!, proceed for every line.
    for(ln=0;ln<nlines;ln++){
      /*
      if(ln!=10000&&ln!=10702&&ln!=10402)
	continue;
      if(ln<9000||ln>11000)
	continue;
      */

      wavn=1.0/wfct/ltwl[ln];
      /* 
	 when out of borders enabled
	 if(wavn<wni||wavn>wnf)
	 continue;
      */
      //if it is beyond lower limit
      if(wavn<iniwn)
	continue;
      /* 
	 This is when out of borders is enabled
	 w=-(long)((iniwn-wavn)/dwn+1)
      */
      else
	w=(wavn-iniwn)/dwn;
      transitDEBUG(21,verblevel,
		   "wavn:%g lgf:%g\n"
		   ,wavn,ltgf[ln]);
      //If it is beyond the last then just skip that line
      /* out of borders enabled =>change following */
      if(w>=nwn)
	continue;

      subw=ex->vf*(wavn-w*dwn-iniwn)/dwn;
      i=ltisoid[ln];
      k=kiso[i];

      //If this isotope is marked as ignore (no density info) continue
      //with the next transition.
      if(iso->isodo[i]==ignore)
	continue;

      transitASSERT(wa[i]!=-1&&wa[i]<w,
		    "Database is not ordered!, previous wavenumber was\n"
		    "at index %i, new one at %i (it should have been smaller)\n"
		    ,wa,w);
      //if $'w'<='wrc'$ then recalcute Voigt
      if(w<=wrc[i]){
	//Find number of wavenumbers until the next recalculation, and
	//store as in \lin{voigtrec}
	j=(int)(ex->maxratio*wn[w]/dwn+0.5);
	if(!j) j=1;
	wrc[i]=w-j;
	transitDEBUG(22,verblevel,
		     "Recalculating voigt for isotope %i... current\n"
		     "wavenumber %i, next wavenumber %i/%i\n"
		     ,i,w,wrc[i],nwn);

	//free old profile and recalculate voigt
	free(profile[i][0]);
	if((nwnh[i]=newprofile(profile[i],ex->vf,dwn,
			       wn[w]*alphad[i],alphal[i],ex->ta)
	    )<1)
	  transiterror(TERR_CRITICAL,
		       "newprofile() returned error code %i for\n"
		       "isotope %i\n"
		       ,nwnh[i],i);
      }

      //Calculate opacity coefficient less the voigt spread
      /* CAVEATS: _mass_ densitty 
                  _log_ gf */
      propto_k=densiso[i]	//mass density
	*SIGCTE			//Constant in sigma
	*ltgf[ln]		//gf
	*exp(-EXPCTE*efct*ltelow[ln]/temp) //Level population
	*(1-exp(-EXPCTE*wavn/temp)) //induced emission
	/mass[i]		//mass
	/ziso[i];		//Partition function

      transitDEBUG(21,verblevel,
		   "i=%i   temp=%g   Elow=%g\n"
		   "aD=%.7g   aL=%.7g\n"
		   "wl=%.10g  wn=%.10g\n"
		   "k= %12.5g  //densiso[i] \n"
		   "  *%12.5g  //SIGCTE\n"
		   "  *%12.5g  //ltgf[ln]\n"
		   "  *%12.5g  //exp(-EXPCTE*ltelow[ln]/temp)\n"
		   "  *%12.5g  //(1-exp(-EXPCTE*wavn/temp))\n"
		   "  /%12.5g  //mass[i]\n"
		   "  /%12.5g  //ziso[i]\n"
		   " = %12.5g   //extinction\n\n"
		   ,i,temp,ltelow[ln]
		   ,alphad[i]*wavn,alphal[i]
		   ,ltwl[ln],1.0/wfct/ltwl[ln]/tr->wavs.fct
		   ,densiso[i]
		   ,SIGCTE
		   ,ltgf[ln]
		   ,exp(-EXPCTE*ltelow[ln]/temp)
		   ,(1-exp(-EXPCTE*wavn/temp))
		   ,mass[i]
		   ,ziso[i]
		   ,propto_k);

      //set 'profwn' such that the index mimic wavenumber's array
      profwn=profile[i][subw]+nwnh[i]-w;

      //set upper and lower limits for Voigt spread
      minj=w-nwnh[i];
      if(minj<0)
	minj=0;
      maxj=w+nwnh[i];
      if(maxj>nwn)
	maxj=nwn;

      //distribute the oscillator strength according to the voigt
      //profile
      for(j=minj;j<maxj;j++)
	k[j]+=propto_k
	  *profwn[j];

      //'wa[i]' is just the last wavelength per isotope.
      wa[i]=w;
    }
    //Free the profiles of every non-ignored isotopes
    for(i=0;i<niso;i++)
      if(iso->isodo[i]!=ignore)
	free(profile[i][0]);
  }


  //free memory that is no longer needed.
  freemem_lineinfotrans(tr->ds.li,&tr->pi);
  freemem_isotopes(tr->ds.iso,&tr->pi);

   //save current status if requested.
  savefile_extwn(tr);

  //Set porogress indicator, and print and output extinction if one P,T
  //was desired, otherwise return success
  tr->pi|=TRPI_EXTWN;
  if(tr->rads.n==1)
    printone(tr);
  return 0;
}


/* \fcnfh
   calculates a new voigt profile

   @returns number of points to center wavelength
*/
inline int newprofile(PREC_VOIGT **pr, /* output 2d profile */
		      int vf, /* 1st already allocated dimension */
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
  nvgt=2*wvgt/dwn+1;

  //Basic check that 'lor' or 'dop' are within sense
  if(nvgt<0)
    transiterror(TERR_CRITICAL,
		 "Number of voigt bins are not positive!. Possibly out of\n"
		 "a nonsensical doppler(%g) or Lorenz(%g) profile.\n"
		 ,dop,lor);

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



/* \fcnfh
   Printout for one P,T conditions
*/
void
printone(struct transit *tr)
{
  int rn;
  FILE *out=stdout;
  struct isotopes *iso=tr->ds.iso;

  //open file
  if(tr->f_out&&tr->f_out[0]!='-')
    out=fopen(tr->f_out,"w");

  transitprint(1,verblevel,
	       "\nPrinting extinction for one radius (at %gcm) in '%s'\n"
	       ,tr->rads.v[0],tr->f_out?tr->f_out:"standard output");

  //print!
  fprintf(out,
	  "#wavenumber[cm-1]\twavelength[nm]\textinction[cm-1]\t"
	  "cross-section[cm2]\n");
  for(rn=0;rn<tr->wns.n;rn++)
    fprintf(out,"%12.6f%14.6f%17.7g%17.7g\n"
	    ,tr->wns.fct*tr->wns.v[rn],WNU_O_WLU/tr->wns.v[rn]/tr->wns.fct,
	    tr->ds.ex->e[0][0][rn],
	    AMU*tr->ds.ex->e[0][0][rn]*iso->isof[0].m/iso->isov[0].d[0]);

  exit(EXIT_SUCCESS);
}


/* \fcnfh
   Frees structure allocated by extwn

   @returns 0 on success
*/
int
freemem_extinction(struct extinction *ex, /* Extinciton info */
		   long *pi)	/* progress indicator flags from which
				   to clear */
{
  //free arrays
  free(ex->e[0][0]);
  free(ex->e[0]);
  free(ex->e);

  //clear indicator and return success
  *pi&=!(TRPI_EXTWN);
  return 0;
}


/* \fcnfh 
   Saves all the memory that is going to be used after running extwn()

   @returns 0 on success
*/
int
savefile_extwn(struct transit *tr)
{


  return 0;
}
