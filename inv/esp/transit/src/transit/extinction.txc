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
      
  wnmar=tr->wnm;
  iniwn=tr->wns.i;
  wn=tr->wns.v;
  wni=wn[0]-wnmar;
  wnf=wn[nwn-1]+wnmar;
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
  csiso=  (PREC_CS *)  calloc(niso,sizeof(PREC_CS  ));
  ziso=   (PREC_ZREC *)calloc(niso,sizeof(PREC_ZREC));
  mass=   (PREC_ZREC *)calloc(niso,sizeof(PREC_ZREC));
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
    kiso=ex->k[r]=*ex->k+r*niso;
    alphal=ex->al[r]=*ex->al+r*niso;
    alphad=ex->ad[r]=*ex->ad+r*niso;

    //set some auxiliary variables.
    temp=tr->atm.t[r];

    //'propto\_adop' is proportional to the Doppler width, which in its
    //total splendor is
    //\[
    //\alpha_D = \frac{\Wn}{\sqrt{m}}\underbrace{
    //\frac{\sqrt{2k_B T \ln 2}}{c}}_{\mathrm{propto\_adop}}
    //\label{dopbr}
    //\]
    propto_adop=sqrt(2*KB*temp)*SQRTLN2/LS;

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
      for(j=0;j<neiso;j++)
	alphal[i]+=csiso[i]*tr->isov[j].d[r]/tr->isof[j].m
	  *sqrt(1/mass[i] + 1/tr->isof[j].m);
      alphal[i]*=propto_alor;

      //Following calculates doppler divided by central wavenumber.
      alphad[i]=propto_adop/sqrt(mass[i]);

      //Now get a new profile, 'profile[i]' has dimensions
      //'ex->vf'x'ex->lw[0]'. The latter is calculated by newprofile,
      //though.
      if((nwnh[i]=newprofile(profile[i],ex->vf,&ex->lw[0],dwn,
			     wn[0]*alphad[i],alphal[i],ex->ta)
	  )<1)
	transiterror(TERR_CRITICAL,
		     "newprofile() returned error code %i on its\n"
		     "first try for isotope %i.\n"
		     ,nwnh[i],i);
      //set 'w' to the last wavenumber index (because the lines are
      //sorted by wavelength.
      //And set 'cp[w]' to how many other wavenumbers bins to
      //recalculate the Voigt profile, it has to be at least 1.
      //verbline{voigtrec}
      w=nwn-1;
      cp=ex->recalc[i];
      cp[w]=(int)(ex->maxratio*wn[w]/dwn+0.5);
      if(!cp[w])
	cp[w]=1;

      //now set 'wrc[i]' to the next wavenumber that Voigt needs to be
      //recalculated, and put a 1 in that position in the 'cp[]' array.
      wrc[i]=w-cp[w];
      if(wrc[i]>=0)
	cp[wrc[i]]=1;

      //initialize 'wa[i]' to last wavenumber
      wa[i]=w;
    }

    //Compute the spectra!, proceed for every line.
    for(ln=0;ln<tr->n_l;ln++){
      /*
      if(ln!=10000&&ln!=10702&&ln!=10402)
	continue;
      if(ln<9000||ln>11000)
	continue;
      */

      wavn=WNU_O_WLU/line[ln].wl;
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
      transitDEBUG(20,verblevel,
		   "wavn:%g lgf:%g\n"
		   ,wavn,line[ln].lgf);
      //If it is beyond the last then just skip that line
      /* out of borders enabled =>change following */
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
      //if $'w'<='wrc'$ then recalcute Voigt
      if(w<=wrc[i]){
	//Find number of wavenumbers until the next recalculation, and
	//store as in \lin{voigtrec}
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

	//free old profile and recalculate voigt
	free(profile[i][0]);
	if((nwnh[i]=newprofile(profile[i],ex->vf,&ex->lw[w],dwn,
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

      //set upper and lower limits for Voigt spread
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

      //'wa[i]' is just the last wavelength per isotope.
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

