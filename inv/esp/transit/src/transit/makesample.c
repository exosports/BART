/*
 * makesample.c - create array sampling after initial, final and spacing
 *                parameters. Component of the Transit program.
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
  prop_samp *lin=&(tr->ds.li->wavs);

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
  res=makesample(&tr->wavs,samp,lin,
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
  fromwav.i=WNU_O_WLU/(wsamp->f+tr->m);

  //convert from wavelength minimum
  fromwav.f=WNU_O_WLU/(wsamp->i-tr->m);

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

  //'.wnm' is the margin of not trusted wavenumbers, if it is not given
  if(trh->na&=TRH_WNM && trh->wnm>0)
    transitaccepthint(tr->wnm, trh->wnm, trh->na, TRH_WNM);
  else
    tr->wnm=tr->m*fromwav.f*fromwav.f/WNU_O_WLU;

  //make the sampling
  res=makesample(&tr->wns,nsamp,&fromwav,trh->na,
		 TRH_WNSFT,tr->wnm);

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
  int nrad,neiso=tr->n_e,niso=tr->n_i,ndb=tr->n_db;

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
  transitASSERT(atms->rads.n<1||!ndb||!niso||!neiso,
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
  isovt->d=(PREC_ATM  *)calloc(nrad*neiso,sizeof(PREC_ATM ));
  isovt->q=(PREC_ATM  *)calloc(nrad*neiso,sizeof(PREC_ATM ));
  isovt->c=(PREC_CS   *)calloc(nrad*neiso,sizeof(PREC_CS  ));
  for(i=1;i<neiso;i++){
    isovt[i].z=isovt->z+i*nrad;
    isovt[i].d=isovt->d+i*nrad;
    isovt[i].q=isovt->q+i*nrad;
    isovt[i].c=isovt->c+i*nrad;
  }
  atmt->t=(PREC_ATM *)calloc(nrad,sizeof(PREC_ATM));

  /* TD: interpolation */
  //interpolate temperature values according to radius
  resamplex(tr->fl,nrad,rad->v,rsamp->n,rsamp->v);
  resampley(tr->fl,1,
	    atms->atm.t,atmt->t);

  //Now for the isotope. 
  //First, for the isotopes that were added by getatm() (extended). We
  //can use same x axis as for atmospheric sample, because they were set
  //in the same file!
  //Density for all the isotopes
  /* TD-BUG: Find out why it fails if I take the brackets away. */
  for(i=0;i<neiso;i++){
    resampley(tr->fl,2,
	      atms->isov[i].d,isovt[i].d,
	      atms->isov[i].q,isovt[i].q);
  }

  //and cross section of only the extended isotopes. non extended go
  //below. This is used in calculating lorenz width
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
      transitASSERT(iso1db+j>niso-1,
		    "trying to reference an isotope (%i) outside\n"
		    "the extended limit (%i)\n"
		    ,iso1db+j,niso-1);
      resampley(tr->fl,2,
		isovs[j].z,isovt[iso1db+j].z,
		isovs[j].c,isovt[iso1db+j].c);
    }
  }

  if(res>=0)
    tr->pi|=TRPI_MAKERAD;
  return res;
}

