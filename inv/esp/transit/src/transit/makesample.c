/*
 * makesample.c
 * makesample.txc - create array sampling after initial, final and spacing
 *                  parameters. Component of the Transit program.
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

   @returns 1|2  for modified initial|final value
             0 if nothing was changed but there is a sampled array
	    -1 if hinted initial is bigger than maximum allowed.
	    -2 if hinted final is smaller than minimum allowed.
	    -3 if accepted initial value is greater or equal to final one.
	    -4 if both spacing and number of elements were hinted.
	    -5 if none or both of spacing or number of elements were in
               the referenced.
	    -6 Not valid oversampling was given by ref when requested.
*/
int
makesample(prop_samp *samp,	/* Resulting sampled data */
	   prop_samp *hint,	/* Proposed sampling */
	   prop_samp *ref,	/* Reference values */
	   const long fl,
	   const float margini,
	   const float marginf)
{
  //'res' is the returned value
  //'n' and 'v' are auxiliary variables to produced sampling array
  int res=0,n;
  PREC_RES *v;
  double osd,si;
  _Bool nhint,dhint;

  //The ratio of what is acceptable to exceed as upper value and not
  //removing the last bin.
  const double okfinalexcess=1e-8;

  //Check multiplicative factor
  if(hint->fct<=0)
    samp->fct=ref->fct;
  else
    samp->fct=hint->fct;

  //check initial value
  if(hint->i<=0||hint->i<ref->i+margini){
    samp->i=ref->i+margini;
    res|=0x1;
  }
  else if(hint->i>ref->f-marginf){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Hinted initial value for %s sampling, is bigger than\n"
		 " maximum allowed final value %.8g. Consider final\n"
		 "margin %.8g\n"
		 ,TRH_NAME(fl),ref->f-marginf,marginf);
    return -1;
  }
  else
    samp->i=hint->i;
  si=samp->i;

  //check final value
  if(hint->f<=0||hint->f>ref->f-marginf){
    samp->f=ref->f-marginf;
    res|=0x2;
  }
  else if(hint->f<ref->i+margini){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Hinted final value for %s sampling is smaller than\n"
		 "minimum allowed initial value %.8g.\n"
		 "Consider initial margin %.8g\n"
		 ,TRH_NAME(fl),ref->i+margini,margini);
    return -2;
  }
  else
    samp->f=hint->f;

  //check that resultant range makes sense
  if(samp->f<=si){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Initial accepted sampling value (%g) is greater or\n"
		 "equal than final accepted sample value(%g).\n"
		 "%s was being hinted\n"
		 ,si,samp->f,TRH_NAME(fl));
    return -3;
  }

  nhint=hint->n>0;
  dhint=hint->d>0;

  transitprint(21,verblevel,
	       "Flags: 0x%lx    hint.d:%g   hint.n:%li\n"
	       ,fl,hint->d,hint->n);
  //check that only one of spacing or number of elements field have been
  //hinted
  if(nhint&&dhint){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Both spacing(%g) and number of elements(%i) has\n"
		 "been hinted. That doesn't makes sense!. %s was being sampled.\n"
		 ,hint->d,hint->n,TRH_NAME(fl));
    return -4;
  }

  //if none has been hinted then use ref's
  if(!nhint&&!dhint){
    //If none of the ref exists then error
    if((ref->d<=0&&ref->n<=0)){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "Spacing(%g) and number of elements(%i) were\n"
		   "either both or none in the reference for %s sampling.\n"
		   "And yes, none was hinted.\n"
		   ,ref->d,ref->n,TRH_NAME(fl));
      return -5;
    }
    //if spacing exists then trust it
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
		     "values MIGHT have been modified.\n"
		     ,ref->n,TRH_NAME(fl)
		     ,ref->i,si,ref->f,samp->f);
      }
      samp->n=ref->n;
      samp->d=-1;
      samp->v=(PREC_RES *)calloc(ref->n,sizeof(PREC_RES));
      memcpy(samp->v,ref->v,ref->n*sizeof(PREC_RES));
      if(ref->o!=0)
	transiterror(TERR_WARNING,
		     "Fixed sampling array of length %i was referenced\n"
		     "But also oversampling was given (%li). Ignoring in\n"
		     "%s sampling\n"
		     ,samp->n,ref->o,TRH_NAME(fl));

      //return any possible modification
      return res;
    }
  }
  //else if spacing was hinted, then it has to be positive at this point
  else if(dhint){
    transitASSERT(hint->d<=0,
		  "OOPS!, logic test 1 failed in %s's makesample()!!\n"
		  ,TRH_NAME(fl));
    samp->d=hint->d;
  }
  //otherwise we have a positive hinted n
  else{
    transitASSERT(hint->n<=0,
		  "OOPS!, logic test 2 failed in %s's makesample()!!\n"
		  ,TRH_NAME(fl));
    samp->n=hint->n;
    samp->d=-1;
    samp->v=(PREC_RES *)calloc(hint->n,sizeof(PREC_RES));
    memcpy(samp->v,hint->v,hint->n*sizeof(PREC_RES));
    if(hint->o!=0)
      transiterror(TERR_WARNING,
		   "Fixed sampling array of length %li was hinted\n"
		   "But also oversampling was given (%li). Ignoring in\n"
		   "%s sampling\n"
		   ,samp->n,hint->o,TRH_NAME(fl));
    return res;
  }

  //At this points make the sampling if it was not given as an
  //array.
  n=samp->n=((1.0+okfinalexcess)*samp->f - si)/samp->d+1;

  //if there is an oversampling, check whether a value is not hinted
  if(hint->o<=0){
    //if so, check if we have a valid ref or error
    if(ref->o<=0){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "Not valid oversampling in the reference for\n"
		   "%s sampling.\n"
		   ,TRH_NAME(fl));
      return -6;
    }
    samp->o=ref->o;
  }
  else
    samp->o=hint->o;

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
		 ,TRH_NAME(fl),samp->d);

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
  prop_samp *samp=&tr->ds.th->wavs;
  prop_samp *lin=&(tr->ds.li->wavs);

  transitcheckcalled(tr->pi,"makewavsample",1,
		     "chkrange",TRPI_CHKRNG);

  if(samp->d<=0&&samp->n<=0){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Spacing or number must be hinted for wavelength,\n"
		 "cannot just guess them.\n"
		 );
    return -10;
  }

  //make the sampling
  res=makesample(&tr->wavs,samp,lin,
		 TRH_WAV,tr->m,tr->m);

  //set progress indicator if sampling was successful and return status
  if(res>=0)
    tr->pi|=TRPI_MAKEWAV;
  return res;
}


/* \fcnfh
 Calls makesample with the appropiate parameters and set the flags for a
 wavenumber sampling

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

  //use wavelength's oversampling and multiplicative factor
  fromwav.o=wsamp->o;
  fromwav.fct=wsamp->fct;

  //don't give a fixed array.
  fromwav.n=-1;

  //convert from wavelength maximum
  fromwav.i=WNU_O_WLU/wsamp->f;

  //convert from wavelength minimum
  fromwav.f=WNU_O_WLU/wsamp->i;

  //set margin. If not given take it from wavelength's (the bigger side)
  if(trh->wnm>0)
    tr->wnmf=tr->wnmi=trh->wnm;
  else{
    tr->wnmf=tr->m*fromwav.f*fromwav.f/WNU_O_WLU;
    tr->wnmi=tr->m*fromwav.i*fromwav.i/WNU_O_WLU;
  }

  //set spacing such that the wavenumber grid has the same number of
  //points as the wavelength grid. Then change reference initial and
  //final point to include margin
  fromwav.d=(fromwav.f-fromwav.i)/((wsamp->n-1)/wsamp->o);
  fromwav.f+=tr->wnmf;
  fromwav.i-=tr->wnmi;

  //make the sampling
  res=makesample(&tr->wns,nsamp,&fromwav,
		 TRH_WN,tr->wnmi,tr->wnmf);

  //set progress indicator if sampling was successful and return status
  if(res>=0)
    tr->pi|=TRPI_MAKEWN;
  return res;
}


/* \fcnfh
 Calls makesample with the appropiate parameters and set the flags for
 an impact parameter sampling

 @returns makesample() output
*/
int makeipsample(struct transit *tr)
{
  //'res' will be the result status
  int res;
  struct transithint *trh=tr->ds.th;
  prop_samp *usamp=&trh->ips;
  prop_samp *rsamp=&tr->rads;

  transitcheckcalled(tr->pi,"makeipsample",1,
		     "makeradsample",TRPI_MAKERAD);

  //make the sampling taking as reference the radius sampling
  res=makesample(&tr->ips,usamp,rsamp,
		 TRH_IPRM,0,0);

  //set progress indicator if sampling was successful and return status
  if(res>=0)
    tr->pi|=TRPI_MAKEIP;
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
  struct isotopes *iso=tr->ds.iso;
  int nrad,neiso=iso->n_e,niso=iso->n_i,ndb=iso->n_db;

  prop_isov *isovs;

  struct atm_data *atms=tr->ds.at;
  prop_samp *rsamp=&atms->rads;

  struct lineinfo *li=tr->ds.li;
  prop_isov *isovt=iso->isov;
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
    rad->n=1;
    rad->i=rsamp->i;
    rad->f=rsamp->f;
    rad->fct=rsamp->fct;
    rad->d=-1;
    rad->v=(PREC_RES *)calloc(1,sizeof(PREC_RES));
    rad->v[0]=rsamp->v[0];
    //return succes as it would have makesample()
    res=0;
    /* TD: warn that hinted values are going to be useless */
  }
  //If there is more than one atmospheric point
  else{
    //If no initial value is hinted, we take atmosphere minimum and
    //maximum.\par
    //do the sampling
    res=makesample(rad,&tr->ds.th->rads,rsamp,
		   TRH_RAD,0,0);
  }
  nrad=rad->n;

  //Allocate arrays that will receive the interpolated data
  isovt->d=(PREC_ATM  *)calloc(nrad*neiso,sizeof(PREC_ATM ));
  isovt->q=(PREC_ATM  *)calloc(nrad*neiso,sizeof(PREC_ATM ));
  isovt->c=(PREC_CS   *)calloc(nrad*niso, sizeof(PREC_CS  ));
  isovt->z=(PREC_ZREC *)calloc(nrad*niso, sizeof(PREC_ZREC));
  for(i=1;i<neiso;i++){
    isovt[i].d=isovt->d+i*nrad;
    isovt[i].q=isovt->q+i*nrad;
    if(i<niso){
      isovt[i].c=isovt->c+i*nrad;
      isovt[i].z=isovt->z+i*nrad;
    }
  }
  atmt->t=(PREC_ATM *)calloc(nrad,sizeof(PREC_ATM));
  atmt->p=(PREC_ATM *)calloc(nrad,sizeof(PREC_ATM));

  /* TD: interpolation */
  //interpolate temperature values according to radius
  resamplex(tr->fl,rsamp->n,rsamp->v,nrad,rad->v);
  resampley(tr->fl,2,
	    atms->atm.t,atmt->t,
	    atms->atm.p,atmt->p);

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

  //Second, non-extended isotopes (isotopes coming from the line database):
  //We have to go to each database separately
  for(i=0;i<ndb;i++){
    //position in the first isotope of the database
    iso1db=iso->db[i].s;
    isovs=li->isov+iso1db;

    //interpolate variable isotope info respect to temperature
    resamplex(tr->fl,li->db[i].t,li->db[i].T,nrad,atmt->t);
    for(j=0;j<iso->db[i].i;j++){
      transitASSERT(iso1db+j>niso-1,
		    "trying to reference an isotope (%i) outside\n"
		    "the extended limit (%i)\n"
		    ,iso1db+j,niso-1);
      resampley(tr->fl,2,
		isovs[j].z,isovt[iso1db+j].z,
		isovs[j].c,isovt[iso1db+j].c);
    }
  }

  //set progress indicator and return
  if(res>=0)
    tr->pi|=TRPI_MAKERAD;
  return res;
}


/* \fcnfh
   Prints sample info for a structure
*/
static void
printsample(FILE *out,	/* File pointer to write out */
	    prop_samp *samp,	/* Sample strucuture to print */
	    char *desc,		/* Description */
	    long fl)		/* Various flag */
{
  int i;

  fprintf(out,
	  "############################\n"
	  "   %-12s Sampling\n"
	  "----------------------------\n"
	  ,desc);

  fprintf(out,"Factor to cgs units: %g\n",samp->fct);

  fprintf(out,"Initial value: %g\nFinal value: %g\n",samp->i,samp->f);

  fprintf(out,"Spacing: %g\n",samp->d);
  if(!(fl&TRF_NOOVERSAMP))
    fprintf(out,"Oversample: %i\n",samp->o);

  fprintf(out,"Number of elements: %li\n",samp->n);

  if(!(fl&TRF_NOVALUE)){
    fprintf(out,"Values: ");
    for(i=0;i<samp->n;i++)
      fprintf(out," %8.3g",samp->v[i]);
    fprintf(out,"\n");
  }

}

/* \fcnfh
   Prints each of the sample data
*/
void
outsample(struct transit *tr)
{
  FILE *out=stdout;

  char *filename=tr->f_outsample;

  if(strcmp(filename,"-")==0)
    if((out=fopen(filename,"w"))!=NULL){
      transiterror(TERR_WARNING,
		   "Cannot open file '%s' for writing sampling data...\n"
		   ,filename);
      return;
    }

  printsample(out,&tr->wns,"Wavenumber",TRF_NOVALUE);
  printsample(out,&tr->wavs,"Wavelength",TRF_NOVALUE);
  printsample(out,&tr->rads,"Radius",TRF_NOOVERSAMP);
  printsample(out,&tr->ips,"Impact parameter",0);

  fclose(out);

}


#ifdef DBGSAMPLE
int
main(int argc, char *argv[])
{
  prop_samp lim,hint={0,0,0,0,0},res;
  float mar=0;
  int i;

  if(argc<5){
    fprintf(stderr,"Syntax:\n"
	    "\tdbgsample <ini> <fin> <delt> <oversampling> [<margin>]\n");
    exit(0);
  }

  lim.i=atof(argv[1]);
  lim.f=atof(argv[2]);
  lim.d=atof(argv[3]);
  lim.o=atof(argv[4]);
  lim.n=0;
  if(argc==6)
    mar=atof(argv[5]);

  i=makesample(&res,&hint,&lim,0,0,mar);

  fprintf(stderr,
	  "Makesample returned %i\n"
	  "Initial %g, final %g, delta %g, oversamp %i, number %li, margin %g\n"
	  ,i,res.i,res.f,res.d,res.o,res.n,mar);

  if(i<0)
    exit(1);

  for(i=0;i<res.n;i++)
    fprintf(stderr," rad(%i): %g\n",i,res.v[i]);

}

#endif /* debug */
