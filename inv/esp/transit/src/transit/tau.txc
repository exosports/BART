/*
 * tau.c
 * tau.txc - Finds the optical depth. Component of the Transit program.
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

static inline PREC_RES
totaltau(PREC_RES b, PREC_RES *rad, PREC_RES *refr, PREC_RES ***ex,
	 PREC_RES *dt, long nrad, short iso, long wn);

/* \fcnfh
   Calculates optical depth as a function of radii for a spherical
   symmetric planet

   @returns 0 on success
 */
int
tau(struct transit *tr)
{
  /* Different structures */
  static struct optdepth tau;
  tr->ds.tau=&tau;
  prop_samp *rad=&tr->rads;
  prop_samp *wn=&tr->wns;
  prop_samp *ip=&tr->ips;
  PREC_RES ***e=tr->ds.ex->e;

  //index, initial and final values
  long ii,wi;
  long inn,rnn,wnn;
  //'bb' is the impact parameter, while 'n' is the index of refraction,
  //'w' is the wavenumber, 't' is the tau as function of impact
  //parameter, 'r' is the radius
  PREC_RES *bb=ip->v;
  PREC_RES *n=tr->ds.ir->n;
  PREC_RES *r=rad->v;
  PREC_RES *t;

  transitcheckcalled(tr->pi,"tau",2,
		     "idxrefrac",TRPI_IDXREFRAC,
		     "extwn",TRPI_EXTWN
		     );

  transitacceptflag(tr->fl,tr->ds.th->fl,TRU_TAUBITS);

  //set tau structures' value
  tau.toomuch=50;
  if(tr->ds.th->na&TRH_TOOMUCH&&tr->ds.th->toomuch>0)
    tau.toomuch=tr->ds.th->toomuch;
  tau.iso=0;
  if(tr->ds.th->na&TRH_TAUISO&&tr->ds.th->tauiso>=0&&tr->ds.th->tauiso<tr->n_i)
    tau.iso=tr->ds.th->tauiso;
  tau.first=(long *)calloc(wn->n,sizeof(long));
  tau.t=(PREC_RES **)calloc(wn->n,sizeof(PREC_RES *));
  tau.t[0]=(PREC_RES *)calloc(wn->n*ip->n,sizeof(PREC_RES));
  for(ii=1;ii<wn->n;ii++)
    tau.t[ii]=tau.t[0]+ii*ip->n;

  //final index or number of elements
  wnn=wn->n;
  inn=ip->n;
  rnn=rad->n;

  //to temporarily store a per radius info
  PREC_RES dt[wnn];

  //for each wavenumber
  for(wi=0;wi<wnn;wi++){
    t=tau.t[wi];

    //For each resultant impact parameter
    for(ii=inn-1;ii>=0;ii--){
      if((t[ii]=totaltau(bb[ii],r,n,e,dt,rnn,tau.iso,wi))>tau.toomuch){
	tau.first[wi]=ii;
	break;
      }
    }
  }

  tr->pi|=TRPI_TAU;
  return 0;
}


/* \fcnfh
 Computes optical depth at a given impact parameter

 @returns optical depth
*/
static inline PREC_RES
totaltau(PREC_RES b,		/* impact parameter */
	 PREC_RES *rad,		/* radius array */
	 PREC_RES *refr,	/* refractivity index */
	 PREC_RES ***ex,	/* extinction[rad][iso][nwn] */
	 PREC_RES *dt,		/* differential optical depth [rad] */
	 long nrad,		/* number of radii elements */
	 short iso,		/* isotope chosen */
	 long wn)		/* wavenumber looked */
{
  PREC_RES r0a=b;
  PREC_RES r0=0;
  int i;
  int maxiterations=50;
  int rs;

  //Look for closest approach radius
  i=0;
  while(1){
    r0=b/lineinterp(r0a,rad,refr,nrad);
    if(r0==r0a)
      break;
    if(i++>maxiterations)
      transiterror(TERR_CRITICAL,
		   "Maximum iterations(%i) reached while looking for\n"
		   "r0. Convergence not reached (%.6g!=%.6g)\n"
		   ,maxiterations,r0,r0a);
    r0a=r0;
  }

  //get bin value 'rs' such that r0 is between rad[i] inclusive and
  //rad[i+1] exclusive
  rs=binsearch(rad,0,nrad-1,r0);
  transitASSERT(rs==-4,
		"Number of radius sample(%i) has to be at least 1!\n"
		,nrad);
  if(rs<0)
    transiterror(TERR_CRITICAL,
		 "Closest approach value(%g) is outside sampled radius\n"
		 "range(%g - %g)\n"
		 ,r0,rad[0],rad[nrad-1]);

  //Now integrate the ray's path from r0 to rad->v[nrad-1] using
  //\[
  //\tau_{\wn}(\rho)=2\int_{r_0}^{\infty}
  //\frac{\extc_{\wn}~n~r}{\sqrt{n^2r^2-\rho^2}}\dd r
  //\]
#ifdef _USE_GSL
  for(i=rs;i<nrad;i++){
    r0a=b/refr[i]/rad[i];
    dt[i]=2*ex[i][iso][wn]/sqrt(1-r0a*r0a);
  }

  //This part depends on a proper installation of GSL library
  gsl_interp_accel *acc=gsl_interp_accel_alloc();
  gsl_spline *spl=gsl_spline_alloc(gsl_interp_cspline,nrad-rs);
  gsl_spline_init(spl,rad+rs,dt+rs,nrad-rs);
  r0a=gsl_spline_eval_integ(spl,r0,rad[nrad-1],acc);
  gsl_spline_free(spl);
  gsl_interp_accel_free(acc);
  return r0a;

#else
  //THIS IS NOT IMPLEMENTED, CURRENT OUTPUT IS DEPENDENT IN AN
  //APPROPIATE GSL INSTALLATION
#error computation of tau() without GSL is not implemented
#endif
}
