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
	 PREC_RES *dt, long nrad, short iso, long wn,gsl_interp_accel *acc);

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

  //Need at least three radius to calculate a spline interpolation.
  if(inn<3)
    transiterror(TERR_SERIOUS,
		 "tau(): At least three impact parameters points are required!.\n"
		 );

  gsl_interp_accel *acc=gsl_interp_accel_alloc();
  //for each wavenumber
  for(wi=0;wi<wnn;wi++){
    t=tau.t[wi];

    //For each resultant impact parameter
    for(ii=inn-1;ii>=0;ii--){
      if((t[ii]=totaltau(bb[ii],r,n,e,dt,inn,tau.iso,wi,acc))>tau.toomuch){
	tau.first[wi]=ii;
	break;
      }
    }
  }
  gsl_interp_accel_free(acc);

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
	 long wn,		/* wavenumber looked */
	 gsl_interp_accel *acc)	/* accelerating pointer */
{
  PREC_RES res;
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

  //get bin value 'rs' such that r0 is between rad[rs-1] inclusive and
  //rad[rs] exclusive.
  rs=binsearch(rad,0,nrad-1,r0)+1;

  //return 0 optical depth if it goes to less than three layers (this is
  //for the spline integration to work, an alternative method could be
  //installed instead).
  if(rs>nrad-3){
    /* TD: Install an alternative integration method for less than 3
       points */
    return 0;
  }

  if(rs<0)
    transiterror(TERR_CRITICAL,
		 "Closest approach value(%g) is outside sampled radius\n"
		 "range(%g - %g)\n"
		 ,r0,rad[0],rad[nrad-1]);

  //A fraction 'analiticfrac' of the integration near the closest
  //appraoach is calcualated analitically, otherwise, I get a division
  //by zero. In formula\par
  //\[
  //\tau_{\wn}(\rho)=
  //\underbrace{
  //2\extc_{\wn}r_0\sqrt{2r_0\delta r} 
  //}_{\mathrm{analitic}} +
  //\underbrace{
  //2\int_{r_0+\delta r}^{\infty}
  //\frac{\extc_{\wn}~n~r}{\sqrt{n^2r^2-\rho^2}}\dd r
  //}_{\mathrm{numerical}}
  //\]\par
  //First for the analitical part of the integral
  PREC_RES analiticfrac=(rad[rs]-r0);
  res=ex[i][iso][wn]*sqrt(2*r0*analiticfrac);

  //And now for the numerical integration.\par
  //This part currently depends on a proper installation of GSL library
#ifdef _USE_GSL
  for(i=rs;i<nrad;i++){
    r0a=b/refr[i]/rad[i];
    transitASSERT(r0a>1,
		  "Oops! assert condition not met, b/(nr)=%g",r0a);

    dt[i]=ex[i][iso][wn]/sqrt(1-r0a*r0a);
  }


  acc->cache = 0;
  acc->hit_count = 0;
  acc->miss_count = 0;
  gsl_spline *spl=gsl_spline_alloc(gsl_interp_cspline,nrad-rs);
  gsl_spline_init(spl,rad+rs,dt+rs,nrad-rs);
  res+=gsl_spline_eval_integ(spl,rad[rs],rad[nrad-1],acc);
  gsl_spline_free(spl);
  return 2*res;

  //Without {\bf GSL} is currently not implemented. Output is dependent in an
  //appropiate installation of those libraries
#else
#error computation of tau() without GSL is not implemented
#endif
}
