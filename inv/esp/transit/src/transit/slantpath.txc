/*
 * slantpath.c
 * slantpath.txc - Functions to handle a slant path problem. Component
 *                 of the transit program
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


/***** Warning:
       To speed up, slantpath.txc assumed version 1.4 of gsl
       utilities. If not, structure gsl_interp_accel might have
       different components, which will produce erroneus assumptions in
       totaltau(). 
***************/

/* \fcnfh
 Computes optical depth at a given impact parameter, note that b needs
 to be given in units of 'rad' and the result needs to be multiplied by
 the units 'rad' to be real.

 @returns $\frac{tau}{units_{rad}}$ returns optical depth divided by units
                                    of 'rad'
*/
static inline PREC_RES
totaltau(PREC_RES b,		/* impact parameter, in units of rad. */
	 PREC_RES *rad,		/* radius array */
	 PREC_RES *refr,	/* refractivity index */
	 PREC_RES ***ex,	/* extinction[rad][iso][nwn] */
	 long nrad,		/* number of radii elements */
	 short iso,		/* isotope chosen */
	 long wn,		/* wavenumber looked */
	 PREC_RES *dt,		/* differential optical depth [rad].
				   Auxiliary array */
	 gsl_interp_accel *acc)	/* accelerating pointer. Auxiliary array
				   */
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

  //get bin value 'rs' such that r0 is between rad[rs-1] inclusive
  //and rad[rs] exclusive.
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
  res=ex[rs-1][iso][wn]*sqrt((analiticfrac+2*r0)*analiticfrac)/refr[rs-1];

  //And now for the numerical integration. Set the variables
  for(i=rs;i<nrad;i++){
    r0a=b/refr[i]/rad[i];
    transitASSERT(r0a>1,
		  "Oops! assert condition not met, b/(nr)=%g",r0a);

    dt[i]=ex[i][iso][wn]/sqrt(1-r0a*r0a);
  }

  //and integrate!.\par
  //This part currently depends on a proper installation of GSL library
#ifdef _USE_GSL
  acc->cache = 0;
  acc->hit_count = 0;
  acc->miss_count = 0;
  gsl_spline *spl=gsl_spline_alloc(gsl_interp_cspline,nrad-rs);
  gsl_spline_init(spl,rad+rs,dt+rs,nrad-rs);
  res+=gsl_spline_eval_integ(spl,rad[rs],rad[nrad-1],acc);
  gsl_spline_free(spl);

#else
  //Without {\bf GSL} is currently not implemented. Output is dependent in an
  //appropiate installation of those libraries
# error computation of totaltau() without GSL is not implemented
#endif

  return 2*res;

}

/* \fcnfh 
   observable information as it would be seen before any telescope
   interaction

   @returns modulation obtained
*/
static inline PREC_RES
modulationperwn (PREC_RES *tau,
		 long first,
		 double toomuch,
		 prop_samp *ip,
		 struct geometry *sg,
		 gsl_interp_accel *acc)
{
  //general variables
  double srad=sg->starrad*sg->starradfct;
  double *integ,tt;
  PREC_RES res;

  //Impact parameter and azimuth variables
  long ipn=ip->n,azn;
  PREC_RES ipd=ip->d,ipv;
  PREC_RES azz,*az,azd;

  //star centered coordinates and counters
  double starx,stary;
  long i,j;
  
  //allocate enough azimuthal bins for the outermost (biggest) radii,
  //others will necessarily be smaller.
  azd=2*PI*ip->v[ipn-1]*ip->fct/ipd;
  azn=(long)azd+1;
  integ=(double *)alloca(azn*sizeof(double));
  az=(double *)alloca(azn*sizeof(double));

  //for each of the planet's layer starting from the outermost until the
  //closest layer
  first--;
  for(i=ipn-1;i>first;i--){
    //take azimuthal spacing equal to the radial spacing. Add one bin so
    //that integration can be done until 2PI
    ipv=ip->v[i]*ip->fct;
    azd=ipd/ipv;
    azn=2*PI/azd +1;
    tt=tau[i];

    //fill azimuthal integrand
    for(j=0;j<azn;j++){
      azz=az[j]=j*azd;
      starx=(sg->x+ipv*sin(azz))/sg->starrad;
      stary=(sg->y+ipv*cos(azz))/sg->starrad;
      integ[j]=starvariation(starx,stary,srad)*
	exp(-tt)*ipv;
    }

    //feed gsl.
#ifdef _USE_GSL
    acc->cache = 0;
    acc->hit_count = 0;
    acc->miss_count = 0;
    gsl_spline *spl=gsl_spline_alloc(gsl_interp_cspline,azn);
    gsl_spline_init(spl,az,integ,azn);
    res=gsl_spline_eval_integ(spl,0,2*PI,acc);
    gsl_spline_free(spl);

#else
    //Without {\bf GSL} is currently not implemented. Output is dependent in an
    //appropiate installation of those libraries
# error computation of modulation() without GSL is not implemented
#endif

  }

  /* TD: Add normalized unblocked area of the star */
  res+=0;

  return 0;
}


const transit_ray_solution slantpath =
  {
    "Slant Path",		/* Name of the solution */
    "slantpath.c",		/* This file name */
    "1.4",			/* GSL version */
    1,				/* Equispaced impact parameter requested? */
    &totaltau,			/* per impact parameter and per
				   wavenumber value computation */
    &modulationperwn		/* per wavenumber value computation */
  };
