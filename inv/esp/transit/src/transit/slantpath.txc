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

static inline PREC_RES
modulation1 (PREC_RES *tau,
	     long first,
	     double toomuch,
	     prop_samp *ip,
	     struct geometry *sg,
	     gsl_interp_accel *acc) __attribute__((always_inline));

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
   interaction. Calculates ratio in-transit over out-transit for the
   simplest expression of modulation.

   @returns modulation obtained
*/
static PREC_RES
modulationperwn (PREC_RES *tau,
		 long first,
		 double toomuch,
		 prop_samp *ip,
		 struct geometry *sg,
		 int exprlevel,
		 gsl_interp_accel *acc)
{
  PREC_RES res;

  switch(exprlevel){
  case 1:
    res=modulation1(tau,first,toomuch,ip,sg,acc);
    break;
  default:
    res=0;
    transiterror(TERR_CRITICAL,
		 "slantpath:: modulationperwn:: Level %i of detail\n"
		 "has not been implemented to compute modulation\n"
		 ,exprlevel);
  }
  return res;
}


/* \fcnfh
   Computes most basic modulation scheme, obtained when there is no limb
   darkening or emitted flux 

   @returns modulation obtained
*/
static inline PREC_RES
modulation1 (PREC_RES *tau,
	     long first,
	     double toomuch,
	     prop_samp *ip,
	     struct geometry *sg,
	     gsl_interp_accel *acc)
{
  //general variables
  PREC_RES res;
  double srad=sg->starrad*sg->starradfct;
  double *rinteg,*azinteg,tt;

  //Impact parameter and azimuth variables
  long ipn=ip->n,azn;
  PREC_RES ipd=ip->d*ip->fct,*ipv,ipvv;
  PREC_RES azz,*az,azd;

  //star centered coordinates and counters
  double starx,stary;
  long i,j;

  //allocate enough azimuthal bins for the outermost (biggest) radii,
  //others will necessarily be smaller.
  azd=2*PI*ip->v[ipn-1]*ip->fct/ipd;
  azn=(long)azd+1;
  azinteg=(double *)alloca(azn*sizeof(double));
  rinteg=(double *)alloca(ipn*sizeof(double));
  ipv=(double *)alloca(ipn*sizeof(double));
  az=(double *)alloca(azn*sizeof(double));

  //this function calculates 1 minus the ratio of in-transit over out-of-transit
  //expresion for the simplest case, which is given by
  //\[
  //M_{\lambda}^{(0)}=
  //\frac{1}{\pi R_M^2}\int_{R<R_M}\ee^{-\tau( b,\xi)} R \dd R\dd\phi
  //\]\par
  //Let's integrate; for each of the planet's layer starting from the
  //outermost until the closest layer
  first--;
  for(i=ipn-1;i>first;i--){
    //take azimuthal spacing equal to the radial spacing. Add one bin so
    //that integration can be done until $2\pi$
    ipvv=ipv[i]=ip->v[i]*ip->fct;
    azd=ipd/ipvv;
    azn=2*PI/azd +1;
    tt=tau[i];

    //fill azimuthal integrand
    for(j=0;j<azn;j++){
      azz=az[j]=j*azd;
      starx=(sg->x+ipvv*sin(azz))/sg->starrad;
      stary=(sg->y+ipvv*cos(azz))/sg->starrad;
      azinteg[j]=exp(-tt)*ipvv;
    }

    //feed gsl.
#ifdef _USE_GSL
    acc->cache = 0;
    acc->hit_count = 0;
    acc->miss_count = 0;
    gsl_spline *spl=gsl_spline_alloc(gsl_interp_cspline,azn);
    gsl_spline_init(spl,az,azinteg,azn);
    rinteg[i]=gsl_spline_eval_integ(spl,0,2*PI,acc);
    gsl_spline_free(spl);

    //Without {\bf GSL} is currently not implemented. Output is dependent in an
    //appropiate installation of those libraries
#else
# error computation of modulation() without GSL is not implemented
#endif
  }
  //fill two more lower part bins with 0. Only two to have a nice ending
  //spline and not unnecessary values.
  first-=2;
  if(first<-1) first=-1;
  for(;i>first;i--){
    ipv[i]=ip->v[i]*ip->fct;
    rinteg[i]=0;
  }

  //set right position from the array beginning and number of items,
  //check that we have enough layers to integrate, I can't imagine now
  //how can that check be failed at this point.
  first++;
  ipn-=first;
  if(ipn<3)
    transiterror(TERR_CRITICAL,
		 "Condition failed, less than 3 items (only %i) for radial\n"
		 "integration.\n"
		 ,first);

  //integrate in radii
#ifdef _USE_GSL
  acc->cache = 0;
  acc->hit_count = 0;
  acc->miss_count = 0;
  gsl_spline *spl=gsl_spline_alloc(gsl_interp_cspline,ipn);
  gsl_spline_init(spl,ipv+first,rinteg+first,ipn);
  res=gsl_spline_eval_integ(spl,ipv[first],ipv[first+ipn-1],acc);
  gsl_spline_free(spl);

  //or err without GSL
#else
# error computation of modulation() without GSL is not implemented
#endif

  /* TD: Add  unblocked area of the star */
  //substract the total area blocked by the planet. This is from the
  //following
  //\begin{align}
  //1-M=&1-\frac{\int_0^{r_p}\int\ee^{-\tau}\dd \theta r\dd r
  //             +\int_{r_p}^{R_s}\dd A}
  //            {\pi R_s^2}\\
  //   =&\frac{\int_0^{r_p}\int\ee^{-\tau}\dd \theta r\dd r
  //           -Area_{planet}}
  //          {\pi R_s^2}
  //\end{align}
  res-=PI*ipv[first+ipn-1]*ipv[first+ipn-1];

  //divide by area of the star
  res/=PI*srad*srad;

  return res;
}



const transit_ray_solution slantpath =
  {
    "Slant Path",		/* Name of the solution */
    "slantpath.c",		/* This file name */
    "1.4",			/* GSL version */
    1,				/* Equispaced impact parameter requested? */
    &totaltau,			/* per impact parameter and per
				   wavenumber value computation */
    &modulationperwn,		/* per wavenumber value computation in
				   its simplest expression */
    1				/* Number of levels of modulation detail
				   as it can be handled by the above
				   fucntion */
  };
