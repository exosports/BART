/*
 * slantpath.c   - Functions to handle a slant path problem. Component
 *                 of the transit program
 *
 * Copyright (C) 2003-2006 Patricio Rojo (pato@astro.cornell.edu)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of version 2 of the GNU General 
 * Public License as published by the Free Software Foundation.
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

/*
List of functions defined here:
static inline PREC_RES 
totaltau(PREC_RES b, PREC_RES *rad, PREC_RES *refr, PREC_RES *ex, long nrad,
         int exprlevel) 

static PREC_RES
totaltau1(PREC_RES b, PREC_RES *rad, PREC_RES refr, PREC_RES *ex, long nrad)

static PREC_RES
totaltau2(PREC_RES b, PREC_RES *rad, PREC_RES *refr, PREC_RES *ex, long nrad)

static PREC_RES
modulationperwn(PREC_RES *tau, long last, double toomuch, prop_samp *ip,
                struct geometry *sg, int exprlevel)

static PREC_RES
modulation1(PREC_RES *tau, long last, double toomuch, prop_samp *ip,
            struct geometry *sg)

static inline PREC_RES
modulationm1(PREC_RES *tau, long last, double toomuch, prop_samp *ip,
             struct geometry *sg)                                            */

#include <transit.h>


/* \fcnfh
   Computes optical depth at a given impact parameter for a medium
   with constant index of refraction (no ray bending).

   The impact parameter, b, and radii rad should have the same units.
   The result needs to be multiplied by the units of 'rad' to be real.
   It can take nrad values of 1 or bigger. However, if 2 is given then
   'ex' and 'rad' need to have a referenceable element at position
   -1; i.e. rad[-1] and ex[-1] need to exist.

   Return: Optical depth divided by rad.fct:  \frac{tau}{units_{rad}} */
static PREC_RES
totaltau1(PREC_RES b,    /* Impact parameter         */
          PREC_RES *rad, /* Equispaced radius array  */
          PREC_RES refr, /* Refractivity index       */
          PREC_RES *ex,  /* Extinction[rad]          */
          long nrad){    /* Number of radii elements */
  int rs, i;     /* Auxilliary radius and for-loop indices    */
  PREC_RES res;  /* Half-path optical depth                   */
  PREC_RES x3[3], r3[3]; /* Auxiliary interpolation variables */

  /* Closest approach radius: */
  PREC_RES r0 = b/refr;
  /* FINDME: isn't this r0 always less than rad[0]??  Because this routine
             receives the array from r[last], not the entire array.        */

  /* Get the index rs, of the sampled radius immediately below or equal
     to r0 (i.e. rad[rs] <= r0 < rad[rs+1]):                            */
  if((rs = binsearch(rad, 0, nrad-1, r0)) == -5)
    return 0;  /* If it is the outmost layer                 */
  /* If some other error occurred:                           */
  else if(rs < 0)
    transiterror(TERR_CRITICAL,
                 "Closest approach value (%g) is outside sampled radius "
                 "range (%g - %g).\n", r0, rad[0], rad[nrad-1]);

  /* Move extinction and radius pointers to the rs-th element: */
  rad  += rs;
  ex   += rs;
  nrad -= rs;

  /* Calculate the extinction coefficient at r0 using a parabolic
     interpolation.  Temporarily store extinction(r0) and r0 at the rs-th
     position of ex and rad arrays: */
  const PREC_RES tmpex  = *ex;  /* Get a constant copy of the arrays */
  const PREC_RES tmprad = *rad;
  /* Interpolate and modify arrays values: */
  if(nrad==2) *ex = interp_parab(rad-1, ex-1, r0); 
  else        *ex = interp_parab(rad,   ex,   r0);
  *rad = r0;
  /* Handle special case of only two elements in array, make an intermediate
     point in between values: */
  if(nrad == 2){
    x3[0] = ex[0];
    x3[2] = ex[1];
    x3[1] = (ex[1]+ex[0])/2.0;
    r3[0] = rad[0];
    r3[2] = rad[1];
    r3[1] = (rad[0]+rad[1])/2.0;
    *rad = tmprad;
    *ex  = tmpex;
    rad  = r3;
    ex   = x3;
    nrad++;
  }
  const PREC_RES dr = rad[1]-rad[0];

  /* Convert to s spacing, i.e. distance along the path. Radius needs
     to be equispaced: */
  PREC_RES s[nrad]; /* Distance along ray path */
  const PREC_RES Dr  = rad[2]-rad[1];
  const PREC_RES cte = dr*(dr + 2*r0);
  s[0] = 0;

  for(i=1; i<nrad; i++)
    s[i] = sqrt(cte + (i-1)*Dr*(2.0*(r0+dr) + (i-1)*Dr) );

  /* Integrate Extinction along ray path: */
  /* Use spline if GSL is available along with at least 3 points: */
#ifdef _USE_GSL
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp       *spl = gsl_interp_alloc(gsl_interp_cspline, nrad);
  gsl_interp_init(spl, s, ex, nrad);
  res = gsl_interp_eval_integ(spl, s, ex, 0, s[nrad-1], acc);
  gsl_interp_free(spl);
  gsl_interp_accel_free(acc);
#else
#error non equispaced integration is not implemented without GSL
#endif /* _USE_GSL */

  /* Reset original values of arrays: */
  *ex  = tmpex;
  *rad = tmprad;

  /* Return: */
  return 2*(res);
}


/* \fcnfh
   Computes optical depth at a given impact parameter for a medium with
   variable index of refraction.

 The impact parameter b needs to be given in units of 'rad' and the result
 needs to be multiplied by the units 'rad' to be real.

 Return: $\frac{tau}{units_{rad}}$ returns optical depth divided by units
                                   of 'rad'                               */
static PREC_RES
totaltau2(PREC_RES b,     /* Impact parameter         */
          PREC_RES *rad,  /* Radius array             */
          PREC_RES *refr, /* Refractivity index       */
          PREC_RES *ex,   /* Extinction[rad]          */
          long nrad){     /* Number of radii elements */
  PREC_RES dt[nrad];
  PREC_RES r0a = b;
  PREC_RES r0 = 0;
  int rs, i=0;
  const int maxiterations=50;

  transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
               "This routine has not been successfully tested yet. Be "
               "critic of the result.\n");

  /* Look for closest approach radius: */
  while(1){
    r0 = b/lineinterp(r0a, rad, refr, nrad);
    if(r0==r0a)
      break;
    if(i++ > maxiterations)
      transiterror(TERR_CRITICAL,
                   "Maximum iterations (%i) reached while looking for "
                   "r0. Convergence not reached (%.6g!=%.6g).\n",
                   maxiterations, r0, r0a);
    r0a = r0;
  }

  /* Get bin value 'rs' such that r0 is between rad[rs-1] inclusive
     and rad[rs] exclusive: */
  /*If we are looking at the outmost layer, then return: */
  if((rs=binsearch(rad, 0, nrad-1, r0))==-5)
    return 0;
  /* If some other error occurred: */
  else if(rs<0)
    transiterror(TERR_CRITICAL,
                 "Closest approach value(%g) is outside sampled radius "
                 "range(%g - %g).\n", r0, rad[0], rad[nrad-1]);
  /* Move pointer to rs-th element: */
  rs++;

  /* A fraction 'analiticfrac' of the integration near the closest
     appraoach is calcualated analitically, otherwise, I get a division
     by zero in formula:
     \tau_{\nu}(\rho) =
     \underbrace{\frac{2e_{\nu}\rho}{n}
                  \left(\sqrt{\left(\frac{nr_1}{\rho}\right)^2-1} -
                        \sqrt{\left(\frac{nr_0}{\rho}\right)^2-1}\right) 
     }_{\rm{analitic}} +
     \underbrace{2\int_{r_1=r_0+\delta r}^{\infty}
                 \frac{e_{\nu}~n~r}{\sqrt{n^2r^2-\rho^2}}{\rm d} r
     }_{\rm{numerical}} */

  /* First for the analitical part of the integral: */
  PREC_RES res;
  if(ex[rs-1]==ex[rs])
    res = ex[rs] * r0 * (sqrt( rad[rs]*rad[rs] / (r0*r0) - 1));
  else{
    PREC_RES alpha = ( ex[rs] - ex[rs-1] ) / ( rad[rs] - rad[rs-1] );
    PREC_RES rm    = rad[rs];
    if(alpha<0)
      res = - alpha * (rm * sqrt( rm*rm -  r0*r0) - r0*r0 * 
                        log(sqrt( rm*rm / (r0*r0) - 1) + rm/r0 )) / 2.0;
    else
      res =   alpha * (rm * sqrt( rm*rm -  r0*r0) + r0*r0 * 
                        log(sqrt( rm*rm / (r0*r0) - 1) + rm/r0 ))  / 2.0;
  }

  /* And now for the numerical integration. Set the variables: */
  for(i=rs; i<nrad; i++){
    r0a = b / (refr[i]*rad[i]);
    transitASSERT(r0a>1, "Condition could not be asserted, b/(nr)=%g > 1.\n",
                  r0a);

    dt[i] = ex[i]/sqrt(1-r0a*r0a);
  }

  /* Integrate: */
  /* Use spline if GSL is available along with at least 3 points: */
#ifdef _USE_GSL
  if(nrad-rs>2){
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spl = gsl_spline_alloc(gsl_interp_cspline, nrad-rs);
    gsl_spline_init(spl, rad+rs, dt+rs, nrad-rs);
    res += gsl_spline_eval_integ(spl, rad[rs], rad[nrad-1], acc);
    gsl_spline_free(spl);
    gsl_interp_accel_free(acc);
  }
  /* Only integrate Trapezium if there are only two points available: */
  else
#endif /* _USE_GSL */
  /* Integrate Simpson-Trapezium if enough(w/o GSL) or not enough (w/ GSL)
     elements: */
  if(nrad-rs>1)
    res += integ_trasim(rad[1]-rad[0], dt+rs, nrad-rs);

  return 2*(res);
}


/* \fcnfh
 Wrapper function to calculate the optical depth at a given impact parameter

 Return: $\frac{tau}{units_{rad}}$ returns optical depth divided by units
                                    of 'rad'                           */
static inline PREC_RES
totaltau(PREC_RES b,      /* Differential impact parameter WRT maximum value */
         PREC_RES *rad,   /* Radius array               */
         PREC_RES *refr,  /* Refractivity index         */
         PREC_RES *ex,    /* Extinction[rad]            */
         long nrad,       /* Number of radii elements   */
         int exprlevel){  /* Expression level of detail */
  switch(exprlevel){
  case 1:
    return totaltau1(b, rad, *refr, ex, nrad);
    break;
  case 2:
    return totaltau2(b, rad,  refr, ex, nrad);
    break;
  default:
    transiterror(TERR_CRITICAL,
                 "slantpath:: totaltau:: Level %i of detail has not been "
                 "implemented to compute optical depth.\n", exprlevel);
    return 0;
  }
}


/* \fcnfh
   Calculate the transit's modulation at a given wavenumber for
   no-limb darkening nor emitted flux.

   Return: the transit's modulation:
   1 - in-transit/out-of-transit flux ratio (Equation 3.12):
    M_{\lambda} = \frac{1}{R_\star^2}\left(R^2 - 
                   2\int_{0}^{R} \exp^{-\tau_\lambda(r)} r\,{\rm d}r\right)  */
static PREC_RES
modulation1(PREC_RES *tau,        /* Optical depth array              */
            long last,            /* Index where tau = toomuch        */
            double toomuch,       /* Maximum optical depth calculated */
            prop_samp *ip,        /* Impact parameter                 */
            struct geometry *sg){ /* Geometry struct                  */
  /* General variables: */
  PREC_RES res;
  /* Stellar radius:    */
  double srad = sg->starrad*sg->starradfct;

  /* Impact parameter variables: */
  long ipn  = ip->n;
  long ipn1 = ipn-1;
  long i;

  /* Max overall tau, for the tr.ds.sg.transparent=True case */
  const PREC_RES maxtau = tau[last]>toomuch?tau[last]:toomuch;
  /* Integral part: */
  PREC_RES rinteg[ipn], /* Integrand */
           ipv[ipn];    /* Impact paramter where to integrate */

  /* Integrate for each of the planet's layer starting from the
     outermost until the closest layer  */
  for(i=0; i<=last; i++){
    ipv   [ipn1-i] = ip->v[i] * ip->fct;
    rinteg[ipn1-i] = exp(-tau[i]) * ipv[ipn1-i];
  }
  /* Add one more layer with 0. Only two to have a nice ending
    spline and not unnecessary values: */
  last += 1;
  if(last>ipn1) last = ipn1;
  for(; i<=last; i++){
    ipv   [ipn1-i] = ip->v[i] * ip->fct;
    rinteg[ipn1-i] = 0;
  }
  /* FINDME: there is no need to use a for-loop here, the first two lines
     are confusing also. This could have been written much simpler        */

  /* Increment last to represent the number of elements, check that we
     have enough: */
  last++;
  if(last<3)
    transiterror(TERR_CRITICAL,
                 "Condition failed, less than 3 items (only %i) for radial "
                 "integration.\n", last);

  /* Integrate along radius: */
#ifdef _USE_GSL
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp *spl       = gsl_interp_alloc(gsl_interp_cspline, last);
  gsl_interp_init(spl, ipv+ipn-last, rinteg+ipn-last, last);
  res = gsl_interp_eval_integ(spl, ipv+ipn-last, rinteg+ipn-last,
                               ipv[ipn-last], ipv[ipn1], acc);
  gsl_interp_free(spl);
  gsl_interp_accel_free (acc);
#else
# error computation of modulation() without GSL is not implemented
#endif

  /* TD: Add real unblocked area of the star, considering geometry */
  /* Substract the total area blocked by the planet. This is from the
     following:
     \begin{eqnarray}
     1-M = &1-\frac{\int_0^{r_p}\int e^{-\tau}r{\rm d}\theta {\rm d}r
                  \ +\ \int_{r_p}^{R_s}{\rm d}A} {\pi R_s^2}       \\
         = & -\frac{\int_0^{r_p}\int e^{-\tau}r{\rm d}\theta {\rm d}r
                  \ +\ Area_{p}} {\pi R_s^2}                       \\
         = & -\frac{2\int_0^{r_p} e^{-\tau}r{\rm d}r
                \ +\ r_p^2} {\pi R_s^2}   
     \end{eqnarray}                                                */

  res = ipv[ipn1]*ipv[ipn1] - 2.0*res;

  /* If the planet is going to be transparent with its maximum optical
     depth given by toomuch then: */
  if(sg->transpplanet)
    res -= exp(-maxtau) * ipv[ipn-last] * ipv[ipn-last];

  /* Normalize by the stellar radius: */
  res *= 1.0 / (srad*srad);

  return res;
}


/* \fcnfh
   Calculate the modulation at a given wavenumber, considering the planet
   as an opaque disc of radius r = r(tau=toomuch), for no-limb darkening
   nor planet emission.

   Return: transit's modulation, else
           -1 if toomuch was not reached                               */
static inline PREC_RES
modulationm1(PREC_RES *tau,        /* Optical depth array              */
             long last,            /* Index where tau = toomuch        */
             double toomuch,       /* Maximum optical depth calculated */
             prop_samp *ip,        /* Impact parameter                 */
             struct geometry *sg){ /* Geometry struct                  */
  long i,         /* Auxiliary for-loop index         */
       ini;       /* Impact parameter index           */
  double ipv[2];  /* Impact parameter around toomuch  */
  double muchrad; /* Radius where tau reaches toomuch */
  double srad = sg->starrad * sg->starradfct; /* Stellar radius */

  /* toomuch was not reached: */
  if(tau[last] < toomuch)
    return -1;

  /* Get impact parameter before and after tau=toomuch: */
  ini = ++last-2;
  if(ini<0) ini = 0;
  for(i=ini; i<last; i++)
    ipv[i-ini] = ip->v[i]*ip->fct;

  /* Find the planet radius through a linear interpolation: */
  muchrad = interp_line(tau+ini, ipv, toomuch);

  /* Return the modulation: */
  return muchrad * muchrad / (srad*srad);
}


/* \fcnfh
   Wrapper function to calculate the modulation in/out-of-transit ratio.

   Return: modulation                                      */
static PREC_RES
modulationperwn(PREC_RES *tau,       /* Optical depth                    */
                long last,           /* Radius index where tau = toomuch */
                double toomuch,      /* Maximum optical depth calculated */
                prop_samp *ip,       /* Impact parameter array           */
                struct geometry *sg, /* Geometry struct                  */
                int exprlevel){      /* Modlevel                         */
  switch(exprlevel){
  case 1:
    return modulation1(tau,  last, toomuch, ip, sg);
    break;
  case -1:
    return modulationm1(tau, last, toomuch, ip, sg);
    break;
  default:
    transiterror(TERR_CRITICAL,
                 "slantpath:: modulationperwn:: Level %i of detail "
                 "has not been implemented to compute modulation.\n",
                 exprlevel);
    return 0;
  }
}


const transit_ray_solution slantpath = {
       "Slant Path",     /* Name of the solution                     */
       "slantpath.c",    /* Source code file name                    */
       1,                /* Equispaced impact parameter requested?   */
       &totaltau,        /* Per impact parameter and per wavenumber 
                            value computation                        */
       &modulationperwn, /* Per wavenumber value computation         */
       1                 /* Number of levels of modulation detail as
                            it can be handled by the above fucntion  */
       };
