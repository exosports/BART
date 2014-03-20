/*
 * voigt.c - Functions to compute the Voigt profile
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

/*
  voigt.c: Functions to return a Voigt profile according to the numerical
  approximation described by Pierlusi et al. in
  J. Quant. Spectrosc. Radiat. Transfer., Vol 18 pp.555

   Normalized line shape is given by:
    \Psi(x,y) = \frac{y}{\pi}
                \int_{-\infty}^{\infty} \frac{\exp(-t^2)}{y^2+(x-t)^2} {\rm d}t
   However, all the functions to be used will be approximations in
   different regions to the function:
     Psi(x,y) = Re[w(z=x+iy)]
              = Re[exp(-z^2)(1-erf(-iz))]

   (c) Patricio Rojo 2003                                                  */

  /* TD: Set iteration limit from fitted function and not convergence
     check!  */
  /* TD: check that function using voigt array hanndles correctly if the
     item happens to be exactly in between two bins: it should be 0.5
     the contribution from its of its boundary bin in the first and
     latter \emph{center shift position}*/
//\omitfh
#include <pu/profile.h>

#define SQRTLN2 0.83255461115769775635
#define TWOOSQRTPI 1.12837916709551257389
#define SQRTLN2PI  0.46971863934982566689

#define A1  0.46131350
#define A2  0.19016350
#define A3  0.09999216
#define A4  1.78449270
#define A5  0.002883894
#define A6  5.52534370

#define B1  0.51242424
#define B2  0.27525510
#define B3  0.05176536
#define B4  2.72474500

#define MAXCONV 61
/* ferf[n]= 1/(n!(2n+1)) */
static double ferf[MAXCONV] = {
  1.000000000000000000000,    //0
  0.333333333333333333333,    //1
  0.100000000000000000000,    //2
  2.38095238095238095238e-2,  //3
  4.62962962962962962963e-3,  //4
  7.57575757575757575758e-4,  //5
  1.06837606837606837607e-4,  //6
  1.32275132275132275132e-5,  //7
  1.45891690009337068161e-6,  //8
  1.45038522231504687645e-7,  //9
  1.31225329638028050726e-8,  //10
  1.08922210371485733805e-9,  //11
  8.35070279514723959168e-11, //12
  5.94779401363763503681e-12, //13
  3.95542951645852576340e-13, //14
  2.46682701026445692771e-14, //15
  1.44832646435981372650e-15, //16
  8.03273501241577360914e-17, //17
  4.22140728880708823303e-18, //18
  2.10785519144213582486e-19, //19
  1.00251649349077191670e-20, //20
  4.55184675892820028624e-22, //21
  1.97706475387790517483e-23, //22
  8.23014929921422135684e-25, //23
  3.28926034917575173275e-26, //24
  1.26410789889891635220e-27, //25
  4.67848351551848577373e-29, //26
  1.66976179341737202699e-30, //27
  5.75419164398217177220e-32, //28
  1.91694286210978253077e-33, //29
  6.18030758822279613746e-35, //30
  1.93035720881510785656e-36, //31
  5.84675500746883629630e-38, //32
  1.71885606280178362397e-39, //33
  4.90892396452342296700e-41, //34
  1.36304126177913957635e-42, //35
  3.68249351546114573519e-44, //36
  9.68728023887076175384e-46, //37
  2.48306909745491159104e-47, //38
  6.20565791963739670594e-49, //39
  1.51310794954121709805e-50, //40
  3.60157930981012591661e-52, //41
  8.37341968387228154283e-54, //42
  1.90254122728987952724e-55, //43
  4.22678975419355257584e-57, //44
  9.18642950239868569596e-59, //45
  1.95410258232417110410e-60, //46
  4.07013527785325672298e-62, //47
  8.30461450592911058168e-64, //48
  1.66058051345108993284e-65, //49
  3.25539546201302778914e-67, //50
  6.25918411694871134025e-69, //51
  1.18076183891157008800e-70, //52
  2.18621042295388572103e-72, //53
  3.97425272266506578576e-74, //54
  7.09571739181805357327e-76, //55
  1.24466597738907071213e-77, //56
  2.14564844309633852739e-79, //57
  3.63615636540051474579e-81, //58
  6.05939744697137480783e-83, //59
  9.93207019544894768776e-85
};
int _voigt_maxelements=99999;
int _voigt_computeeach=10;

#define NFCN(x,y) (x<1?15:(int)(6.842*x+8.0))

static inline int
meanintegSimp(PREC_VOIGT *in,
              PREC_VOIGT *out,
              int ni,
              int no,
              int ipo,
              PREC_VOIGT d);


static inline int
meanintegTrap(PREC_VOIGT *in,
              PREC_VOIGT *out,
              int ni,
              int no,
              int ipo,
              PREC_VOIGT d);


static inline void
voigtxy(double x,
        double y,
        PREC_VOIGT *res,
        PREC_VOIGT eps,
        PREC_VOIGTP alphaD){
  int i, n;
  long double or, nr, oi, ni, ar, ai; 
  const long double x2y2 = x*x-y*y;
  const long double xy2 = 2*x*y;
  const long double cosxy = cosl(xy2);
  const long double sinxy = sinl(xy2);

#ifdef _VOIGT_USE_CONVERG
  if(eps>0) n = MAXCONV;
  else 
#endif /* Precision convergence */
    n = NFCN(x, y) + 1;

  if(x<3 && y<1.8){ /* Region I */
    ar = or =  y;
    ai = oi = -x;
    i = 1;
    while(1){
      ni = or*xy2  + oi*x2y2;
      nr = or*x2y2 - oi*xy2 ;
#ifdef _VOIGT_USE_CONVERG
      if(n == MAXCONV){
        if(i > n){
          fprintf(stderr, "%s:: No convergence after %i iterations in "
                          "VOIGTXY() calculations (voigt.c file)\n",
                          __FILE__, MAXCONV);
          break;
        }
        if(fabs(ferf[i]*(cosxy*nr + sinxy*ni))<eps)
          break;
      }
      else 
#endif /* Precision convergence */
        if(i > n) break;
      ai += ni*ferf[i];
      ar += nr*ferf[i];

      oi = ni;
      or = nr;
      i++;
    }
    (*(res)) = SQRTLN2PI/alphaD*exp(-x2y2)*
               (cosxy*(1-ar*TWOOSQRTPI)-sinxy*ai*TWOOSQRTPI);
  }
  else if(x<5 && y<5){    /* Region II  */
    ar = xy2*xy2;
    nr = xy2*x;
    ni = x2y2-A2;
    ai = x2y2-A4;
    oi = x2y2-A6;
    (*(res)) = SQRTLN2PI/alphaD*(A1*((nr-ni*y)/(ni*ni+ar)) +
                                 A3*((nr-ai*y)/(ai*ai+ar)) +
                                 A5*((nr-oi*y)/(oi*oi+ar)) );
  }
  else{                   /* Region III */
    ar = xy2*xy2;
    nr = xy2*x;
    ni = x2y2-B2;
    ai = x2y2-B4;
    (*(res)) = SQRTLN2PI/alphaD*(B1*((nr-ni*y)/(ni*ni+ar)) +
                                 B3*((nr-ai*y)/(ai*ai+ar)) );
  }
} __attribute__((always_inline))


//\fcnfh
//Computes Voigt Profile 
inline int
voigtf(int nwn,            /* Number of wavenumber bins where the function is
                              to be computed. i.e., the size of wn         */
       PREC_VOIGT *wn,     /* Wavenumber array bins where to evaluate the
                              Voigt function                               */
       PREC_VOIGT wn0,     /* Profile wavenumber center                    */
       PREC_VOIGTP alphaL, /* Lorentz width                                */
       PREC_VOIGTP alphaD, /* Doppler width                                */
       PREC_VOIGT *vpro,   /* Profile values array                         */
       PREC_VOIGTP eps){   /* Precision to which the profile is calculated.
                              If negative, a fixed number of iterations is
                              calculated. */
  double y, x;
  int i;

  y = SQRTLN2*alphaL/alphaD;

  for(i=0; i<nwn; i++){
    x = SQRTLN2*fabs(wn[i]-wn0)/alphaD;
    voigtxy(x, y, vpro+i, eps, alphaD);
  }
  return 1;
}


/*\fcnfh
  Compute Voigt profile on equispaced grid.

  Return: 1 on success                           */
inline int 
voigtn(int m,              /* Number of fine resolution bins (deviation
                              of line center from bin center)             */
       int nwn,            /* Number of wavenumber bins of the profile    */
       PREC_VOIGTP dwn,    /* Wavenumber half-width sample                */
       PREC_VOIGTP alphaL, /* Lorentz width                               */
       PREC_VOIGTP alphaD, /* Doppler width                               */
       PREC_VOIGT **vpro,  /* Array (m by nwn) where to store the profile */
       PREC_VOIGTP eps,    /* Precision to which the profile is to be
                              calculated. If negative, a fixed number of
                              iterations are performed.                     */
       int flags){         /* Miscellaneous flags, so far there is support
                              for: 'VOIGT_QUICK' that performs a quick
                              integration, i.e., the height multiplied
                              by the bin width                              */

  /* The calculation is done filling the array from the centerpoint
     outwards.  The wavenumber value of each bin is considered to be the
     lower value of the bin (not the center as usual).                    */

  double y,      /* Normalized width:    y = sqrt(ln 2) * alphaL/alphaD   */
         x,      /* Normalized position: x = sqrt(ln 2) * |nu-nu0|/alphaD */
         ddwn,   /* Spacing between wavenumbers                           */
         dcshft; /* Centershift spacing                                   */
  int i, j;

  /* The following variables are used to make sure that we are
     obtaining a fine integrated value for each bin:
     dint is the maximum spacing thus implied by the value of nint */
  int nint;    /* Minimum number of points calculated per Doppler width */
  double dint, /* Maximum spacing (?) */
         shft; /* comment me          */
  PREC_VOIGT *aint; /* Auxiliary array to put the fine-spaced profile  */

  /* Initialization:  */
  y      = SQRTLN2 * alphaL / alphaD;
  ddwn   = 2.0 * dwn / (nwn - 1);
  dcshft = ddwn / m;
  nint   = 50;
  dint   = alphaD / (nint - 1);

  /* If only a quick integration is requested, no need to have more than
     the requested bin.  Note that this may cause a binning smaller than
     the minimum. Or the same will happen if the requested spacing is
     smaller than integration minimum.
     If so, we need to allocate aint only to one item bigger
     than what was given in order to calculate the last point needed
     for integration. Also is important to set the \vr{flags} so that the
     program knows later that it has to do 2 point simple trapezoid
     integration. Set the fine-binned values to the given
     ones. Temporarily enlarge by one the array in order to be able to
     integrate the last bin. */
  if( ddwn<dint || (flags&VOIGT_QUICK) ){
    dint = ddwn;
    nint = nwn + 1;
  }

  /* If integration is required then calculate binning so that for each
     wavenumber bin there is an odd number of fine-integration bins (in
     order to be able to make simpson integration.
     Note that \vr{nint} includes a set of (int)(\vr{ddwn}/\vr{dint})
     fine-points for every wavenumber point and one extra point for the
     one required to close the last bin.  */
  else{
    nint = (int)(ddwn / dint) + 1;
    if(nint & 1)
      nint++;
    nint = nwn * nint + 1;
    dint = 2.0 * dwn / (nint - 1);
  }

  /* Initialize aint array: */
  if((aint=(PREC_VOIGT *)calloc(nint, sizeof(PREC_VOIGT)))==NULL){
    fprintf(stderr, "\nUnable to allocate memory in voigtxy for %i double "
                    "elements (%.3g MB).\n", nint,
                    (sizeof(PREC_VOIGT)*nint)/1024.0/1024.0);
    fprintf(stderr, "nwn: %i\n", nwn);
    exit(EXIT_FAILURE);
  }

  /* Work sequentially (index \vr{j}) for each centershift: */
  for(j=0; j<m; j++){
    /* So that the array with the central value of m have the most
       centered array: */
    shft = dwn + (j-m/2)*dcshft;

    /* Even though Voigt is symmetric, a symmetric filling of the answer
       array is not possible because the center of
       the profile won't always coincide with the center of the array;
       only when \vr{m}=0. Then, for each element of the array.          */
    for(i=0; i<nint; i++){
      x = SQRTLN2*fabs(dint*i-shft)/alphaD;
      voigtxy(x, y, aint+i, eps, alphaD);
    }

    /* Integration: */
    /* If user wants a quick integration, return the value at the
       beginning of each bin: */
    if(flags & VOIGT_QUICK){ 
      for(i=0; i<nwn; i++)
        vpro[j][i] = aint[i];
    }

    /* If a slow integration is preferred, then calculate the number of
       fine-point per regular points, note that is one more than the ratio. */
    else{
      /* \linelabel{aux.ddwn} */
      ddwn = (float)(nint-1)/nwn;
      i = (int)ddwn+1;
      if(ddwn+1 != i){
        fprintf(stderr, "%s:: There are not an integer number(%f!=%i) of "
                        "fine-bins(%i) for each bin(%i).\n", 
                        __FILE__, ddwn+1, i, nint-1, nwn);
        exit(EXIT_FAILURE);
      }

      /* Is this number an odd number?
         If so, use Simpson's integration, otherwise use trapezoidal: */
      if(i & 1)
        meanintegSimp(aint, vpro[j], nint, nwn, i, dint);
      else
        meanintegTrap(aint, vpro[j], nint, nwn, i, dint);
    }

  }
  free(aint);

  /* On success return 1.  So far, it is the only answer this function is
     giving. */
  return 1;
}


/*\fcnfh
  Simpson integration
  Return: 1 on success     */
static inline int
meanintegSimp(PREC_VOIGT *in,  /* Input array                 */
              PREC_VOIGT *out, /* Output array                */
              int ni,          /* Number of input elements    */
              int no,          /* Number of output elements   */
              int ipo,         /* Number of inputs per output */
              PREC_VOIGT d){   /* bin width                   */

  //Indexes for the in array. To speed up, \vr{ipo} will be the
  //last index in the sub-arrays instead of the number of elements
  int i;
  ipo--;

  //For each result...
  for(;no;no--){
    //do the integration
    *out=0;
    for(i=1;i<ipo;i+=2)
      *out+=in[i];
    *out*=2;
    for(i=2;i<ipo;i+=2)
      *out+=in[i];
    *out*=2;
    *out+=*in+in[ipo];
    *out/=(ipo*3.0);
    //advance the input array to the next bin and advance in one the
    //return array
    in+=ipo;
    out++;
  }
  //returns 1 if success... It always does this
  return 1;
}

/*\fcnfh
  Trapezoid integration
  @returns 1 if success.
*/
static inline int
meanintegTrap(PREC_VOIGT *in,        /* Input Array */
              PREC_VOIGT *out,        /* Output Array */
              int ni,                /* Number of input elements */         
              int no,                /* Number of output elements */         
              int ipo,                /* Number of inputs per output */
              PREC_VOIGT d)        /* bin width */
{
  //Index for the in array. To speed up, \vr{ipo} will be the
  //last index in the sub-arrays instead of the number of elements
  int i;
  ipo--;

  //For each result...
  for(;no;no--){
    //do the integration
    *out=0;
    for(i=1;i<ipo;i++)
      *out+=in[i];
    *out=(*out+(in[0]+in[ipo])/2.0)/(double)ipo;
    //advance the input array to the next bin, and advance in one the
    //return array
    in+=ipo;
    out++;
  }
  //returns 1 if success... It always does this
  return 1;
}

#undef SQRTLN2
#undef TWOOSQRTPI
#undef MAXCONV
#undef NFCN
#undef VOIGTXY

#undef A1
#undef A2
#undef A3
#undef A4
#undef A5
#undef A6

#undef B1
#undef B2
#undef B3
#undef B4

