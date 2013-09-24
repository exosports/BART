#ifndef _SAMPLING_H
#define _SAMPLING_H

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#define SAMP_BITS     0x0000000f
#define SAMP_LINEAR   0x00000001
#define SAMP_SPLINE   0x00000002

#define INTERP_LINEAR 0x0000001


#define resamplex(fl,nrefx,refx,ndatx,datx)      \
          resample(1,fl,nrefx,refx,ndatx,datx,0)
#define resampley(fl,ny,...)                     \
          resample(0,fl,0,NULL,0,NULL,ny,__VA_ARGS__)
#define resample_free()                          \
          resample(2,0,0,NULL,0,NULL,0)


//function definition (proto_sampling.h)
#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* sampling.c */
extern inline int resample P_((const short getingx, const long flags, long nrefx, double *refx, long noutx, double *outx, long ny, ...));
extern int lineinterpol P_((int ndat, double *x, double *y, int n, long *indx, float *t, double *yout, double *dbgout));
extern void natcubspline P_((int ndat, double *x, double *y, int n, long *indx, float *t, double *yout, double *dout));
extern inline void natcubsplinecoef P_((long n, double *x, double *y, double *h, double *D));
extern inline double interp P_((double refx, double *x, double *y, long n, int intkind));
extern inline double lineinterp P_((double refx, double *x, double *y, long n));

#undef P_


#endif /* _SAMPLING_H */
