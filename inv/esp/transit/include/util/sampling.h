#ifndef _SAMPLING_H
#define _SAMPLING_H

#ifdef TRANSIT
#define SAMP_BITS   TRU_SAMPBITS
#define SAMP_LINEAR TRU_SAMPLIN
#define SAMP_SPLINE TRU_SAMPSPL
#else
#define SAMP_BITS   0x0000000f
#define SAMP_LINEAR 0x00000001
#define SAMP_SPLINE 0x00000002
#endif

#define resamplex(fl,nrefx,refx,ndatx,datx)      \
          resample(1,fl,nrefx,refx,ndatx,datx,0)
#define resampley(fl,ny,...);                    \
          resample(0,fl,0,NULL,0,NULL,ny,__VA_ARGS__)

/*
inline void natcubsplinecoef(long n, double *x, double *y, double *h, double *D);
inline int resample(const int getingx, const long flag, long nrefx, double *refx, long ndatx, double *datx, long ny, ...);
*/

#include <sampling_proto.h>


#endif
