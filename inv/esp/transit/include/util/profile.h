#ifndef _PROFILE_H
#define _PROFILE_H

#define PREC_VOIGT float
#define PREC_VOIGTP double

#include <P/Putil.h>
#include <transit.h>
#include <math.h>

#define VOIGT_QUICK 0x00001   //Quick integration.

/*
inline int voigtn(int m, int nwn, PREC_VOIGTP dwn, 
		  PREC_VOIGTP alphaL, PREC_VOIGTP alphaD,
		  PREC_VOIGT **vpro, PREC_VOIGTP eps, int flags);
*/

#include <voigt_proto.h>

#endif /* _PROFILE_H */
