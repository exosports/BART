#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* sampling.c */
extern inline int resample P_((const int getingx, const long flags, long nrefx, double *refx, long noutx, double *outx, long ny, ...));
extern int lineinterpol P_((int ndat, double *x, double *y, int n, long *indx, float *t, double *yout, double *dbgout));
extern void natcubspline P_((int ndat, double *x, double *y, int n, long *indx, float *t, double *yout, double *dout));
extern inline void natcubsplinecoef P_((long n, double *x, double *y, double *h, double *D));
extern inline double interp P_((double refx, double *x, double *y, long n, int intkind));
extern double lineinterp P_((double refx, double *x, double *y, long n));

#undef P_
