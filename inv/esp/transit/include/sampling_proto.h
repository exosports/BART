#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* sampling.c */
extern __inline__ int stat P_((__const char *__path, struct stat *__statbuf));
extern __inline__ int lstat P_((__const char *__path, struct stat *__statbuf));
extern __inline__ int fstat P_((int __fd, struct stat *__statbuf));
extern __inline__ int mknod P_((__const char *__path, __mode_t __mode, __dev_t __dev));
extern inline int resample P_((const int getingx, const long flags, long nrefx, double *refx, long noutx, double *outx, long ny, ...));
extern int lineinterpol P_((int ndat, double *x, double *y, int n, long *indx, float *t, double *yout, double *dbgout));
extern void natcubspline P_((int ndat, double *x, double *y, int n, long *indx, float *t, double *yout, double *dout));
extern inline void natcubsplinecoef P_((long n, double *x, double *y, double *h, double *D));

#undef P_
