#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* voigt.c */
extern __inline__ int stat P_((__const char *__path, struct stat *__statbuf));
extern __inline__ int lstat P_((__const char *__path, struct stat *__statbuf));
extern __inline__ int fstat P_((int __fd, struct stat *__statbuf));
extern __inline__ int mknod P_((__const char *__path, __mode_t __mode, __dev_t __dev));
extern inline int voigtf P_((int nwn, float *wn, float wn0, double alphaL, double alphaD, float *vpro, double eps));
extern inline int voigtn P_((int m, int nwn, double dwn, double alphaL, double alphaD, float **vpro, double eps, int flags));

#undef P_
