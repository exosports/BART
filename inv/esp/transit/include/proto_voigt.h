#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* voigt.c */
extern inline int voigtf P_((int nwn, float *wn, float wn0, double alphaL, double alphaD, float *vpro, double eps));
extern inline int voigtn P_((int m, int nwn, double dwn, double alphaL, double alphaD, float **vpro, double eps, int flags));

#undef P_
