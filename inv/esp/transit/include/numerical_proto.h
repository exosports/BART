#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* numerical.c */
extern inline int binsearchie P_((double *arr, long i, long f, double val));
extern inline int binsearchei P_((double *arr, long i, long f, double val));
extern inline int binsearch P_((double *arr, long i, long f, double val));

#undef P_
