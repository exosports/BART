#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* numerical.c */
extern inline int binsearchie P_((double *arr, long i, long f, double val));
extern inline int binsearchei P_((double *arr, long i, long f, double val));
extern inline int binsearch P_((double *arr, long i, long f, double val));
extern inline double integ_trasim P_((double dx, double *y, long n));
extern inline double interp_parab P_((double *x, double *y, double xr));
extern inline double interp_line P_((double *x, double *y, double xr));

#undef P_
