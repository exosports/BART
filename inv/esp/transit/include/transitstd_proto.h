#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* transitstd.c */
extern inline void transitdot P_((int thislevel, int verblevel));
extern int transiterror P_((int flags, char *str, ...));
extern int fileexistopen P_((char *in, FILE **fp));
extern int verbfileopen P_((char *in, FILE **fp, char *desc));
extern void transitcheckcalled P_((const long pi, const char *fcn, const int n, ...));

#undef P_
