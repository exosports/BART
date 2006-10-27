#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/messagep.c */
extern inline void msgpdot P_((int thislevel, int verblevel, ...));
extern int mperror_fcn P_((int flags, const char *file, const long line, const char *str, ...));
extern int vmperror_fcn P_((int flags, const char *file, const long line, const char *str, va_list ap));
extern int fileexistopen P_((char *in, FILE **fp));
extern FILE *verbfileopen P_((char *in, char *desc));
extern void linetoolong P_((int max, char *file, int line));

#undef P_
