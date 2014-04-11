#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/transitstd.c */
extern inline void transitdot P_((int thislevel, int verblevel, ...));
extern int transiterror_fcn P_((int flags, const char *file, const long line, const char *str, ...));
extern int vtransiterror_fcn P_((int flags, const char *file, const long line, const char *str, va_list ap));
extern int fileexistopen P_((char *in, FILE **fp));
extern FILE *verbfileopen P_((char *in, char *desc));
extern void transitcheckcalled P_((const long pi, const char *fcn, const int n, ...));
extern void error P_((int exitstatus, int something, const char *fmt, ...));
extern void free_isov P_((prop_isov *isov));
extern void freemem_molecules P_((struct molecules *mol, long *pi));
extern void free_isof P_((prop_isof *isof));
extern void free_mol P_((prop_mol *molec));
extern void free_db P_((prop_db *db));
extern void free_dbnoext P_((prop_dbnoext *db));
extern void free_samp P_((prop_samp *samp));
extern void free_atm P_((prop_atm *atm));
extern void savestr P_((FILE *out, char *str));
extern int reststr P_((FILE *in, char **str));
extern void linetoolong P_((int max, char *file, int line));
extern double timecheck P_((int verblevel, long iter, long index, char *str, struct timeval tv1, double t0));

#undef P_
