#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* transitstd.c */
extern __inline__ int stat P_((__const char *__path, struct stat *__statbuf));
extern __inline__ int lstat P_((__const char *__path, struct stat *__statbuf));
extern __inline__ int fstat P_((int __fd, struct stat *__statbuf));
extern __inline__ int mknod P_((__const char *__path, __mode_t __mode, __dev_t __dev));
extern inline void transitdot P_((int thislevel, int verblevel));
extern int transiterror P_((int flags, char *str, ...));
extern int fileexistopen P_((char *in, FILE **fp));
extern int verbfileopen P_((char *in, FILE **fp, char *desc));
extern void transitcheckcalled P_((const long pi, const char *fcn, const int n, ...));

#undef P_
