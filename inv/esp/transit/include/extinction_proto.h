#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* extinction.c */
extern __inline__ int stat P_((__const char *__path, struct stat *__statbuf));
extern __inline__ int lstat P_((__const char *__path, struct stat *__statbuf));
extern __inline__ int fstat P_((int __fd, struct stat *__statbuf));
extern __inline__ int mknod P_((__const char *__path, __mode_t __mode, __dev_t __dev));
extern int extwn P_((struct transit *tr));
extern inline int newprofile P_((float **pr, int vf, int *lw, double dwn, float dop, float lor, float ta));

#undef P_
