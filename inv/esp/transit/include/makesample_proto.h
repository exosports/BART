#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* makesample.c */
extern __inline__ int stat P_((__const char *__path, struct stat *__statbuf));
extern __inline__ int lstat P_((__const char *__path, struct stat *__statbuf));
extern __inline__ int fstat P_((int __fd, struct stat *__statbuf));
extern __inline__ int mknod P_((__const char *__path, __mode_t __mode, __dev_t __dev));
extern int makesample P_((prop_samp *samp, prop_samp *hint, prop_samp *ref, const long fl, const int bitsshift, const float margini, const float marginf));
extern int makewavsample P_((struct transit *tr));
extern int makewnsample P_((struct transit *tr));
extern int makeradsample P_((struct transit *tr));

#undef P_
