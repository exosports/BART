#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* readlineinfo.c */
extern __inline__ int stat P_((__const char *__path, struct stat *__statbuf));
extern __inline__ int lstat P_((__const char *__path, struct stat *__statbuf));
extern __inline__ int fstat P_((int __fd, struct stat *__statbuf));
extern __inline__ int mknod P_((__const char *__path, __mode_t __mode, __dev_t __dev));
extern int readtwii_bin P_((FILE *fp, struct transit *tr, struct lineinfo *li));
extern int readtwii_ascii P_((FILE *fp, struct transit *tr, struct lineinfo *li));
extern int getinifinasctwii P_((double *ini, double *fin, FILE *fp, char *file));
extern int checkrange P_((struct transit *tr, struct lineinfo *li));
extern int readinfo_twii P_((struct transit *tr, struct lineinfo *li));
extern int readdatarng P_((struct transit *tr, struct lineinfo *li));
extern int readlineinfo P_((struct transit *transit));

#undef P_
