#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* readatminfo.c */
extern __inline__ int stat P_((__const char *__path, struct stat *__statbuf));
extern __inline__ int lstat P_((__const char *__path, struct stat *__statbuf));
extern __inline__ int fstat P_((int __fd, struct stat *__statbuf));
extern __inline__ int mknod P_((__const char *__path, __mode_t __mode, __dev_t __dev));
extern int getatm P_((struct transit *tr));
extern double checkaddmm P_((double *mm, long r, prop_isov *isov, prop_isof *isof, int n, _Bool mass, enum isodo *isodo));
extern double askforposd P_((char *fmt, ...));
extern long askforposl P_((char *fmt, ...));
extern int getmnfromfile P_((FILE *fp, struct atm_data *at, struct transit *tr, int nmb));
extern int readatmfile P_((FILE *fp, struct transit *tr, struct atm_data *at, prop_samp *rads, int nrad));
extern void storename P_((struct atm_data *at, char *line));

#undef P_
