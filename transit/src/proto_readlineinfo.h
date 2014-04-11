#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/readlineinfo.c */
extern int readtli_bin P_((FILE *fp, struct transit *tr, struct lineinfo *li));
extern int readtli_ascii P_((FILE *fp, struct transit *tr, struct lineinfo *li));
extern int getinifinasctli P_((double *ini, double *fin, FILE *fp, char *file));
extern int setimol P_((struct transit *tr));
extern int checkrange P_((struct transit *tr, struct lineinfo *li));
extern int readinfo_tli P_((struct transit *tr, struct lineinfo *li));
extern int readdatarng P_((struct transit *tr, struct lineinfo *li));
extern int readlineinfo P_((struct transit *tr));
extern int freemem_isotopes P_((struct isotopes *iso, long *pi));
extern int freemem_lineinfotrans P_((struct lineinfo *li, long *pi));
extern void saveline P_((FILE *fp, struct lineinfo *li));

#undef P_
