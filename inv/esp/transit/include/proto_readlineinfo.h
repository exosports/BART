#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* readlineinfo.c */
extern int readtwii_bin P_((FILE *fp, struct transit *tr, struct lineinfo *li));
extern int readtwii_ascii P_((FILE *fp, struct transit *tr, struct lineinfo *li));
extern int getinifinasctwii P_((double *ini, double *fin, FILE *fp, char *file));
extern int checkrange P_((struct transit *tr, struct lineinfo *li));
extern int readinfo_twii P_((struct transit *tr, struct lineinfo *li));
extern int readdatarng P_((struct transit *tr, struct lineinfo *li));
extern int readlineinfo P_((struct transit *transit));
extern int free_lineinfotrans P_((struct line_transition *lt, struct lineinfo *li, long *pi));

#undef P_
