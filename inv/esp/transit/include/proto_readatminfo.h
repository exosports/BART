#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* readatminfo.c */
extern int getatm P_((struct transit *tr));
extern double checkaddmm P_((double *mm, long r, prop_isov *isov, prop_isof *isof, int n, _Bool mass, enum isodo *isodo));
extern int getmnfromfile P_((FILE *fp, struct atm_data *at, struct transit *tr, int nmb));
extern int readatmfile P_((FILE *fp, struct transit *tr, struct atm_data *at, prop_samp *rads, int nrad));
extern void storename P_((struct atm_data *at, char *line));
extern void sethcdef P_((struct transit *tr, struct atm_data *at, prop_samp *rads));

#undef P_
