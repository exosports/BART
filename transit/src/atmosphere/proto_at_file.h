#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/atmosphere/at_file.c */
extern int getmnfromfile P_((FILE *fp, struct atm_data *at, struct transit *tr, int nmb));
extern int readatmfile P_((FILE *fp, struct transit *tr, struct atm_data *at, prop_samp *rads, int nrad));
extern void storename P_((struct atm_data *at, char *line));

#undef P_
