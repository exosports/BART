#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/atmosphere/at_file.c */
extern int getmnfromfile P_((FILE *fp, struct atm_data *at, struct transit *tr, PREC_ZREC *f_remainder));
extern int readatmfile P_((FILE *fp, struct transit *tr, struct atm_data *at, prop_samp *rads, int nrad, PREC_ZREC *f_remainder));
extern void storename P_((struct atm_data *at, char *line));
extern void getmass P_((struct atm_data *at, struct molecules *mol));
#undef P_
