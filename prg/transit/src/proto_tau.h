#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/tau.c */
extern int tau P_((struct transit *tr));
extern void outdebtauex P_((char *name, double **e, prop_samp *ip, double **t, long rn, long w));
extern void outdebex P_((char *name, double **e, double *r, long rn, long wi, long wf));
extern void outdebtau P_((char *name, prop_samp *ip, double **t, long wi, long wf));
extern void printtoomuch P_((char *file, struct optdepth *tau, prop_samp *wn, prop_samp *rad));
extern void printtau P_((struct transit *tr));
extern int freemem_tau P_((struct optdepth *tau, long *pi));
extern int detailout P_((prop_samp *wn, prop_samp *rad, struct detailfld *det, double **arr, short flag));

#undef P_
