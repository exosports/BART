#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* tau.c */
extern int tau P_((struct transit *tr));
extern void printtoomuch P_((char *file, struct optdepth *tau, prop_samp *wn, prop_samp *rad));
extern void printtau P_((struct transit *tr));
extern int freemem_tau P_((struct optdepth *tau, long *pi));

#undef P_
