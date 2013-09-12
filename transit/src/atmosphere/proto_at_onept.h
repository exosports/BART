#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/atmosphere/at_onept.c */
extern void askonenpt P_((struct onept *onept, struct atm_data *at, int rad));
extern void askonemn P_((struct onept *onept, prop_isof *isof, int n, int nf));
extern void sethcdef P_((struct transit *tr, struct atm_data *at, prop_samp *rads));

#undef P_
