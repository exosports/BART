#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* tau.c */
extern int tau P_((struct transit *tr));
extern void printtau P_((struct transit *tr));

#undef P_
