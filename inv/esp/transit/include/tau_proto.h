#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* tau.c */
extern int tau P_((struct transit *tau));

#undef P_
