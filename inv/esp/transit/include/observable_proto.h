#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* observable.c */
extern int modulation P_((struct transit *tr));

#undef P_
