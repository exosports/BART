#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* idxrefraction.c */
extern int idxrefrac P_((struct transit *tr));

#undef P_
