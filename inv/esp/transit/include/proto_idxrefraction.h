#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* idxrefraction.c */
extern int idxrefrac P_((struct transit *tr));
extern int freemem_idexrefrac P_((struct idxref *ir, long *pi));

#undef P_
