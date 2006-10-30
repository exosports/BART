#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/idxrefraction.c */
extern int idxrefrac P_((struct transit *tr));
extern int freemem_idexrefrac P_((struct idxref *ir, long *pi));
extern int restidxref P_((FILE *in, long nrad, struct idxref *ir));
extern void saveidxref P_((FILE *out, long nrad, struct idxref *ir));

#undef P_
