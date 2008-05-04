#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/geometry.c */
extern int setgeomhint P_((struct transit *tr));
extern int setgeom P_((struct geometry *sg, double time, long *flags));
extern inline double starvariation P_((double x, double y, double radius));

#undef P_
