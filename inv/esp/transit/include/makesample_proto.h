#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* makesample.c */
extern int makesample P_((prop_samp *samp, prop_samp *hint, prop_samp *ref, const long fl, const int bitsshift, const float margini, const float marginf));
extern int makewavsample P_((struct transit *tr));
extern int makewnsample P_((struct transit *tr));
extern int makeipsample P_((struct transit *tr));
extern int makeradsample P_((struct transit *tr));

#undef P_
