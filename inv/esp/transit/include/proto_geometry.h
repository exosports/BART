#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* geometry.c */
extern int setgeom P_((struct transit *tr));

#undef P_
