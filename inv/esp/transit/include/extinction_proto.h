#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* extinction.c */
extern int extwn P_((struct transit *tr));
extern inline int newprofile P_((float **pr, int vf, double dwn, float dop, float lor, float ta));

#undef P_
