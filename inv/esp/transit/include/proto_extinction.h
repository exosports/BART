#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* extinction.c */
extern int extwn P_((struct transit *tr));
extern inline int newprofile P_((float **pr, int vf, double dwn, float dop, float lor, float ta));
extern void printone P_((struct transit *tr));
extern int freemem_extinction P_((struct extinction *ex, long *pi));
extern int restextinct P_((FILE *in, long nrad, short niso, long nwav, struct extinction *ex));
extern void saveextinct P_((FILE *out, long nrad, short niso, long nwn, struct extinction *ex));
extern int savefile_extwn P_((struct transit *tr));

#undef P_
