#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* transit.c */
extern int main P_((int argc, char **argv));
extern void printone P_((struct transit *tr));
extern int processparameters P_((int argc, char **argv, struct transithint *hints));

#undef P_
