#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/drivers.c */
extern int db_drivers P_((struct hints *hint));
extern int find_alldrivers P_((struct hints *hint));
extern FILE *setdriversnoutput P_((struct hints *hint));
extern int readwritetransition P_((unsigned short *niso, FILE *fpout, double ini, double fin, double del));
extern int readwritepartition P_((unsigned short *niso, FILE *fpout));

#undef P_
