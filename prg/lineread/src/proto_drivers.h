#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/drivers.c */
extern int db_drivers P_((struct hints *hint));
extern void drivers_free P_((struct hints *hint));
extern int find_alldrivers P_((struct hints *hint, unsigned short nfcn));
extern int setdriversnoutput P_((struct hints *hint));
extern int readwritepartition P_((unsigned short *acum));
extern int readwritetransition P_((unsigned short *acum, double ini, double fin, double del));

#undef P_
