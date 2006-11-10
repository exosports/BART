#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/dbread_pands.c */
extern driver_func *initdb_pands P_((void));

#undef P_
