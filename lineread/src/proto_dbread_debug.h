#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/dbread_debug.c */
extern driver_func *initdb_debug P_((void));

#undef P_
