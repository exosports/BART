#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/dbread_text.c */
extern void linetoolong_text P_((int max, char *file, int line));
extern driver_func *initdb_text P_((void));

#undef P_
