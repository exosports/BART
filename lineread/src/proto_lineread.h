#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/lineread.c */
extern int main P_((int argc, char **argv));
extern void lineread_free P_((void));

#undef P_
