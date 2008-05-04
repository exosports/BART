#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/argum.c */
extern int argum P_((int argc, char **argv, struct hints *hint));
extern void hints_free P_((struct hints *hint));

#undef P_
