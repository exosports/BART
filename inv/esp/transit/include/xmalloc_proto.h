#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* xmalloc.c */
extern void *xmalloc P_((size_t n));
extern void *xcalloc P_((size_t n, size_t s));
extern void *xrealloc P_((void *p, size_t n));
extern char *xstrdup P_((char *str));

#undef P_
