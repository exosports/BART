#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* iomisc.c */
extern int ncharchg P_((char *str, char car, char chg));
extern char *readstr_sp P_((char *line, char **next, char fspace));
extern inline char fgetupto P_((char *line, int max, FILE *fp, void (*errfcn)(int, char *, int), char *name, int other));
extern int getad P_((int n, char sep, char *str, double **array));
extern int getnd P_((int n, char sep, char *str, ...));
extern int getnl P_((int n, char sep, char *str, ...));
extern void fprintpad P_((FILE *fp, int indent, char *fmt, ...));
extern double readds P_((FILE *fp, char *c, char *string, int maxstring));
extern double getds P_((char *in, char *c, char *string, int maxstring));
extern long readl P_((FILE *fp, char *c));
extern char *linepad P_((char *out, int nc, char *in));
extern double askforposd P_((char *fmt, ...));
extern long askforposl P_((char *fmt, ...));
extern char *fgets_alloc P_((FILE *fp, int *max));

#undef P_
