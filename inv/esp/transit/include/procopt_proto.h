#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* procopt.c */
extern int getprocopt P_((int argc, char **argv, char *paramfile, struct optdocs *opts, struct optcfg *cfg, int *longidxp));
extern void prochelp P_((int status, struct optdocs *opts, struct optcfg *cfg));

#undef P_
