#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* procopt.c */
extern int getprocopt P_((int argc, char **argv, struct optdocs *opts, struct optcfg *cfg, int *longidxp));
extern int getopt_long_def P_((int argc, char **argv, char *shortopts, struct option *getopts, int *longidxp, char *paramfilelist));
extern void prochelp P_((int status, struct optdocs *opts, struct optcfg *cfg));

#undef P_
