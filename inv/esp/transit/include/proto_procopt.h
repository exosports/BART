#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* procopt.c */
extern int getprocopt P_((int argc, char **argv, struct optdocs *opts, struct optcfg *cfg, int *longidxp));
extern int getopt_long_files P_((int argc, char **argv, char *shortopts, struct option *getopts, int *longidxp, char *paramfilelist));
extern int getopt_long_files_old P_((int argc, char **argv, char *shortopts, struct option *getopts, int *longidxp, char *paramfilelist));
extern void prochelp P_((int status));
extern void getprocopt_free P_((void));

#undef P_
