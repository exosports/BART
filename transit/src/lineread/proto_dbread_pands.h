#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/lineread/dbread_pands.c */
extern int databasename P_((char **name));
extern long dbread_pands P_((char *filename, struct linedb **lines, float wlbeg, float wlend, char *Zfilename, double ***Z, double **T, double **isomass, double ***isocs, int *nT, int *nIso, char ***isonames));

#undef P_
