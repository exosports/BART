#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/lineread/dbread_pands.c */
extern long dbread_pands P_((char *filename, struct linedb **lines, float wlbeg, float wlend, char *Zfilename, double ***Z, double **isomass, double ***isocs, double **T, int *nT, int *nIso, char ***isonames));

#undef P_
