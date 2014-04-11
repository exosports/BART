#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/atmosphere/readatm.c */
extern void freemem_onept P_((struct onept *o));
extern int freemem_atmosphere P_((struct atm_data *at, long *pi));
extern void saveonept_arr P_((FILE *out, struct onept *onept));
extern int restonept_arr P_((FILE *in, struct onept *onept));
extern double checkaddmm P_((double *mm, long r, prop_mol *molec, struct molecules *mol, int n, short mass));
extern void telldefaults P_((struct isotopes *iso, struct atm_data *at));
extern int getatm P_((struct transit *tr));

#undef P_
