#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* test_slantpath.c */
extern double calcex P_((double alpha, double rm, double ri));
extern int tau_now P_((double ip, double *rad, double *refr, double *ex, long nrad, double res, int type, char *desc));
extern int tau_dens P_((double rm, double ip, long nrad, double alpha, double res));
extern int tau_ip P_((double rm, double ip, long *nrad, int nr, double alpha, double (*fcn)(void)));
extern int tau_rad P_((double rm, double *ip, int ni, long *nrad, int nr, double alpha, double (*fcn)(void)));
extern double const_ex_anal P_((double alpha, double rm, double ip));
extern double incout_ex_anal P_((double alpha, double rm, double ip));
extern double incin_ex_anal P_((double alpha, double rm, double ip));
extern int tau_ex P_((double *rmax, int nx, double *ip, int ni, long *nrad, int nr, double alpha, char *str));
extern double *mod_ctau P_((double prm, double *res, double star, double ipmax, double first, long nip, double toomuch));
extern double *mod_itau P_((double prm, double *res, double star, double ipmax, double first, long nip, double toomuch));
extern double *mod_dtau P_((double prm, double *res, double star, double ipmax, double first, long nip, double toomuch));
extern int mod_now P_((double *tau, long last, double toomuch, prop_samp *ip, struct geometry *sg, double res, int type, char *desc));
extern int mod_nip P_((struct geometry *sg, double ipmax, double first, long nip, double toomuch, double res, double *(*tauf)(void), double tprm));
extern int mod_first P_((struct geometry *sg, double ipmax, double first, long *nip, int nnip, double toomuch, double *(*tauf)(void), double tprm));
extern int mod_ip P_((struct geometry *sg, double ipmax, double *first, int nfirst, long *nip, int nnip, double toomuch, double *(*tauf)(void), double tprm));
extern int mod_star P_((double star, double *ipmax, int nipmax, double *first, int nfirst, long *nip, int nnip, double toomuch, double *(*tauf)(void), double tprm));
extern int mod_tau P_((double *star, int nstar, double *ipmax, int nipmax, double *first, int nfirst, long *nip, int nnip, double toomuch, double *(*tauf)(void), char *desc, double tprm));
extern int test_mod P_((void));
extern int test_tau P_((void));
extern int main P_((int argc, char *argv[]));

#undef P_
