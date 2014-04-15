#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/argum.c */
extern int processparameters P_((int argc, char **argv, struct transit *tr));
extern int acceptsoltype P_((transit_ray_solution **sol, char *hname));
extern int accepteclipsetype P_((eclipse_ray_solution **ecl, char *hname));
extern int acceptgenhints P_((struct transit *tr));
extern void savehint P_((FILE *out, struct transithint *hints));
extern int resthint P_((FILE *in, struct transithint *hint));
extern void printintro P_((void));
extern void freemem_hints P_((struct transithint *h));
extern void freemem_cloud P_((struct extcloud *c));
extern void freemem_detailout P_((struct detailout *d));
extern void freemem_detailfld P_((struct detailfld *f));

#undef P_
