#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* argum.c */
extern int processparameters P_((int argc, char **argv, struct transit *tr));
extern int acceptsoltype P_((transit_ray_solution **sol, char *hname));
extern int acceptgenhints P_((struct transit *tr));

#undef P_
