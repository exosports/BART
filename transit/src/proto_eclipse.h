#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/eclipse.c */
extern int tau_eclipse P_((struct transit *tr));
extern void printecl P_((struct transit *tr));
extern int emergent_intens P_((struct transit *tr));

#undef P_
