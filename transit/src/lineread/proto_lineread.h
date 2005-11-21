#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/lineread/lineread.c */
extern int main P_((int argc, char *argv[]));

#undef P_
