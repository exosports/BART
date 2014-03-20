/*
 * readlineinfo.c   - reads line info as returned by
 *                    lineread. Component of Transit program
 *
 * Copyright (C) 2003-2006 Patricio Rojo (pato@astro.cornell.edu)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of version 2 of the GNU General 
 * Public License as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.

  List of functions defined here:

static inline void datafileBS(FILE *fp, PREC_NREC offs, PREC_LNDATA target,
                              PREC_NREC *resultp, int reclength)

int readlineinfo(struct transit *tr)

int readinfo_tli(struct transit *tr, struct lineinfo *li)

int readtli_bin(FILE *fp, struct transit *tr, struct lineinfo *li)

int readtli_ascii(FILE *fp, struct transit *tr, struct lineinfo *li)

int getinifinasctli(double *ini, double *fin, FILE *fp, char *file)

int checkrange(struct transit *tr, struct lineinfo *li)

int readdatarng(struct transit *tr, struct lineinfo *li)

static void notyet(int lin, char *file)

static int invalidfield(char *line, char *file, int nmb, int fld, char *fldn)

int freemem_isotopes(struct isotopes *iso, long *pi)

int freemem_lineinfotrans(struct lineinfo *li, long *pi)

void saveline(FILE *fp, struct lineinfo *li)

int main(int argc, char **argv)
 */


#include <transit.h>

#define TLI_WAV_UNITS 1e-4 /* TLI wavelength (microns, as of v4) */
#define TLI_E_UNITS   1    /* TLI Elow units (cm-1, as of v4)    */

static double tli_to_microns = TLI_WAV_UNITS/1e-4;

static void notyet(int lin, char *file);

static int invalidfield(char *line, char *file, int nmb,
                        int fld, char *fldn);

/* Check if pre condition is satisfied (return error if True).
   Advance pointer while omit condition is satisfied.
   Check if post condition is satisfied (return error if True).   */
#define checkprepost(pointer, pre, omit, post) do{                         \
   if(pre)                                                                 \
     transiterror(TERR_SERIOUS,                                            \
                  "Pre-condition failed on line %i(%s)\n while reading:\n" \
                  "%s\n\nTLI_Ascii format most likely invalid\n",          \
                  __LINE__, __FILE__, line);                               \
   while(omit)                                                             \
     pointer++;                                                            \
   if(post)                                                                \
     transiterror(TERR_SERIOUS,                                            \
                  "Post-condition failed on line %i(%s)\n while reading:\n"\
                  "%s\n\nTLI_Ascii format most likely invalid\n",          \
                  __LINE__, __FILE__, line);                               \
                                             }while(0)


/* \fcnfh
  Do a binary search in file pointed by 'fp' between 'off' and 'off+nfields'
  looking for 'target' as the first item of a record of length 'reclength', 
  result index (with respect to offs) is stored in 'resultp'.       */
static inline void 
datafileBS(FILE *fp,            /* File pointer                         */
           PREC_NREC offs,      /* Initial position of data in tli file */
           PREC_NREC nfields,   /* Number of fields to search           */
           PREC_LNDATA target,  /* Target value                         */
           PREC_NREC *resultp,  /* Result index                         */
           int reclength){      /* Total length of record               */
  /* Variable to keep wavelength: */
  PREC_LNDATA temp;
  const int trglength = sizeof(PREC_LNDATA);
  /* Search boundaries: */
  PREC_NREC ini=0,
            fin=nfields-1;

  transitDEBUG(21, verblevel,
               "BS: Start looking from %li in %li fields for %f\n",
               offs, nfields, target);
  do{
    *(resultp) = (fin+ini)/2;                       /* Middle record's index */
    fseek(fp, offs + reclength*(*resultp), SEEK_SET); /* Go there            */
    fread(&temp, trglength, 1, fp);                   /* Read value          */
    transitDEBUG(21, verblevel, "BS: found wl %f microns at position %li\n",
                 temp*tli_to_microns, (*resultp));
    /* FINDME: No reversebytes? */
    /* Re-set lower or higher boundary: */
    if(target>temp)     
      ini = *(resultp);
    else
      fin = *(resultp);
  }while (fin-ini>1);
  *resultp = ini;
}


/* \fcnfh
  Read initial and final wavelength limits and number of databases.
  Allocate pointers to databases, and isotopes arrays.  Get databases
  info: names, number of temperatures, temperatures, number of
  isotopes, isotope names and masses, partition function, and cross
  sections. Get cumulative number of isotopes.
  Returns 0 on success                                                  */
int 
readtli_bin(FILE *fp, 
            struct transit *tr,
            struct lineinfo *li){
  /* Declare varables:                                                      */
  double iniw, finw;   /* Initial and final wavelength of database          */
  unsigned short ndb;  /* Number of databases                               */
  unsigned short rs;   /* FINDME: auxiliary what??                          */
  unsigned short nT,   /* Number of temperatures per database               */
                 nIso; /* Number of isotopes per database                   */
  PREC_ZREC *T, *Z;    /* Auxiliary temperature and part. function pointers */
  PREC_CS *CS;         /* Auxiliary cross section pointer                   */
  int acumiso=0;       /* Cumulative number of isotopes per database        */
  int correliso=0;     /* Isotopes correlative number                       */
  struct isotopes *iso=tr->ds.iso;

  /* Read TLI version, lineread version, and lineread revision number: */
  fread(&li->tli_ver, sizeof(unsigned short), 1, fp);
  fread(&li->lr_ver,  sizeof(unsigned short), 1, fp);
  fread(&li->lr_rev,  sizeof(unsigned short), 1, fp);
  /* Check compatibility of versions:                                  */
  if(li->tli_ver != compattliversion)
    transiterror(TERR_SERIOUS,
                 "The version of the TLI file: %i (lineread v%i.%i) is not "
                 "compatible with this version of transit, which can only "
                 "read version %i.\n", li->tli_ver, li->lr_ver,
                 li->lr_rev, compattliversion);

  /* Read initial wavelength, final wavelength, and number of databases: */
  fread(&iniw, sizeof(double), 1, fp);
  fread(&finw, sizeof(double), 1, fp);
  fread(&rs,   sizeof(unsigned short), 1, fp);
  char undefined_string[rs+1]; /* FINDME: Can declare at the beginning?  */
  fread(undefined_string, sizeof(char), rs, fp);
  undefined_string[rs] = '\0'; /* FINDME: undefined_string is not used */
  fread(&ndb, sizeof(unsigned short), 1, fp);

  /* Allocate pointers according to the number of databases:              */
  iso->db     = (prop_db      *)calloc(ndb, sizeof(prop_db     ));
  li->db      = (prop_dbnoext *)calloc(ndb, sizeof(prop_dbnoext));
  /* Allocate isotope info.  Start with size 1, then reallocate as needed */
  iso->isof   = (prop_isof    *)calloc(1,   sizeof(prop_isof));
  li->isov    = (prop_isov    *)calloc(1,   sizeof(prop_isov));
  li->isov->z = (double       *)calloc(1,   sizeof(double));
  li->isov->c = (double       *)calloc(1,   sizeof(double));

  /* Read info for each database: */
  for(short i=0; i<ndb; i++){
    /* Allocate and get DB's name: */
    fread(&rs, sizeof(unsigned short), 1, fp);         /* Get DB name length */
    iso->db[i].n = (char *)calloc(rs+1, sizeof(char)); /* Allocate           */
    fread(iso->db[i].n, sizeof(char), rs, fp);         /* Read               */
    iso->db[i].n[rs] = '\0';
    
    /* Get number of temperatures and isotopes: */
    fread(&nT,   sizeof(unsigned short), 1, fp);
    fread(&nIso, sizeof(unsigned short), 1, fp);
    li->db[i].t  = nT;
    iso->db[i].i = nIso;

    /* Allocate for the different temperature points and read: */
    T = li->db[i].T = (double *)calloc(nT, sizeof(double));
    fread(T, sizeof(double), nT, fp);

    /* Reallocate memory to account for new isotopes: */
    li->isov  = (prop_isov *)realloc(li->isov,
                                      (correliso+nIso)*sizeof(prop_isov));
    iso->isof = (prop_isof *)realloc(iso->isof,
                                      (correliso+nIso)*sizeof(prop_isof));
    li->isov[correliso].z = (double *)calloc((correliso+nIso)*nT,
                                             sizeof(double));
    li->isov[correliso].c = (double *)calloc((correliso+nIso)*nT,
                                             sizeof(double));
    transitDEBUG(21, verblevel,
                 "So far, cumIsotopes: %i, at databases: %i, position %li.\n",
                 correliso+nIso, i, ftell(fp)); 

    /* Display database general info: */
    transitDEBUG(23, verblevel,
                 "DB %i: \"%s\" has %i (%i) temperatures, %i (%i) isotopes, "
                 "and starts at cumulative isotope %i.\n",
                 iso->isof[correliso].d, iso->db[i].n,
                 li->db[i].t, nT, 
                 iso->db[i].i, nIso, iso->db[i].s);

    for (unsigned int j=0; j<nIso; j++){
      transitDEBUG(22, verblevel, "isotope %i/%i for DB %i.\n", j+1, nIso, i);

      /* Store isotopes'  DB index number:                    */
      iso->isof[correliso].d = i;

      /* Read isotopes' name: */
      fread(&rs, sizeof(unsigned short), 1, fp);
      iso->isof[correliso].n = (char *)calloc(rs+1, sizeof(char));
      fread(iso->isof[correliso].n, sizeof(char), rs, fp);
      iso->isof[correliso].n[rs] = '\0';
      transitDEBUG(21, verblevel, "  Name's length: %i, position: %li, value: "
                                "%s.\n", rs, (long)(ftell(fp)), iso->isof[i].n);

      /* Read mass: */
      fread(&iso->isof[correliso].m, sizeof(double), 1, fp);
      transitDEBUG(21, verblevel,"  Mass read: %g * %g = %g, position: %li, "
                   "size %i.\n", iso->isof[i].m, AMU, iso->isof[i].m*AMU,
                   (long)(ftell(fp)), (int)(sizeof(iso->isof[i].m)));

      /* Read partition function and cross section: */
      Z  = li->isov[correliso].z = li->isov[correliso-j].z+nT*j;
      fread(Z,  sizeof(double), nT, fp);
      CS = li->isov[correliso].c = li->isov[correliso-j].c+nT*j;
      fread(CS, sizeof(double), nT, fp);
      li->isov[correliso].n = nT;

      transitDEBUG(12, verblevel, "Z(%i/%i):%g %g ... %g.\n", 
                                  j+1, nIso, Z[0], Z[1], Z[nT-1]);
      correliso++;
    }

    /* Update cumulative isotope count: */
    iso->db[i].s = acumiso;
    acumiso += nIso;

    /* Read database correlative number: */
    fread(&rs, sizeof(unsigned short), 1, fp);
    if (i != rs)
      transiterror(TERR_SERIOUS,
                   "Problem in TLI file: database correlative number (%i) "
                   "doesn't match information read (%i)\n"
                   "Isotopes read: %i\n"
                   "Last DB #temps: %i\n"
                   "Last DB #iso: %i\n",
                   i, rs, acumiso, nT, nIso);
  }

  /* Read total number of isotopes: */
  fread(&iso->n_i, sizeof(unsigned short), 1, fp);
  if(iso->n_i != acumiso)
    transiterror(TERR_SERIOUS,
                 "Read number of isotopes (%i), doesn't match the "
                 "total number of isotopes (%i).\n",
                 iso->n_i, acumiso);

  /* Update structure values: */
  li->ni  = iso->n_i;      /* Number of isotopes                  */
  li->ndb = ndb;           /* Number of databases                 */
  li->endinfo = ftell(fp); /* Pointer position of first line data */
  li->wi = iniw;           /* Initial wavelength                  */
  li->wf = finw;           /* Final wavelength                    */
  iso->isov = (prop_isov *)calloc(iso->n_i, sizeof(prop_isov));
  iso->n_db = ndb;         /* Number of databases                 */

  return 0;
}


/* \fcnfh
    Get number of databases and names.  Get number of temperatures, 
    number of isotopes.  Allocate pointers to databases and isotopes
    arrays.  Get isotope names and masses.  Get temperatures, partition
    function, and cross sections.  Get initial and final wavelength
    limits of databases.  Count cumulative number of isotopes.
    Returns 0 on success.                                               */
int 
readtli_ascii(FILE *fp, 
              struct transit *tr,
              struct lineinfo *li){
  char rc;
  char line[maxline+1], *lp, *lp2;
  int ndb, db;
  unsigned int nIso, nT;
  int acumiso;
  int rn, i;
  prop_isov *isov;
  PREC_ZREC *T;
  struct isotopes *iso=tr->ds.iso;

  /* Format of the TLI-ASCII file is the following (names should not
     contain spaces, '_' are replaced by spaces)
     \begin{verb}
     <m-database>
     <DATABASE1-name> <n1-iso> <nt1-temp>
     <NAME1> <MASS1> ... <NAMEn1> <MASSn1>
     <TEMP1>    <Z1-1>  ...  <Z1-n1>   <CS1-1>  ...  <CS1-n1>
     ...
     <TEMPnt1> <Znt1-1> ... <Znt1-n1> <CSnt1-1> ... <CSnt1-n1>
     <DATABASE2-name> <n2-iso> <nt2-temp>
     <NAME1> <MASS1> ... <NAMEn2> <MASSn2>
     <MASS1> ... <MASSn2>
     <TEMP1>    <Z1-1>  ...  <Z1-n2>   <CS1-1>  ...  <CS1-n2>
     ...
     <TEMPnt2> <Znt2-1> ... <Znt2-n2> <CSnt2-1> ... <CSnt2-n2>
     ....
     ....
     <DATABASEm-name> <nm-iso> <ntm-temp>
     <NAME1> <MASS1> ... <NAMEnm> <MASSnm>
     <TEMP1>    <Z1-1>  ...  <Z1-nm>   <CS1-1>  ...  <CS1-nm>
     ...
     <TEMPntm> <Zntm-1> ... <Zntm-nm> <CSntm-1> ... <CSntm-nm>
     <CTRWAV1> <ISOID1> <LOWENER1> <LOGGF1>
     <CTRWAV2> <ISOID2> <LOWENER2> <LOGGF2>
     ...
     ...
     \end{verb}                                                      */

  /* Get Number of databases from first line:                        */
  /* Skip comments and empty lines:                                  */
  while((rc=fgetupto_err(line, maxline, fp, &linetoolong, tr->f_line,
                         li->asciiline++))=='#' || rc=='\n');
  if(!rc) /* EOF found */
    notyet(li->asciiline, tr->f_line);
  ndb = strtol(line, &lp, 0);
  fprintf(stderr, "%i", ndb);
  /* Check for errno and ERANGE flags, then check that lp now
     points to \0 (skipping blank and tabs)                          */
  checkprepost(lp, errno&ERANGE, *lp==' ' || *lp=='\t', *lp!='\0');

  /* Allocate pointers according to the number of databases: */
  iso->db = (prop_db      *)calloc(ndb, sizeof(prop_db));
  li->db  = (prop_dbnoext *)calloc(ndb, sizeof(prop_dbnoext));

  /* For each database:                                              */
  for(db=0; db<ndb; db++){
    /* Get name:                                                     */
    while((rc=fgetupto_err(lp=line, maxline, fp, &linetoolong, tr->f_line,
                           li->asciiline++)) == '#' || rc=='\n');
    if(!rc) notyet(li->asciiline, tr->f_line);
    /* Read name and store it, replacing underscores by blank spaces */
    if((iso->db[db].n = readstr_sp_alloc(lp, &lp, '_'))==NULL)
      transitallocerror(0);

    /* Go to next field:                                             */
    checkprepost(lp, 0, *lp==' ' || *lp=='\t', *lp=='\0');
    /* Get number of temperatures and isotopes:                      */
    rn = getnl(2, ' ', lp, &nIso, &nT);
    /* Check two numbers were read, and they are != 0:               */
    checkprepost(lp, rn!=2, 0, nIso==0 || nT==0);
    /* Store the values:                                             */
    li->db[db].t  = nT;
    iso->db[db].i = nIso;

    /* Update acumulated isotope count:                              */
    iso->db[db].s = iso->n_i;
    iso->n_i += nIso;

    /* Allocate for variable and fixed isotope info as well as for
       temperature points:                                           */
    T = li->db[db].T = (PREC_ZREC *)calloc(nT, sizeof(PREC_ZREC));

    /* Allocate structure to receive the isotope info:               */
    if(!db){  /* calloc if this is the first database                */
      isov = li->isov = (prop_isov  *)calloc(iso->n_i, sizeof(prop_isov));
      iso->isof       = (prop_isof  *)calloc(iso->n_i, sizeof(prop_isof));
      //iso->isov       = (prop_isov  *)calloc(iso->n_i, sizeof(prop_isov));
    }
    else{     /* Else, realloc                                       */
      isov = li->isov =
                 (prop_isov  *)realloc(li->isov,  iso->n_i*sizeof(prop_isov));
      iso->isof =(prop_isof  *)realloc(iso->isof, iso->n_i*sizeof(prop_isof));
      //iso->isov =(prop_isov  *)realloc(iso->isov, iso->n_i*sizeof(prop_isov));
    }

    /* Allocate cross section and partition function.  Set isov at
       the beginning of this database:                               */
    isov   += (acumiso=iso->db[db].s);
    isov->z = (PREC_ZREC *)calloc(nIso*nT, sizeof(PREC_ZREC));
    isov->c = (PREC_CS   *)calloc(nIso*nT, sizeof(PREC_CS));

    /* Get isotope name and mass: */
    while((rc=fgetupto_err(lp2=line, maxline, fp, &linetoolong, tr->f_line,
                           li->asciiline++))=='#' || rc=='\n');
    if(!rc) notyet(li->asciiline, tr->f_line);
    /* For each isotope:                                        */
    for(i=0; i<nIso; i++){
      isov[i].z = isov->z+nT*i;
      isov[i].c = isov->c+nT*i;
      /* Get name:                                              */
      if((iso->isof[acumiso+i].n = readstr_sp_alloc(lp2, &lp, '_'))==NULL)
        transitallocerror(0);
      /* Get mass and convert to cgs:                           */
      /* FINDME: where is being converted to cgs?               */
      iso->isof[acumiso+i].m = strtod(lp, &lp2);
      if(i != nIso-1)
        checkprepost(lp2, lp==lp2, *lp2==' ' || *lp2=='\t', *lp2=='\0');
    }
    /* Last isotope has to be followed by an end of string:     */
    checkprepost(lp2, 0, *lp2==' ' || *lp2=='\t', *lp2!='\0');

    /* Get temperatures, partition function, and cross section: */
    for(rn=0; rn<nT; rn++){
      while((rc=fgetupto_err(lp=line, maxline, fp, &linetoolong, tr->f_line,
                             li->asciiline++)) == '#' || rc == '\n');
      if(!rc) notyet(li->asciiline, tr->f_line);
      while(*lp == ' ')  /* Skip empty spaces                   */
        lp++;
      /* Read temperature:                                      */
      T[rn] = strtod(lp, &lp);
      checkprepost(lp, *lp=='\0', *lp==' ' || *lp=='\t', *lp=='\0');
      /* For each isotope in database, read partition function: */
      for(i=0; i<nIso; i++){
        isov[i].z[rn] = strtod(lp, &lp);
        checkprepost(lp, *lp=='\0', *lp==' ' || *lp=='\t', *lp=='\0');
      }
      /* Read cross section:                                    */
      for(i=0; i<nIso-1; i++){
        isov[i].c[rn] = strtod(lp, &lp);
        checkprepost(lp, *lp=='\0', *lp==' ' || *lp=='\t', *lp=='\0');
      }
      /* The last field needs a different check at the end:     */
      isov[i].c[rn] = strtod(lp, &lp);
      checkprepost(lp, 0, *lp==' ' || *lp=='\t', *lp!='\0');
    }
  }

  /* Store number of database and position of first line data:  */
  iso->n_db = ndb;
  li->endinfo = ftell(fp);

  /* Look for the ascii storage range:                          */
  getinifinasctli(&li->wi, &li->wf, fp, tr->f_line);

  return 0;
}


/* \fcnfh
   Find initial and final wavelength of ASCII TLI file
   @returns 0 on success                                            */
int
getinifinasctli(double *ini, /* where initial value would be stored */
                double *fin, /* where initial value would be stored */
                FILE *fp,    /* File pointer to first line record   */
                char *file){ /* Name of file being read             */

  char rc;
  char line[maxline+1], *lp, *lplast;

  /* Go to first line:                                        */
  while((rc=fgetupto_err(lp=line, maxline, fp, &linetoolong,
                         file, 0)) =='#' || rc=='\n');
  /* There should be at least one isotope:                    */
  if(!rc){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "readlineinfo:: There was no transition info in file '%s', "
                 "only general isotope info.\n", file);
    exit(EXIT_FAILURE);
  }
  /* Read first value: */
  *ini = strtod(line, &lp);
  if(line==lp){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "readlineinfo:: First wavelength transitions in "
                 "file '%s' is not valid in line:\n%s\n", file, line);
    exit(EXIT_FAILURE);
  }
  /* Go to end of file: */
  fseek(fp, 0, SEEK_END);
  if(ftell(fp)<maxline){
    transiterror(TERR_WARNING,
                 "readlineinfo:: weird, TLI-Ascii file has less than %i "
                 "bytes.  That looks improbable.\n", maxline);
    fseek(fp, 0, SEEK_SET);
  }
  else
    fseek(fp, 1-maxline, SEEK_END); /* Go to last record: */
  fread(lp=line, sizeof(char), maxline-1, fp);
  line[maxline-1] = '\0';

  /* Set pointer to the end of record: */
  lplast = NULL;
  while(*lp){
    lp++;
    if(*lp!='\0' && lp[-1]=='\n')
      lplast=lp;
  }
  if(!lplast){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Last line in '%s' is longer than %i bytes.\n", file, maxline);
    exit(EXIT_FAILURE);
  }
  /* Read wavelength of last record: */
  *fin = strtod(lplast, &lp);
  if(line==lp){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "readlineinfo:: Last central wavelength of transitions in "
                 "file '%s' is not valid in line:\n%s\n", file, lplast);
    exit(EXIT_FAILURE);
  }

  return 0;
}
/* TD: data info in a structure */

/* Set the molecule's index of the isotopes: */
int
setimol(struct transit *tr){
  struct molecules *mol = tr->ds.mol;
  struct isotopes  *iso = tr->ds.iso;
  int i; /* Auxiliary for-loop indices */

  iso->imol = (int *)calloc(iso->n_i, sizeof(int));
  for (i=0; i<iso->n_i; i++){
    //transitprint(1, verblevel, "Isotope %d is '%s', from DB %d: '%s'.\n",
    //                  i, iso->isof->n, iso->isof->d, iso->db[iso->isof->d].n);
    /* If the database is P&S, it's a water isotope: */
    if (strcmp(iso->db[iso->isof->d].n, "Partridge & Schwenke (1997)\0")==0){
      /* Search water molecule's index:              */
      iso->imol[i] = findstring("H2O\0", mol->name, mol->nmol);
      transitprint(30, verblevel, "Isotope '%s', is mol %d: '%s'.\n",
                      iso->isof[i].n, iso->imol[i], mol->name[iso->imol[i]]);
    }
  }
  return 0;
}

int
getisoratio(struct transit *tr){
  struct isotopes *iso = tr->ds.iso;
  int i, j,
      maxlinelen=501,
      niratio=4;  /* Number of isotope ratios in list       */
  double *iratio; /* Isotope number-abundance ratio in list */
  char **iname;   /* Isotope name in list                   */
  char line[maxlinelen], *lp; /* string pointers            */
  char *filename = "../inputs/molecules.dat";
  FILE *elist;

  iso->isoratio = (double *)calloc(iso->n_i, sizeof(double));

  iratio   = (double *)calloc(niratio,             sizeof(double));
  iname    = (char  **)calloc(niratio,             sizeof(char *));
  iname[0] = (char   *)calloc(niratio*maxeisoname, sizeof(char));
  for (i=1; i<niratio; i++)
    iname[i] = iname[0] + i*maxeisoname;

  /* Open Molecules file: */
  if((elist=verbfileopen(filename, "Molecular info ")) == NULL)
    exit(EXIT_FAILURE);

  do{
    lp = fgets(line, maxlinelen, elist);
  }while (lp[0] == '\0' || lp[0] == '\n' || lp[0] == '#');
  do{
    lp = fgets(line, maxlinelen, elist);
  }while (lp[0] != '\0' && lp[0] != '\n' && lp[0] != '#');

  do{
    lp = fgets(line, maxlinelen, elist);
  }while (lp[0] == '\0' || lp[0] == '\n' || lp[0] == '#');
  do{
    lp = fgets(line, maxlinelen, elist);
  }while (lp[0] != '\0' && lp[0] != '\n' && lp[0] != '#');

  do{
    lp = fgets(line, maxlinelen, elist);
  }while (lp[0] == '\0' || lp[0] == '\n' || lp[0] == '#');
  do{
    lp = fgets(line, maxlinelen, elist);
  }while (lp[0] != '\0' && lp[0] != '\n' && lp[0] != '#');

  do{
    lp = fgets(line, maxlinelen, elist);
  }while (lp[0] == '\0' || lp[0] == '\n' || lp[0] == '#');

  /* Read list of isotope ratios:                       */
  for (i=0; i<niratio; i++){
    /* FINDME: For now skip molecule name:              */
    lp = nextfield(lp);
    getname(lp, iname[i]);        /* Read isotope name  */
    //transitprint(1, verblevel, "word: %s\n", iname[i]);
    lp = nextfield(lp);
    iratio[i] = strtod(lp, NULL); /* Read isotope ratio */
    //transitprint(1, verblevel, "number: %.5f\n", iratio[i]);
    lp = fgets(line, maxlinelen, elist);
  }

  /* Find respective isotope name in list */
  for (i=0;  i < iso->n_i;  i++){
    j = findstring(iso->isof[i].n, iname, niratio);
    //transitprint(1, verblevel, "jota: %d, name: %s.\n", j, iso->isof[i].n);
    iso->isoratio[i] = iratio[j];
  }
  return 0;
}

/* \fcnfh
   Initialize wavelength sample struct.  Set margin.  
   Set initial and final wavelengths to use.  
   Check that margin leaves a non-zero wavelength range.

   Return:  0   if all hinted values were accepted, else
           (positive if something was changed):
           0x20 if hinted final wavelength was changed
           0x10 if hinted initial wavelength was changed
           (negative: If something bad happen):
           -1   if hinted initial wavelength is larger than sug. final
           -2   if hinted initial is larger than largest allowed final
           -3   if hinted final is shorter than smallest allowed initial
           -4   if margin value was too big                             */
int
checkrange(struct transit *tr,   /* General parameters and  hints    */
           struct lineinfo *li){ /* Values returned by readinfo_tli  */

  int res=0;                           /* Return value                   */
  PREC_RES margin;                     /* Margin value                   */
  struct transithint *th = tr->ds.th;  /* transithint                    */
  prop_samp *msamp = &li->wavs;        /* transit wavelength sampling    */
  prop_samp *hsamp = &th->wavs;        /* hint    wavelength sampling    */
  PREC_LNDATA dbini = li->wi*TLI_WAV_UNITS, /* Minimum db wavelength     */
              dbfin = li->wf*TLI_WAV_UNITS; /* Maximum db wavelength     */
  double extra;    /* FINDME */

  /* FINDME: hack prints: */
  //transitprint(1, verblevel, " hsamp->f: %g,  hsamp->i: %g\n",
  //               hsamp->f, hsamp->i);
  //transitprint(1, verblevel, " hsamp->fct: %g\n", hsamp->fct);
  //transitprint(1, verblevel, " db i: %g,  db f: %g\n", dbini, dbfin);

  /* Initialize modified hints: */
  msamp->n = -1;
  msamp->d = -1;
  msamp->v = NULL;
  msamp->fct = 1;

  /* Check that factor is positive and non-zero, set it:        */
  if(hsamp->fct<0)
    transiterror(TERR_SERIOUS, "User specified wavelength factor is "
                               "negative (%g).\n", hsamp->fct);
  /* Set lineinfo wavelength units factor equal to transithint factor: */
  if(hsamp->fct>0)
    msamp->fct = hsamp->fct;

  /* transit lineinfo.wavs conversion factor to cgs: */
  double fct = msamp->fct;

  /* fct to microns conversion factor:     */
  double fct_to_microns = msamp->fct/1e-4;

  // hack print:
  // transitprint(1, verblevel, " margin: %g\n", th->margin);

  /* Check that the margin leaves a non-zero range in the line dataset,
     set it:                                                             */
  /* FINDME: The range of the lineinfo should not interfere with the 
     transit's wavelength sampling. There should be a warning note if the
     database limits are smaller than the transit range, though.         */
  if(2*th->margin*fct > (dbfin - dbini)){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Margin value (%g microns) is too big for this dataset "
                 "whose range is %g to %g microns. Factor to convert user "
                 "margin (%g) to centimeters is %g.\n",
                 th->margin*fct_to_microns, li->wi*tli_to_microns,
                 li->wf*tli_to_microns, th->margin, msamp->fct);
    return -4;
  }
  /* Set transit margin in cgs units: */
  margin = tr->margin = th->margin*msamp->fct;

  /* If TLI file is ASCII, double the margin limits. This is because,
     as opposed to binary archiving, I'm not guaranteed that I have all
     the lines between 'dbini' and 'dbfin', instead they are the minimum
     and maximum wavelength of the given transitions.                     */
  if(li->asciiline){
    if(margin==0.0)
      transiterror(TERR_WARNING,
                   "Wavelength margin is zero in a TLI-ASCII file. "
                   "There will be no points to the left or "
                   "right of the extreme central wavelengths.\n");
    extra = 2*margin;
  }
  else
    extra = 0;

  transitDEBUG(21, verblevel,
               "Hinted initial and final wavelengths are %g and %g cm.\n"
               "Databse's max and min wavelength are %g and %g cm.\n",
               hsamp->i*fct, hsamp->f*fct, dbini, dbfin);

  /* Set final wavelength:                                         */
  /* If invalid/not set hint final wavelength, default it to zero: */
  if(hsamp->f<0){
    hsamp->f = 0;
    transiterror(TERR_WARNING,
                 "Incorrect upper wavelength limit in hint.  Default: setting "
                 "to %g before extraction.\n", hsamp->f*fct);
  }
  /* If hint is 0, set it to max db wavelength:                    */
  if(hsamp->f <= 0){
      msamp->f = (dbfin + extra)/fct;
  }
  else{  /* Else, hinted f is a positive value:                    */
    transitDEBUG(20, verblevel, "dbini: %g  margin: %g  sampf: %g.\n",
                 dbini, margin, hsamp->f);
    /* Check that it is not below the minimum value:               */
    if(dbini+margin > fct*hsamp->f){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                   "Considering margin, final wavelength (%g * %g) is smaller "
                   "than minimum wavelength in database (%g = %g + %g).\n",
                   hsamp->f, fct, dbini+margin, dbini, margin);
      return -3;
    }
    // hack print:
    //transitprint(1, verblevel, " hsamp->f: %g,  margin: %g,  dbfin: %g\n", 
    //             hsamp->f, margin, dbfin);
    /* Warn if it is above maximum value with information:         */
    if(hsamp->f*fct+margin > dbfin)
      transiterror(TERR_WARNING,
                   "Final requested wavelength (%g microns) is larger "
                   "than the maximum informative value in database "
                   "(%g microns).\n", hsamp->f, dbfin*fct_to_microns);
    /* Set the final wavelength value:                             */
    msamp->f = hsamp->f;
  }
  /* Set initial wavelength:                                       */
  /* If invalid value, default it to 0: */ 
  if(hsamp->i < 0){
    hsamp->i = 0;
    transiterror(TERR_WARNING,
                 "Setting hinted lower wavelength limit before "
                 "extraction as %g cgs. It was not user-hinted.\n",
                 hsamp->i*fct);
  }
  /* If default value, set it to min db wavelength:                */
  if(hsamp->i<=0)
    msamp->i = (dbini - extra)/fct;
  else{
    transitDEBUG(20, verblevel, "dbfin: %g  margin: %g  sampi: %g.\n",
                 dbfin, margin, fct*hsamp->i);
    /* Check that it is not larger than the maximum db wavelength: */
    if(dbfin < margin+fct*hsamp->i){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                   "Initial wavelength (%g cm) is larger than maximum "
                   "wavelength in database (%g cm = %g + %g cm).\n",
                   fct*hsamp->i, dbfin-margin, dbfin, margin);
      return -2;
    }
    if(fct*hsamp->i-margin < dbini)
      transiterror(TERR_WARNING,
                   "Initial requested wavelength (%g microns) is smaller "
                   "than the minimum informative value in database "
                   "(%g microns).\n", hsamp->i, dbini*fct_to_microns);
    msamp->i = hsamp->i;
  }

  /* Check that we still have a non-zero wavelength range remaining: */
  if(2*margin > (msamp->f-msamp->i)*fct){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Usable final (%g cm) has to be larger than usable initial "
                 "wavelength (%g cm). Note that those values could have been "
                 "modified according to the database range (%g - %g cm) and "
                 "margin (%g cm).\n", fct*msamp->i+margin,
                 fct*msamp->f-margin, dbini, dbfin, margin);
    return -1;
  }

  /* Set progress indicator and return status:                     */
  tr->pi |= TRPI_CHKRNG;
  return res;
}


/* \fcnfh
    Outputs error and exit the program when EOF is found before
    expected                                                     */
static void
notyet(int lin, char *file){
  transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
               "readlineinfo:: EOF unexpectedly found at line %i in "
               "ascii-TLI linedb info file '%s'.\n", lin, file);
  exit(EXIT_FAILURE);
}


/* \fcnfh
   Check TLI file exists.  Check that machine formating is compatible
   with lineread.  Determine if TLI is ASCII or binary.  Read either
   ASCII or binary TLI file. Declare line_transition.

  TD:  Checks on allocation errors.
  Return: 1 on success
         -1 unavailable file
         -2 Filename not hinted
         -3 TLI format not valid (missing magic bytes)
         -4 Improper TLI-ASCII input                              */
int
readinfo_tli(struct transit *tr,
             struct lineinfo *li){
  int rn;
  FILE *fp;  /* File pointer of info file: */

  /* Decalre and initialize the union sign:                         */
  /* sign.s contains the magic numbers of this machine's and TLI's: */
  union {char sig[4]; int32_t s[2];} sign =
    {.s={0, ((0xff-'T')<<24)|((0xff-'L')<<16)|((0xff-'I')<<8)|(0xff)}};
  char line[maxline+1];

  /* Pointer to hint: */
  struct transithint *th = tr->ds.th;

  /* Get TLI file name from hint:                              */
  if(!th->f_line){  /* Check if it was defined in hint         */
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT, "Undefined TLI file name.\n");
    return -2;
  }
  /* Check that the file exists and make a pointer to read it: */
  if((rn=fileexistopen(th->f_line, &fp)) != 1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Line info file '%s' is not found. "
                 "fileexistopen() error code %i.\n", th->f_line, rn);
    return -1;
  }
  /* Set transit TLI file pointer and TLI file name: */
  tr->fp_line = fp;
  tr->f_line  = th->f_line;

  /* Read first four bytes, they should be either
  `(0xff-T)(0xff-L)(0xff-I)(0xff)' or '\#TLI'. They are stored as integer.  
  This checks whether the machine where the TLI file and the one this
  program is being run have the same endian order.  If the first two are
  '\#TLI', then the first line should also start as '\#TLI-ascii'          */
  fread(sign.s, sizeof(int32_t), 1, fp);

  /* Determine if TLI is a binary (asciiline=0) or ASCII (asciiline=1) file: */
  li->asciiline = 0;
  transitDEBUG(13, verblevel,
               "Comparing %i and %i for Magic Number (len: %li)\n",
               sign.s[0], sign.s[1], sizeof(sign.s[0]));
  
  if(sign.s[0] != sign.s[1]){
    /* Does it look like an ASCII TLI?, if so check it: */
    rn = strncasecmp(sign.sig, "#TLI", 4);  /* FINDME: strncasecmp */
    if(!rn){
      strcpy(line, "#TLI");
      fread(line+4, sizeof(char), 6, fp);
      rn = strncasecmp(line, "#TLI-ascii", 10);
    }
    /* If it wasn't a valid TLI, throw error and exit: */
    if(rn){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                   "The file '%s' has not a valid TLI format. It might be "
                   "because the machine were the file was created have "
                   "different endian order, which is incompatible.\n",
                   tr->f_line);
      return -3;
    }
    li->asciiline = 1;
    /* Ignore the rest of the first line: */
    fgetupto_err(line, maxline, fp, &linetoolong, tr->f_line, 1);
  }

  /* Read ASCII TLI:   */
  if(li->asciiline){
    if((rn=readtli_ascii(fp, tr, li))!=0){
      transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
                   "readtli_ascii() return error code %i.\n", rn);
      return -5;
    }
  }
  /* Read binary TLI:  */
  else
    if((rn=readtli_bin(fp, tr, li))!=0){
      transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
                   "readtli_bin() return error code %i.\n", rn);
      return -6;
    }
  transitprint(3, verblevel, "TLI file read from %g to %g microns.\n",
               li->wi, li->wf);

  /* Declare linetransition struct and set wavelength and lower energy unit 
     factors (As of TLI v4, always in microns and cm-1, respectively):     */
  struct line_transition *lt = &li->lt;
  lt->wfct = TLI_WAV_UNITS;
  lt->efct = TLI_E_UNITS;

  /* Close TLI file pointer, set progres indicator and return success: */
  fclose(fp);
  tr->pi |= TRPI_READINFO;
  return 1;
}


/* \fcnfh
   Print out error details when a line transition record has an invalid
   field 
   @returns -5 always                                                   */
static int
invalidfield(char *line,   /* Contents of the line  */
             char *file,   /* File name             */
             int nmb,      /* File number           */
             int fld,      /* field with the error  */
             char *fldn){  /* Name of the field     */
  transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
               "Line %i of file '%s': Field %i (%s) does not have a valid "
               "value: %s.\n", nmb, file, fld, fldn, line);
  return -5;
}


/*\fcnfh
  Read and store the line transition info (central wavelength, isotope
  ID, lowE, log(gf)) into lineinfo.  Return the number of lines read.

  Return: the number of records read on success, else:
          -1 unexpected EOF
          -2 file non-seekable
          -3 on non-integer number of structure records
          -4 First field is not valid while looking for starting point 
          -5 One of the fields contained an invalid flaoating point     */
int readdatarng(struct transit *tr,   /* transit structure    */
                struct lineinfo *li){ /* lineinfo structure   */

  FILE *fp;          /* Data file pointer                           */
  int rn;            /* Return IDs                                  */
  long i, j;         /* Auxiliary for indices                       */
  PREC_NREC offs=0;  /* Number of lines since the first transition  */
  PREC_NREC alloc=8; /* Size for linetransition arrays allocation   */
  /* Line transition arrays:                          */
  PREC_LNDATA *ltgf;   /* Pointer to log(gf) array    */
  PREC_LNDATA *ltelow; /* Pointer to low-energy array */
  PREC_LNDATA *ltwl;   /* Pointer to wavelength array */
  short *ltisoid;      /* Pointer to isotope ID array */
  /* Auxiliary variables to keep wavelength limits:   */
  PREC_LNDATA iniw = li->wavs.i*li->wavs.fct/TLI_WAV_UNITS;
  PREC_LNDATA finw = li->wavs.f*li->wavs.fct/TLI_WAV_UNITS;
  PREC_LNDATA wltmp;   /* Auxiliary variable to store wavelength */
  PREC_NREC nfields;   /* Number of fields                       */
  char line[maxline+1], rc, *lp, *lp2;

  /* Open line data file : */
  if((rn=fileexistopen(tr->f_line, &fp)) != 1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Data file '%s' not found.  fileexistopen() error "
                 "code: %i.\n", tr->f_line, rn);
    return -1;
  }

  /* Initial allocation for line transition structures: */
  ltgf    = (PREC_LNDATA *)calloc(alloc, sizeof(PREC_LNDATA));
  ltwl    = (PREC_LNDATA *)calloc(alloc, sizeof(PREC_LNDATA));
  ltelow  = (PREC_LNDATA *)calloc(alloc, sizeof(PREC_LNDATA));
  ltisoid = (short       *)calloc(alloc, sizeof(short));

  /* Check for allocation errors:                       */
  if(!ltgf || !ltwl || !ltelow || !ltisoid)
    transiterror(TERR_CRITICAL|TERR_ALLOC,
                 "Couldn't allocate memory for linetran structure array of "
                 "length %i, in function readdatarng.\n", alloc);

  /* Find starting point in datafile:       */ 
  /* For ASCII TLI, do a sequential search: */
  if(li->asciiline){
    /* Go to line where readinfo_tli ended reading: */
    fseek(fp, li->endinfo, SEEK_SET);
    /* Skip all comments and blank lines:           */
    while((rc=fgetupto_err(lp=line, maxline, fp, &linetoolong, tr->f_line,
                           li->asciiline++)) =='#' || rc=='\n');
    /* Set pointer to first linetransition record:  */
    li->endinfo = ftell(fp)-strlen(line)-1;

    /* Read wavelengths until we reach iniw:                         */
    while(!feof(fp)){
      wltmp = strtod(lp, &lp2);
      /* If reading failed:                                          */
      if(lp==lp2){
         transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                     "First field of line %i in file '%s' is not a valid "
                     "floating point value.\n", li->asciiline+offs, tr->f_line);
        return -4;
      }
      /* Is the current line's central wavelength larger than iniw?: */
      if(wltmp>=iniw)
        break;
      /* Skip following comments and increase offs:                  */
      while((rc=fgetupto_err(lp=line, maxline, fp, &linetoolong, tr->f_line,
                             li->asciiline+offs++)) =='#' || rc=='\n');
    }
    /* Go back to first record:                                      */
    /* FINDME: It should stay at current position! */
    fseek(fp, li->endinfo, SEEK_SET);
  }

  /* For binary TLI, do a binary search:  */
  else{
    /* Check seekability: */
    if(fseek(fp, 0, SEEK_END)){
      transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
                   "File '%s' was not seekable when trying to go to the end.\n",
                   tr->f_line);
      return -2;
    }

    /* Find number of fields, checking that there are an integer
       number of them:                                         */
    offs = li->endinfo;  /* Position of first record           */
    j    = ftell(fp);    /* Position of last  record           */
    rn = sizeof(short) + 3*sizeof(PREC_LNDATA); /* Record size */  
    /* FINDME: hardcoded size: rn, can this be predefined? */
    nfields = ((j-offs)/rn); /* Number of records              */
    if(nfields*rn+offs != j){
      transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
                   "Data file does not have an integer number of records. "
                   "Initial byte %i, final %i, record size %i.\n", offs, j, rn);
      /* FINDME: Haven't we already checked this? */
      return -3;
    }
    /* Do binary search in units of TLI:                       */
    datafileBS(fp, offs, nfields, iniw, &j, rn);
    transitDEBUG(21, verblevel, "Beginning found at position %li ", j);
    /* Check if previous records have the same wavelength and 
       adjust offs if necessary:                               */
    if (j){
      do{
        fseek(fp, offs + --j*rn, SEEK_SET);
        fread(&wltmp, sizeof(PREC_LNDATA), 1, fp);
      }while(wltmp>=iniw);
      j++;
    }
    /* Set file pointer at first record to be read: */
    transitDEBUG(21, verblevel, "and then slide to %li.\n", j);
    fseek(fp, offs+j*rn, SEEK_SET);
  }

  /* Main loop to read all the data:       */
  i = 0;
  while(1){
    rn = rc = 1;
    /* Re-allocate if we need more space:  */
    if(i==alloc){
      alloc<<=1;
      ltisoid = (short       *)realloc(ltisoid, alloc*sizeof(short      ));
      ltgf    = (PREC_LNDATA *)realloc(ltgf,    alloc*sizeof(PREC_LNDATA));
      ltwl    = (PREC_LNDATA *)realloc(ltwl,    alloc*sizeof(PREC_LNDATA));
      ltelow  = (PREC_LNDATA *)realloc(ltelow,  alloc*sizeof(PREC_LNDATA));
      /* Check for allocation errors:      */
      if(!ltgf || !ltwl || !ltelow || !ltisoid)
        transiterror(TERR_CRITICAL|TERR_ALLOC,
                     "Cannot allocate memory for linetransition structure "
                     "to size %i in function readdatarng.\n", alloc);
    }
    /* If ASCII: skip comments and read the four fields:
       central wavelength, isotope ID, lowE, log(gf):    */
    if(li->asciiline){
      while((rc=fgetupto_err(lp=line, maxline, fp, &linetoolong, tr->f_line,
                             li->asciiline+offs++)) =='#' || rc=='\n');
      /* If it is not end of file, read the records:     */
      if(rc){
        /* Read values and report errors if necessary:   */
        ltwl[i] = strtod(lp, &lp2);
        if(lp==lp2) 
          return invalidfield(line, tr->f_line, li->asciiline+offs,
                              1, "central wavelength");
        ltisoid[i] = strtol(lp2, &lp, 0);
        if(lp==lp2)
          return invalidfield(line, tr->f_line, li->asciiline+offs,
                              2, "isotope ID");
        ltelow[i] = strtod(lp, &lp2);
        if(lp==lp2)
          return invalidfield(line, tr->f_line, li->asciiline+offs,
                              3, "lower energy level");
        ltgf[i] = strtod(lp2, &lp);
        if(lp==lp2)
          return invalidfield(line, tr->f_line, li->asciiline+offs,
                              4, "log(gf)");
      }
    }
    /* If binary read a data structure: */
    else{
      rn = fread(ltwl+i, sizeof(PREC_LNDATA), 1, fp);
      fread(ltisoid+i,   sizeof(short),       1, fp);
      fread(ltelow+i,    sizeof(PREC_LNDATA), 1, fp);
      fread(ltgf+i,      sizeof(PREC_LNDATA), 1, fp);
      transitDEBUG(26, verblevel, "Wavelength: %.8f iso: %i.\n", ltwl[i],
                   ltisoid[i]);
    }

    /* Warn if EOF was found.  This should be OK if you want to read
       until the very last line.        */
    if(!rn || !rc){
      transiterror(TERR_WARNING,
                   "End-of-file in datafile '%s'. Last wavelength read (%f) "
                   "was in record %i. If you are reading the whole range, "
                   "you can safely ignore this warning.\n",
                   tr->f_line, ltwl[i-1], i);
      break;
    }
    /* TD: Skip isotope at read time */

    /* End if we reached maximum requested wavelength: */
    if(ltwl[i]>finw)
      break;
    i++;
  }
  transitDEBUG(21, verblevel, "Number of lines just read: %li.\n", i);

  /* Realloc linetransition structure array to its final size and set the
     pointer in the transit structure: */
  alloc = i;
  /* Declare line_transition structure into lineinfo: */
  struct line_transition *lt = &li->lt;
  lt->isoid = (short       *)realloc(ltisoid, alloc*sizeof(short      ));
  lt->gf    = (PREC_LNDATA *)realloc(ltgf,    alloc*sizeof(PREC_LNDATA));
  lt->wl    = (PREC_LNDATA *)realloc(ltwl,    alloc*sizeof(PREC_LNDATA));
  lt->elow  = (PREC_LNDATA *)realloc(ltelow,  alloc*sizeof(PREC_LNDATA));
  /* Check for allocation errors:                       */
  if(!ltgf || !ltwl || !ltelow || !ltisoid){
    transiterror(TERR_CRITICAL|TERR_ALLOC,
                 "Cannot allocate memory for linetransition structure to "
                 "final size of %i in function readdatarng.\n", alloc);
    return -1;
  }

  /* Store number of lines: */
  li->n_l = i;

  fclose(fp);              /* Close file                      */
  tr->pi |= TRPI_READDATA; /* Update progress indicator       */
  return i;                /* Return the number of lines read */
}


/* TD: - Accept isotope list
       - Error return codes
       - Check that strtok doesn't overwrite first string */

/* \fcnfh
    Driver function to read TLI: read isotopes info, check margin
    and ranges, and read line transition information.

    Return: 0 on success.                                              */
int 
readlineinfo(struct transit *tr){
  long rn;  /* Sub-routines returned status */
  struct transithint *th=tr->ds.th; 
  static struct lineinfo li;
  static struct isotopes iso;
  memset(&li,  0, sizeof(struct lineinfo));
  memset(&iso, 0, sizeof(struct isotopes));
  tr->ds.li  = &li;  /* lineinfo */
  tr->ds.iso = &iso; /* isotopes */

  /* Read hinted info file: */
  transitprint(1, verblevel, "Reading info file '%s' ...\n", th->f_line);
  if((rn=readinfo_tli(tr, &li)) != 1)
    transiterror(TERR_SERIOUS, "readinfo_tli() returned an error "
                 "code %i.\n", rn);
  transitprint(1, verblevel, " Done.\n\n");

  /* Get the molecule index for the isotopes: */
  /* FINDME: Move this out of readline later. */
  rn = setimol(tr);
  rn = getisoratio(tr);

  /* Check the remainder (margin and range) of the hinted values
     related to line database reading:                           */
  if((rn=checkrange(tr, &li))<0)
    transiterror(TERR_SERIOUS,
                 "checkrange() returned error code %i!.\n", rn);
  /* Output status so far if the verbose level is enough:        */
  if(rn>0 && verblevel>1)
    transiterror(TERR_WARNING,
                 "checkrange() modified the suggested parameters, "
                 "it returned code 0x%x.\n\n", rn);

  /* Scale factors: */
  double fct = li.wavs.fct;
  double fct_to_microns = fct/1e-4;

  transitprint(2, verblevel,
               "After checking limits, the wavelength range to be "
               "used is %g to %g cm, including a margin of %g cm.\n",
               fct*tr->ds.li->wavs.i, fct*tr->ds.li->wavs.f, tr->margin);

  /* Read data file:            */
  transitprint(1, verblevel, "\nReading data ...\n");
  if((rn=readdatarng(tr, &li))<1)
    transiterror(TERR_SERIOUS,
                 "readdatarng() returned an error code %li\n", rn);
  transitprint(1, verblevel, "Done.\n\n");

  /* Status so far:             */
  transitprint(2, verblevel,
               "Status so far:\n"
               " * I read %li records from the datafile.\n"
               " * The wavelength range read was %.8g to %.8g microns.\n"
               " * Current margin is %.4g microns.\n"
               " * Usable range is thus %.8g to %.8g microns.\n",
               li.n_l,
               li.wavs.i*fct_to_microns, li.wavs.f*fct_to_microns,
               tr->margin*1e4,
               li.wavs.i*fct_to_microns + tr->margin*1e4,
               li.wavs.f*fct_to_microns - tr->margin*1e4);

/* FINDME: Why this indentation? */
#ifndef NODEBUG_TRANSIT
  rn = 1; /* A random index to test */
  struct line_transition *lt = &tr->ds.li->lt;
  transitDEBUG(21, verblevel,
               " * And the record %li has the following info\n"
               "Wavelength: %.10g\n"
               "Lower Energy Level: %.10g\n"
               "Log(gf): %.10g\n"
               "Isotope: %i\n",
               rn, lt->wl[rn], lt->elow[rn], lt->gf[rn], lt->isoid[rn]);
#endif

  transitDEBUG(21, verblevel,
               "Database min and max: %.10g(%.10g) and %.10g(%.10g)\n",
               li.wi, tr->ds.li->wi, li.wf, tr->ds.li->wf);
  return 0;
}


/* \fcnfh
   Frees lineinfo structure 
   Return: 0 on success                                      */
int
freemem_isotopes(struct isotopes *iso,
                 long *pi){
  int i;

  /* Free structures:                                         */
  for(i=0; i < iso->n_i; i++){      /* Allocated in readlineinfo */
    free_isof(iso->isof+i);
    free_isov(iso->isov+i);
  }
  for(i=0; i < iso->n_db; i++)
    free_db(iso->db+i);

  /* Free arrays:                                             */
  free(iso->isov);
  free(iso->isof);
  free(iso->db);
  free(iso->imol);
  free(iso->isoratio);

  /* Unset flags:                                             */
  *pi &= ~(TRPI_READINFO | TRPI_READDATA | TRPI_CHKRNG | TRPI_GETATM);
  return 0;
}



/* \fcnfh
   Free lineinfo structure.
   Return: 0 on success                */
int
freemem_lineinfotrans(struct lineinfo *li,
                      long *pi){
  int i;

  /* Free the four arrays of lt:        */
  struct line_transition *lt = &li->lt;
  free(lt->wl);
  free(lt->elow);
  free(lt->gf);
  free(lt->isoid);

  /* Free isov, dbnoext and samp in li: */
  free_isov(li->isov);
  free(li->isov);

  for(i=0; i<li->ndb; i++)
    free_dbnoext(li->db+i);
  free(li->db);

  free_samp(&li->wavs);

  /* Zero all the structure:            */
  memset(li, 0, sizeof(struct lineinfo));

  /* Unset appropiate flags:            */
  *pi &= ~(TRPI_READDATA | TRPI_READINFO | TRPI_CHKRNG);
  return 0;
}


/* \fcnfh
   Saves line information */
void
saveline(FILE *fp,
         struct lineinfo *li){
}


#ifdef DBGREADLINEINFO
/* \fcnfh
   main function for debugging only */
int main(int argc, char **argv){
  struct transit tr;
  struct transithint th;
  struct lineinfo *li;
  int i,ti1,nbins,ans;
  PREC_LNDATA *ltgf;
  PREC_LNDATA *ltelow;
  PREC_LNDATA *ltwl,*twl;
  short *ltisoid,*tisoid;

  tr.ds.th = &th;
  th.na = 0;

  verblevel = 20;

  th.m = 0.001;
  th.na |= TRH_WM;
  char defile_line[] = "./res/lineread.tli";
  th.f_line = (char *)calloc(strlen(defile_line)+1, sizeof(char));
  strcpy(th.f_line, defile_line);

  th.na |= TRH_FL;

  nbins=20;
  Pprintf(2, "Number of bins[%i]?: ", nbins);
  if(Pgeti(0, &ti1, 6)>0)
    nbins = ti1;

  if((i=readlineinfo(&tr))!=0)
    transiterror(TERR_CRITICAL, "Error code: %i.\n", i);
  transitDEBUG(20, verblevel, "range: %.10g to %.10g.\n",
               tr.ds.li->wi, tr.ds.li->wf);
  li = tr.ds.li;
  ltgf = tr.ds->lt.gf;
  ltwl = tr.ds->lt.wl;
  ltisoid = tr.ds->lt.isoid;
  ltelow = tr.ds->lt.elow;

  ti1 = (int)(log10(li->n_l)+1);

  printf("Done reading the file.\n\n"
         "dbread_pands() test results:\n");
  printf("Chosen wavelength range was from %.10g to %.2f [nm].\n"
         " %*li lines read.\n"
         " Choosing %i equal-sized bins, the result is\n",
         li->wi, li->wf, ti1, li->n_l, nbins);

  long qb[tr.n_i];
  float szb = (li->wf-li->wi)/nbins;
  double endb;
 
  twl = ltwl;
  tisoid = ltisoid;
  if(!nbins)
    /* FINDME: is this a typo? */
    Pprintf(1, "  hmmm, you chose 0 bins!.\n");
  for(i=0; i<nbins; i++){
    memset(qb, 0, sizeof(*qb)*4);
    endb = li->wi+(i+1)*szb;
    //    PPprintf(1,2,"KK %g %f\n",lp->wl,endb);
    while(*twl<endb && twl-ltwl<li->n_l){
      qb[*tisoid++]++;
      twl++;
    }

    Pprintf(1, " %*i = %i + %i + %i + %i lines shorter than %.3f\n",
            ti1, qb[0]+qb[1]+qb[2]+qb[3], qb[0], qb[1], qb[2], qb[3], endb);
  }

  Pprintf(1, "\nWanna know the value of a single record?\n"
             "If so, write record number (range 0 - %i), else "
             "press ^C: ", li->n_l-1);

  while(Pgeti(0,&ans,(int)(log10(li->n_l))+1)>=0){
    if(ans<li->n_l&&ans>=0){
      Pprintf(1, "Wavelength: %.10g\n", ltwl[ans]);
      Pprintf(1, "Lower Energy Level: %.10g\nLog(gf): %.10g\n",
              ltelow[ans], ltgf[ans]);
      printf("Isotope Name: %s\n", tr.isof[ltisoid[ans]].n);
    }
    else
      Pprintf(1, "\nInvalid record number, so ...");

    Pprintf(1, "\nWanna know the value of another single record?\n"
               "If so, write the record number (range 0 - %i), else just "
               "press ^C: ", li->n_l-1);
  }
  
}

#undef checkprepost

#endif

