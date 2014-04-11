/*
 * at_file.c - Read atmospheric info from file. Component of the
 *             Transit program.
 *
 * Copyright (C) 2004-2006 Patricio Rojo (pato@astro.cornell.edu)
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
 */

//#include "at_common.c"
#include <readatm.h>

#define ROUNDOFF 1e7
#define ROUNDTHRESH 1e-5

static double zerorad=0;

char *atmfilename;
_Bool *isolineinatm;

/* FINDME, Comment me: */
struct fonly {
  char n[maxeisoname];
  PREC_ATM q;
} *fonly;

static int nfonly=0;


/* \fcnfh
   Store info about the atmopshere file */
void
storename(struct atm_data *at,
          char *line){
  /* Get the length of the line char array: */
  int len;
  while(*line==' ' || *line=='\t')
    line++;
  len = strlen(line);

  /* Allocate and store info in atm_data.info:            */
  if(!at->info){ /* Only if it hasn't been stored before: */
    at->info = calloc(len+1, sizeof(char));
    strcpy(at->info, line);
  }
}


/* \fcnfh
   Find the isotope in 'isof'  with name 'iso' and get its
   abundance at radius given by index 'r'

   Return: abundance at requested radius                 */
static inline double
findfactq(char *iso,       /* Name of isotope searched   */
          prop_isof *isof, /* Fixed isotope info         */
          prop_isov *isov, /* Variable isotope info      */
          int n,           /* Number of isotopes         */
          PREC_NREC r){    /* Radius index               */

  /* Search molecule: */
  while (--n)
    if(strcasecmp(iso, isof[n].n)==0)
      //return isov[n].q[r];
      return 0;
  /* n == 1 case:     */
  if(strcasecmp(iso, isof->n)==0)
    //return isov->q[r];
    return 0;

  /* If it wasn't in isof, then search in fonly: */
  for(n=0; n<nfonly; n++)
    if(strcasecmp(iso, fonly[n].n)==0)
      return fonly[n].q;

  transiterror(TERR_SERIOUS,
               "Isotope you want to reference(%s) was not found among "
               "those whose abundance was given.\n", iso);
  return -1;
}


/* \fcnfh
   Add the (fractional) abundance from all isotopes except for the
  'other'-factor isotopes

   Return: total abundance                                               */
static inline double
addq(prop_isov *isov,   /* Variable isotope info (abundance among other) */
     enum isodo *isodo, /* Just add those that are ignored               */
     _Bool *other,      /* Don't add those that are going to take care
                           of the remainder to unity                     */ 
     int n,             /* Number of isotopes                            */
     PREC_NREC r){      /* Radius at which to compute everything         */
  double res = 0;  /* Sum of abundances */

  while(--n)
    if(!other[n])
      //res += isov[n].q[r];
      res += 0;
  if(!*other)
    //res += isov->q[r];
    res += 0;

  if(res>1.001 || res<0)
    transiterror(TERR_SERIOUS,
                 "Without processing 'other' molecules, abundance "
                 "addition(%g) is either bigger than 1 or negative!\n", res);
  return res;
}


/* \fcnfh
   Print error message when a line of 'file' is longer than 'max' characters */
static void
atmerr(int max,     /* Maximum length of an accepted line */
       char *file,  /* File from which we were reading     */
       int line){   /* Line being read                     */
  transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
               "Line %i of file '%s' has more than %i characters, "
               "that is not allowed\n", file, max);
  exit(EXIT_FAILURE);
}


/* \fcnfh
   Print error message when a field with transition info is invalid */
static void
invalidfield(char *line,   /* Contents of the line */
             int nmb,      /* File number          */
             int fld,      /* Field with the error */
             char *fldn){  /* Name of the field    */
  transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
               "Line %i of file '%s': Field %i (%s) does not have a valid "
               "value:\n%s.\n", nmb, atmfilename, fld, fldn, line);
  exit(EXIT_FAILURE);
}


/* \fcnfh
   Check that val is positive. Throw error message if not. */
static inline void 
checkposvalue(PREC_RES val, /* Value to check              */
              int field,    /* Field where it was read     */
              long line){   /* Line from which it was read */

  if(val<0)
    transiterror(TERR_SERIOUS,
                 "While reading the %i-th field in line %li of atmosphere "
                 "file %s, a negative value was found (%g)\n", field, 
                 line-1, atmfilename, val);
}


/* \fcnfh
    Get keyword variables from atmosphere file (mass/number abundance bool;
    zero-radius offset; radius, temperature, and pressure units factor;
    atmfile name/info; list isotopes; list of proportional-abundance isotopes).
    Store molecules and proportional isotopes in atm_data struct. 
    Determine which linedb isotope corresponds to such atm_data isotope.
    Solve non-matched linedb isotope cases.
    Put all non-ignore isotopes in transit.ds.iso structure.

    Return: Number of lines read                                     */
int
getmnfromfile(FILE *fp,                /* Pointer to atmospheric file    */
              struct atm_data *at,     /* atmosphere structure           */
              struct transit *tr,      /* transit structure              */
              PREC_ZREC *f_remainder){ /* Remainder molecules' factor    */
  struct molecules *mol=tr->ds.mol;
  char line[maxline], *lp;
  int nimol=0, /* Number of molecules with abundance profile */
      nmol=0,  /* Total number of molecules                  */
      i;       /* Auxiliary for-loop index                   */
  double cumulother = 0; /* Cumulative remainder-molecules' factor */
  int ipi = 0;    /* Number of remainder molecules   */

  /* Is the isotope defined in the atm file?: */
  //isoprop = (struct atm_isoprop *)calloc(ipa, sizeof(struct atm_isoprop));

  at->begline = 0; /* Line where the info begins      */

  /* Read and store the keyword atmospheric variables: */ 
  while(1){
    switch(fgetupto_err(line, maxline, fp, &atmerr, atmfilename,
                        at->begline++)){
    /* Ignore comments and blank lines: */
    case '\n':
    case '#':
      continue;
    case 0:     /* Throw error if EOF   */
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                   "readatm :: EOF unexpectedly found at line %i "
                   "of file %s while no t,p data points have been read.\n",
                   at->begline, atmfilename);
      exit(EXIT_FAILURE);
      continue;

    /* Determine whether abundance is by mass or number:     */
    case 'q':  
      lp = line + 1;
      while(*lp++ == ' '); /* Skip blank spaces              */
      lp--;
      switch(*lp|0x20){
      case 'n':
        at->mass = 0;  /* Number abundance (mixing ratio)    */
        break;
      case 'm':
        at->mass = 1;  /* Mass abundance (mass mixing ratio) */
        break;
      default:
        transiterror(TERR_SERIOUS,
                     "'q' option in the atmosphere file can only be followed "
                     "by 'm' (for abundances by mass) or 'n' (for abundances "
                     "by number). '%s' is invalid.\n", line);
        break;
      }
      continue;

    /* Zero radius value: */
    case 'z':  
      zerorad = atof(line+1);
      continue;

    /* Radius, temperature, or pressure units factor: */
    case 'u':
      switch(line[1]){
      case 'r':
        at->rads.fct = atof(line+2);
        break;
      case 'p':
        at->atm.pfct = atof(line+2);
        break;
      case 't':
        at->atm.tfct = atof(line+2);
        break;
      default:
        transiterror(TERR_SERIOUS, "Invalid unit factor indication in "
                                   "atmosphere file.\n");
        exit(EXIT_FAILURE);
      }
      continue;

    case 'n':  /* Name or identifier for file data */
      storename(at, line+1);
      continue;

    case 'i':  /* Molecule names with an abundance profile: */
      /* Count the number of wrds (molecules) in line:      */
      nimol = countfields(line+1, ' ');
      transitprint(15, verblevel, "The number of molecules is %d.\n", nimol);

      /* Allocate Molecules names:                          */
      mol->name    = (char **)calloc(nimol,             sizeof(char *));
      mol->name[0] = (char  *)calloc(nimol*maxeisoname, sizeof(char));
      for(i=1; i<nimol; i++)
        mol->name[i] = mol->name[0] + i*maxeisoname;

      transitprint(1, verblevel, "Molecules with abundance profile:\n  ");
      lp = line;
      lp = nextfield(lp); /* Skip keyword                   */
      /* Read and store names:                              */
      for (i=0; i<nimol; i++){
        getname(lp, mol->name[i]);
        lp = nextfield(lp);
        transitprint(1, verblevel, "%s, ", mol->name[i]);
      }
      transitprint(1, verblevel, "\b\b.\n");
      continue;

    /* Molecules with abundance proportional to the remainder: */
    case 'f':
      lp = line;
      lp = nextfield(lp); /* Skip keyword                      */

      /* Current total number of molecules:                    */
      nmol = ++ipi + nimol;
      /* Re-allocate to add the new molecule:                  */
      mol->name    = (char **)realloc(mol->name, nmol*sizeof(char *));
      mol->name[0] = (char  *)realloc(mol->name[0],
                                                 nmol*maxeisoname*sizeof(char));
      for (i=1; i<nmol; i++)
        mol->name[i] = mol->name[0] + i*maxeisoname;

      /* Re-allocate remainder factors:                        */
      f_remainder = (PREC_ZREC *)realloc(f_remainder, ipi*sizeof(PREC_ZREC));

      /* Read and store the molecule's name:                   */
      getname(lp, mol->name[nmol-1]);

      lp = nextfield(lp);   /* Move pointer to next field      */
      if(*lp == '=')        /* Skip an optional equal '=' sign */
        lp++;

      /* Read and store factor:                                */
      f_remainder[ipi-1] = strtod(lp, NULL);
      transitprint(30, verblevel, "%s remainder factor: %.3f\n",
                                  mol->name[nmol-1], f_remainder[ipi-1]);
      if(f_remainder[ipi-1] < 0)
        transiterror(TERR_CRITICAL,
                     "Abundance ratio has to be positive in atmosphere "
                     "file '%s' in line: '%s'.\n", atmfilename, line);
      continue;

    /* End of keyword variables: */
    default:   
      break;
    }
    break;
  }
  transitprint(1, verblevel, "Molecules with abundance proportional to "
                             "remainder:\n  ");
  for(i=nimol; i<nmol; i++)
    transitprint(1, verblevel, "%s, ", mol->name[i]);
  transitprint(1, verblevel, "\b\b.\n");

  transitprint(3, verblevel, "Read all keywords in atmosphere file without "
                             "problems.\n");

  /* Set total number of molecules in atmosphere: */
  mol->nmol = at->n_aiso = nmol;

  /* Check that there was at least one isotope defined and re-allocate 
     array sizes to their final size:                                */
  if(!nimol)
    transiterror(TERR_SERIOUS, "No isotopes were found in atmosphere file, "
                               "make sure to specify them in a line starting "
                               "with the letter 'i'. First non-comment line "
                               "read:\n%s\n", line);

  /* Set position of beginning of data: */
  at->begpos = ftell(fp) - strlen(line) - 1;

  /* Calculate cumulative fraction of remainder molecules: */
  for(i=0;  i < nmol-nimol;  i++)
    cumulother += f_remainder[i];

  transitprint(30, verblevel, "Cumulative remainder fraction: %.4f.\n",
                               cumulother);
  /* Check that cumulother sums to 1.0 (within allowed errors):  */
  if(nmol>nimol  &&  abs(1.0 - cumulother) > ROUNDTHRESH)
    transiterror(TERR_SERIOUS, "Sum of remainder-molecules fractional "
           "abundance (%g) must add to 1.0 +/- %g.\n", cumulother, ROUNDTHRESH);

  /* Resolve what to do with those isotopes that appear in the
     line transition database, but not in the atmosphere file. Get
     the number of non-ignored isotopes in atm_data without linelist: */
  //at->n_niso = checknonmatch(tr, at, isodo);
  /* FINDME: This will be a task in readline (if actually needed). */

  return at->begline;
}


/* \fcnfh
    Read radius, pressure, temperature, and abundances and store it into
    at_data of transit.  Calculate mean molecular mass and densities.

    Detailed:
    Read and store radius, pressure, and temperature from file.
    Read abundances for each (non other-factor) isotope.
    Sum fractional abundances. Calculate ramaining (other-factor) abundances.
    Calculate mean molecular mass per radius.
    Calculate densities per isotope at each radius.

    Returns: number of sample radius                                         */
int
readatmfile(FILE *fp,                /* Atmospheric file               */
            struct transit *tr,      /* transit struct                 */
            struct atm_data *at,     /* Atmosphere struct              */
            prop_samp *rads,         /* Radius sampling                */
            int nrad,                /* Size of allocated radius array */
            PREC_ZREC *f_remainder){ /* Remainder molecules' factor    */

  transitprint(1, verblevel, "Start reading abundances.\n");
  /* Find abundance related quantities for each radius */
  int lines = at->begline;
  PREC_NREC r = 0; /* Radius index (number of radii being read) */
  char rc;         /* File reading output */
  float allowq = 1 - tr->allowrq;
  int nabundances;  /* Number of abundances in list */
  double sumq;      /* Sum of abundances per line   */
  char line[maxline], *lp, *lp2;
  prop_mol *molec = at->molec;
  struct molecules *mol = tr->ds.mol;
  int i, j;            /* Auxiliary for-loop indices */
  /* Variables to be used by factor (except ieq which is general): */

  /* Count the number of abundances in each line:                      */
  fseek(fp, at->begpos, SEEK_SET); /* Go to position where data begins */
  /* Skip comments:                    */
  while((rc=fgetupto_err(lp=line, maxline, fp, &atmerr, atmfilename, lines++))
        =='#' || rc=='\n');
  /* Count values per line:            */
  nabundances = countfields(lp, ' ') - 3; /* Subtract rad, p, and T columns */
 
  fseek(fp, at->begpos, SEEK_SET); /* Go to position where data begins */
  while(1){
    /* Reallocate if necessary: */
    if(r==nrad){
      nrad <<= 1;
      rads->v     = (PREC_ATM *)realloc(rads->v,   nrad*sizeof(PREC_ATM));
      at->atm.t   = (PREC_ATM *)realloc(at->atm.t, nrad*sizeof(PREC_ATM));
      at->atm.p   = (PREC_ATM *)realloc(at->atm.p, nrad*sizeof(PREC_ATM));
      at->mm      = (double   *)realloc(at->mm,    nrad*sizeof(double));
      for(i=0; i<at->n_aiso; i++){
        molec[i].d = (PREC_ATM *)realloc(molec[i].d, nrad*sizeof(PREC_ATM));
        molec[i].q = (PREC_ATM *)realloc(molec[i].q, nrad*sizeof(PREC_ATM));
        molec[i].n = nrad;
      }
    }

    /* Skip comments and read next line: */
    while((rc=fgetupto_err(lp=line, maxline, fp, &atmerr, atmfilename, lines++))
          =='#' || rc=='\n');
    /* If it is end of file, stop loop: */
    if(!rc)
      break;

    /* Read and store radius, pressure, and temperature from file: */
    rads->v[r] = strtod(lp, &lp2) + zerorad; /* Radius       */
    checkposvalue(rads->v[r], 1, lines);       /* Check value is positive */
    if(lp==lp2) 
      invalidfield(line, lines, 1, "radius");
    at->atm.p[r] = strtod(lp2, &lp);         /* Pressure     */
    checkposvalue(at->atm.p[r], 2, lines); 
    if(lp==lp2)
      invalidfield(line, lines, 2, "pressure");
    at->atm.t[r] = strtod(lp, &lp2);         /* Temperature  */
    checkposvalue(at->atm.t[r], 3, lines);
    if(lp==lp2)
      invalidfield(line, lines, 3, "temperature");

    /* Read abundances for each isotope.  Keep reading-in values
       while there are numbers in line:                            */
    for(i=0, sumq=0; i<nabundances; i++){
      lp = lp2;
      /* Read the abundance of the isotope:                        */
      molec[i].q[r] = strtod(lp, &lp2);
      if (r==0)
        transitprint(30, verblevel, "density[%d, %li]: %.9f.\n",
                                    i, r, molec[i].q[r]);
      sumq += molec[i].q[r]; /* Add the abundances */
      checkposvalue(molec[i].q[r], i+4, lines); /* Check that tmp is positive */
      if(lp==lp2)
        invalidfield(line, lines, 4+i, "isotope abundance");
    }

    /* Remainder of the sum of abundances:     */
    /* Set abundance of remainder molecules:   */
    for(j=0; i < at->n_aiso; i++, j++)
      molec[i].q[r] = f_remainder[j]*(1-sumq);

    transitASSERT(i!=at->n_aiso, "The line %s of file %s contains %d abundance "
                                 "values, when there were %d expected.\n",
                                 __LINE__, __FILE__, i, at->n_aiso);
    
    /* Calculate mean molecular mass and check whether abundances add up
       to one (within roundoff error): */
    sumq = checkaddmm(at->mm+r, r, molec, mol, at->n_aiso, at->mass);
    if((int)(sumq*ROUNDOFF+0.5)<(int)(allowq*ROUNDOFF+0.5))
      transiterror(TERR_WARNING,
                   "In radius %g (%i: %g in file), abundances "
                   "don't add up to 1: %.9g\n",
                   at->rads.v[r], r, at->rads.v[r]-zerorad, sumq);

    /* Calculate densities using ideal gas law: */
    if (r==0){
      transitprint(30, verblevel, "Abund: %.9f, mmm: %.3f, mass: %.3f, "
                                "p: %.3f, T: %.3f.\n", molec[2].q[r], at->mm[r],
                                   mol->mass[2], at->atm.p[r]*at->atm.pfct,
                                   at->atm.t[r]*at->atm.tfct);
    }
    for(i=0; i<at->n_aiso; i++)
      molec[i].d[r] = stateeqnford(at->mass, molec[i].q[r], at->mm[r],
                                   mol->mass[i], at->atm.p[r]*at->atm.pfct,
                                   at->atm.t[r]*at->atm.tfct);
    transitprint(30, verblevel, "dens[%2li]: %.14f,   ", r, molec[2].d[r]);
    r++;
  }

  /* Re-allocate arrays to final size (nrad):  */
  rads->n = nrad = r;
  rads->v   = (PREC_ATM *)realloc(rads->v,   nrad*sizeof(PREC_ATM));
  at->atm.t = (PREC_ATM *)realloc(at->atm.t, nrad*sizeof(PREC_ATM));
  at->atm.p = (PREC_ATM *)realloc(at->atm.p, nrad*sizeof(PREC_ATM));
  at->mm    = (double   *)realloc(at->mm,    nrad*sizeof(double));
  for(i=0; i<at->n_aiso; i++){
    molec[i].d = (PREC_ATM *)realloc(molec[i].d, nrad*sizeof(PREC_ATM));
    molec[i].q = (PREC_ATM *)realloc(molec[i].q, nrad*sizeof(PREC_ATM));
    molec[i].n = nrad;
  }

  /* Free arrays that were used only to get the factorizing elements: */
  free(fonly);
  nfonly = 0;

  return nrad;
}


void getmass(struct atm_data *at, struct molecules *mol){
  int nmol = at->n_aiso;
  /* FINDME: De-hardcode filename, put it in tr.ds.at: */
  char *filename = "../inputs/molecules.dat";
  FILE *elist;

  /* Atomic masses, names, alias names, alias molecules, and sizes: */
  double *amass,    /* Atomic masses form list          */
         *radius;   /* Molecular radii from list        */
  char **aname,     /* Atomic symbol names              */
       **rname,     /* Molecules names for listed radii */
       **alias,     /* Alias of names given in atmfile  */
       **amol,      /* Corresponding molecule for alias */
       **elements;  /* Elements in a molecule           */
  int natoms = 92,  /* Number of listed atoms           */
      nalias =  2,  /* Number of listed alias names     */
      nradii = 14,  /* Number of listed radii           */
      namelen = 3,  /* Atomic symbol name length        */
      maxlinelen = 501,
      molnamelen,   /* Length of a molecule name        */
      elen,         /* Element name length              */
      iatom,        /* Atom's index from list           */
      ielement,     /* Element counter in a molecule    */
      i, j;         /* Auxiliary for-loop index         */
  int *nelements;   /* Number of elements in a molecule */

  char line[maxlinelen], *lp,
       molecule[MAXNAMELEN]; /* Current molecule's name */

  /* Alias names, and corresponding molecules: */
  amass    = (double *)calloc(natoms,         sizeof(double));
  aname    = (char  **)calloc(natoms,         sizeof(char *));
  aname[0] = (char   *)calloc(natoms*namelen, sizeof(char));
  for (i=1; i<natoms; i++)
    aname[i] = aname[0] + i*namelen;

  /* Open Molecules file: */
  if((elist=verbfileopen(filename, "Molecular info ")) == NULL)
    exit(EXIT_FAILURE);

  do{  /* Read lines, skipping comments and blank lines: */
    lp = fgets(line, maxlinelen, elist);
  }while (lp[0] == '\0' || lp[0] == '\n' || lp[0] == '#');

  /* Fill atoms and mass array with info from element list: */
  for (i=0; i<natoms; i++){
    lp += 19; /* The element's symbol starts at the 19th character */
    getname(lp, aname[i]);
    lp = nextfield(lp);
    amass[i] = strtod(lp, NULL);
    lp = fgets(line, maxlinelen, elist);
  }

  /* Allocate alias names and corresponding molecules: */
  alias    = (char  **)calloc(nalias, sizeof(char *));
  amol     = (char  **)calloc(nalias, sizeof(char *));
  alias[0] = (char   *)calloc(nalias*MAXNAMELEN, sizeof(char));
  amol[0]  = (char   *)calloc(nalias*MAXNAMELEN, sizeof(char));
  for (i=1; i<nalias; i++){
    alias[i] = alias[0] + i*MAXNAMELEN;
    amol[i]  = amol[0]  + i*MAXNAMELEN;
  }

  /* Continue reading the file to get the alias names: */
  do{  /* Skip blank and comment lines: */
    lp = fgets(line, maxlinelen, elist);
  }while (lp[0] == '\0' || lp[0] == '\n' || lp[0] == '#');

  /* Get aliases from file:         */
  for (i=0; i<nalias; i++){
    /* Get alias and molecule name: */
    getname(lp, alias[i]);
    lp = nextfield(lp);
    getname(lp, amol[i]);
    lp = fgets(line, maxlinelen, elist);
  }

  /* Allocate names and radii:         */
  radius   = (double *)calloc(nradii,            sizeof(double));
  rname    = (char  **)calloc(nradii,            sizeof(char *));
  rname[0] = (char   *)calloc(nradii*MAXNAMELEN, sizeof(char));
  for (i=1; i<nradii; i++)
    rname[i] = rname[0] + i*MAXNAMELEN;

  /* Go to next block                  */
  do{
    lp = fgets(line, maxlinelen, elist);
  }while (lp[0] == '\0' || lp[0] == '\n' || lp[0] == '#');

  /* Get radii from file:              */
  for (i=0; i<nradii; i++){
    /* Get molecules' name and radius: */
    getname(lp, rname[i]);
    lp = nextfield(lp);
    radius[i] = strtod(lp, NULL)/2.0;
    lp = fgets(line, maxlinelen, elist);
  }

  /* Allocate max number of molecules and max len of molecule name: */
  nelements   = (int   *)calloc(MAXNAMELEN, sizeof(int));
  elements    = (char **)calloc(MAXNAMELEN, sizeof(char *));
  elements[0] = (char  *)calloc((MAXNAMELEN+1)*MAXNAMELEN, sizeof(char));
  for (j=0; j<MAXNAMELEN; j++)
    elements[j] = elements[0] + j*(MAXNAMELEN+1);

  /* For each molecule: */
  for (i=0; i<nmol; i++){
    /* Check if molecule name is an alias: */
    if ((j=findstring(mol->name[i], alias, nalias)) >= 0)
      strcpy(molecule, amol[j]);
    else
      strcpy(molecule, mol->name[i]);

    /* Allocate elements in a molecule: */
    molnamelen  = (int)strlen(molecule);

    /* Break down molecule into its elements: */
    elen     = 0;  /* Element name length            */
    ielement = 0;  /* Elements in a molecule counter */
    for (j=0; j<molnamelen; j++){
      if (isalpha(molecule[j])){
        if (isupper(molecule[j])){  /* Uppercase letter: */
          if (elen > 0){
            /* End last element, advance j, store new letter */
            if (elen <= 2){ /* If name is longer, it's an alias name */
              elen = 0;
              ielement++;  /* Count 1 more element */
              nelements[ielement] = 1;
            }
            elements[ielement][elen++] = molecule[j];
          }
          else{ /* New Atom (elen==0) */
            elements[ielement][elen++] = molecule[j];
            nelements[ielement] = 1;
          }
        }
        else{ /* Lowercase: */
          elements[ielement][elen++] = molecule[j];
        }
      }
      else{  /* A numeric value: */
        nelements[ielement] = (int)strtol(&molecule[j], NULL, 10);
        elements[ielement][elen] = '\0';
        j += (int)log10((double)nelements[ielement++]);
        elen = 0;
      }
    }
    if (elen != 0)
     ielement++;

    /* Calculate molecule's mass: */
    for (j=0; j<ielement; j++){
      /* Find index of atom in list: */
      iatom = findstring(elements[j], aname, natoms);
      transitprint(30, verblevel, "Found %d %2s[%2d] atom(s) with mass "
                 "%9.6f u.\n", nelements[j], aname[iatom], iatom, amass[iatom]);
      /* Get mass and multiply by the number of atoms in molecule: */
      mol->mass[i] += amass[iatom] * nelements[j];
    }

    /* Set the radius: */
    j = findstring(molecule, rname, nradii);
    mol->radius[i] = radius[j] * ANGSTROM;
    transitprint(30, verblevel, "Molecule '%s' has radius %4.2f A and mass "
                      "%4.2f u.\n", mol->name[i], mol->radius[i]/ANGSTROM,
                       mol->mass[i]);
  }
}
