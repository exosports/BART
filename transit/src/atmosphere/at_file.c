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

#include "at_common.c"

#define ROUNDOFF 1e7

static double zerorad=0;

char *atmfilename;
_Bool *isolineinatm;
static struct atm_isoprop *isoprop;

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
   Check whether the isotope with 'name' corresponds to an isotope in
   lineinfofile. If yes, set's its isodo parameter.                     */
static void
isisoline(char *name,         /* Isotope name                           */
          PREC_ATM mass,      /* Isotope mass                           */
          int *isoeq,         /* Isotope lineinfo position to be stored */
          enum isodo atisodo, /* Action for current isotope             */
          prop_isof *isof,    /* Info from lineinfo file                */
          enum isodo *isodo,  /* Action to be taken storage             */
          PREC_NREC nliso){   /* Number of isotopes in the lineinfo     */
  PREC_NREC i;

  /* Note that here, even though the isotope might be ignored, it will
     still have an equivalent lineiso associated if one is found.
     \refline{isodbassoc}  */
  //transitprint(1, verblevel, " FLAG: isisoline 01.\n");
  //transitprint(1, verblevel, " FLAG: name is: %s.\n", name);
  //transitprint(1, verblevel, " FLAG: nliso: %li\n", nliso);
  for(i=0; i<nliso; i++){
    //transitprint(1, verblevel, " FLAG: isisoline: %s.\n", isof[i].n);
    if(strcasecmp(name, isof[i].n)==0){
      //transitprint(1, verblevel, " FLAG: isisoline 02.\n");
      *isoeq = i;
      //transitprint(1, verblevel, " FLAG: isisoline 03.\n");
      if(isolineinatm[i])
        transiterror(TERR_SERIOUS,
                     "Isotope %s has been defined more than once in the "
                     "atmosphere file.\n", name);
      isolineinatm[i] = 1;
      //transitprint(1, verblevel, " FLAG: isisoline 04.\n");
      if(isodo[i] != ignore)
        isodo[i] = atisodo;
      //transitprint(1, verblevel, " FLAG: isisoline 05.\n");
      if(isodo[i] == ignore)
        transitprint(2, verblevel, "Ignoring isotope %s (%g AMU).\n",
                     isof[i].n, isof[i].m);
      else if(isof[i].m!=mass && mass!=0){
        transiterror(TERR_WARNING,
                     "Mass of isotope %s, is not the same in the atmosphere "
                     "file %s(%g) than in the transition info file(%g).\n",
                     name, atmfilename, mass, isof[i].m);
        //transitprint(1, verblevel, " FLAG: isisoline 08.\n"); 
      }
      break;
    }
  }
  //transitprint(1, verblevel, " FLAG: isisoline 10.\n");
  /* Count non-matched and ignored atmosphere-file isotopes:  */
  if(*isoeq == -1  &&  atisodo == ignore){
    //transitprint(1, verblevel, " FLAG: isisoline 11.\n");
    nfonly++;
    //transitprint(1, verblevel, " FLAG: isisoline 12.\n");
  }
}


/* \fcnfh
   Ask what to do with lineinfo isotopes that don't have a match in 
   the atmosphere file

   Return: number of non-matched isotopes in atm_data struct        */
static int
checknonmatch(struct transit *tr,  /* Info about existent isotopes  */
              struct atm_data *at, /* Info about just read isotopes */
              enum isodo *isodo){  /* What the user will want with theisotope */
  int i, j,  /* Auxiliary for indices                              */
      rn,    /* Boolean to check that there are isotopes available
                to match with                                      */
      ison = at->n_aiso; /* Number of isotopes in atmospheric file */
  struct isotopes *iso = tr->ds.iso;

  /* For each of the non-matched isotopes, ask if they want to be
     matched, ignored or given a fixed value:                       */
  for(i=0; i<iso->n_i; i++){
    if(isodo[i]==unclear){
      /* If want to match then ask with what, if ignored or fixed it will
         be dealt with later */
      isodo[i] = askforposl("Isotope %s is not in atmosphere file, what do "
                            "you want to do? (1:match to some isotope, "
                            "2:ignore this isotope, 3:give a fixed abundance "
                            "value)\n", iso->isof[i].n);
      /* FINDME: option 3 is not implemented! */
      /* Invalid value: */
      if(isodo[i]>3){
        fprintf(stderr, "Invalid value, Try again.\n");
        isodo[i--] = 0; /* Set isodo to unclear and start loop again */
        continue;
      }
      /* Match to a molecule in atmfile: */
      if(isodo[i]==atmfile){
        while(1){
          rn = 0;
          /* Print all non-ignored isotopes from atmfile without a match: */
          for(j=0; j<ison; j++)
            if(at->isoeq[j] == -1  &&  at->isodo[j] != ignore){
              rn = 1;
              fprintf(stderr,"  %2i: %s (%gAMU)\n", j+1, at->n[j], at->m[j]);
            }
          /* Ask again if there were no free isotopes to choose from: */
          if(!rn){
            fprintf(stderr, "There are no isotopes to match with. Try "
                            "another option.\n");
            isodo[i--] = 0;
            continue;
          }
          /* Ask to match to one of them: */
          j = askforposd("Select a isotope number from the above list to "
                         "match %s with:", iso->isof[i].n);
          if(at->isoeq[j-1]==-1){
            at->isoeq[j-1] = i;
            break;
          }
          else
            fprintf(stderr, "Invalid value, Try again.\n");
        }
      }
    }
  }

  /* Count how many non-ignored isotopes from atmosphere file don't
      have a linefile:                                               */
  rn = 0;
  for(j=0; j<ison; j++)
    if( ( at->isoeq[j]==-1 && at->isodo[j]!=ignore ) || 
        ( at->isoeq[j]!=-1 && at->isodo[j]==factor   && 
          isoprop[at->isoeq[j]].eq==-1 ) )
      rn++;

  return rn;
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
      return isov[n].q[r];
  /* n == 1 case:     */
  if(strcasecmp(iso, isof->n)==0)
    return isov->q[r];

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
    if(isodo[n] != ignore && !other[n])
      res += isov[n].q[r];
  if(*isodo != ignore && !*other)
    res += isov->q[r];

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
getmnfromfile(FILE *fp,            /* Pointer to atmospheric file    */
              struct atm_data *at, /* atmosphere structure           */
              struct transit *tr,  /* transit structure              */
              int nmb){            /* Number of non-ignored isotopes */
  /* FINDME: nmb is unnecessary as argument */
  char line[maxline], *lp;
  int ison=0, /* Current number of molecules read */
      i;      /* Auxiliary for index              */
  struct isotopes *iso   = tr->ds.iso;
  enum   isodo    *isodo = iso->isodo;

  transitprint(1, verblevel, "FLAG: getmnfromfile 02.\n");
  /* Set variable to handle isotopic abundance porportions: */
  int ipi = 0,           /* Current number of proportional-isotopes read */
      ipa = at->ipa = 4; /* Number of elements in isoprop                */
  /* Is the isotope defined in the atm file?: */
  isolineinatm = (_Bool *)calloc(iso->n_i, sizeof(_Bool));
  isoprop = (struct atm_isoprop *)calloc(ipa, sizeof(struct atm_isoprop));

  transitprint(1, verblevel, "FLAG: getmnfromfile 05.\n");
  /* Allocate atm_data:                               */
  at->begline = 0; /* Line where the info begins      */
  enum isodo atisodo;
  at->isodo = (enum isodo *)calloc(nmb, sizeof(enum isodo));
  at->isoeq = (int *)       calloc(nmb, sizeof(int)); /* Isotope column index */
  at->m     = (PREC_ZREC *) calloc(nmb, sizeof(PREC_ZREC)); /* Isotope mass   */
  at->n     = (char **)     calloc(nmb, sizeof(char *));    /* Isotope name   */
  at->n[0]  = (char *)      calloc(nmb*maxeisoname, sizeof(char));
  /* Default initialization: */
  at->isoeq[0] = -1;
  for(i=1; i<nmb; i++){
    at->n[i] = at->n[0]+i*maxeisoname;
    at->isoeq[i] = -1;
  }

  //transitprint(1, verblevel, "FLAG: getmnfromfile 08.\n");
  /* While t,p data doesn't start, check for the various modifiers */
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

    /* Determine whether abundance is by mass or number: */
    case 'q':  
      //transitprint(1, verblevel, "FLAG: getmnfromfile 10a.\n");
      lp = line+1;
      while(*lp++==' '); /* Skip blank spaces */
      lp--;
      switch(*lp|0x20){
      case 'n':
        at->mass = 0;
        break;
      case 'm':
        at->mass = 1;
        break;
      default:
        transiterror(TERR_SERIOUS,
                     "'q' option in the atmosphere file can only be followed "
                     "by 'm' (for abundances by mass) or 'n' (for abundances "
                     "by number). '%s' is invalid.\n", line);
        break;
      }
      //transitprint(1, verblevel, "FLAG: getmnfromfile 10b.\n");
      continue;

    /* Zero radius value: */
    case 'z':  
      //transitprint(1, verblevel, "FLAG: getmnfromfile 11a.\n");
      zerorad = atof(line+1);
      //transitprint(1, verblevel, "FLAG: getmnfromfile 11b.\n");
      continue;

    /* Change factorization of radius, temperature, or pressure: */
    case 'u':
      //transitprint(1, verblevel, "FLAG: getmnfromfile 12a.\n");
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
        transiterror(TERR_SERIOUS,
                     "Invalid unit factor indication in atmosphere file.\n");
        exit(EXIT_FAILURE);
      }
      //transitprint(1, verblevel, "FLAG: getmnfromfile 12b.\n");
      continue;

    case 'n':  /* Name or identifier for file data */
      //transitprint(1, verblevel, "FLAG: getmnfromfile 13a.\n");
      storename(at, line+1);
      //transitprint(1, verblevel, "FLAG: getmnfromfile 13b.\n");
      continue;

    case 'i':  /* Isotope information */
      //transitprint(1, verblevel, "FLAG: getmnfromfile 14a.\n");
      lp = line+1;
      while(*lp==' ' || *lp=='\t') lp++;
      /* 'i' isotopes must be defined before 'f': */
      if(ipi)
        transiterror(TERR_CRITICAL,
                     "In line '%s': 'f' lines have to come after all the 'i' "
                     "lines in atmosphere file '%s'.\n", line, atmfilename);
      /* For each field: */
      //transitprint(1, verblevel, "FLAG: getmnfromfile 14b.\n");
      while(*lp){
        atisodo = atmfile;
        /* Allocate if necessary: */
        if(ison==nmb){
          //transitprint(1, verblevel, "FLAG: getmnfromfile 14c.\n");
          nmb <<= 1;
          at->isodo = (enum isodo *)realloc(at->isodo, nmb*sizeof(enum isodo));
          at->isoeq = (int *)       realloc(at->isoeq, nmb*sizeof(int));
          at->m     = (PREC_ZREC *) realloc(at->m,     nmb*sizeof(PREC_ZREC));
          at->n     = (char **)     realloc(at->n,     nmb*sizeof(char *));
          at->n[0]  = (char *)  realloc(at->n[0], nmb*maxeisoname*sizeof(char));
          for(i=1; i<nmb; i++)
            at->n[i] = at->n[0]+i*maxeisoname;
          for(i=nmb/2; i<nmb; i++)
            at->isoeq[i] = -1;
          //transitprint(1, verblevel, "FLAG: getmnfromfile 14d.\n");
        }

        /* Get mass and name, checking that is correct. First see if this
           isotope wants to be ignored. */
        //transitprint(1, verblevel, "FLAG: getmnfromfile 14e.\n");
        if(*lp == '!'){
          lp++;
          atisodo = ignore;
        }
        //transitprint(1, verblevel, "FLAG: getmnfromfile 14f.\n");
        /* Get mass and name: */
        at->m[ison] = getds(lp, 0, at->n[ison], maxeisoname-1);
        if(at->m[ison]<0 || at->n[ison]=='\0'){
          transiterror(TERR_SERIOUS,
                       "Invalid field in file %s, line %i, while reading "
                       "isotope info at:\n%s\n", atmfilename, at->begline, lp);
        }

        //transitprint(1, verblevel, "FLAG: getmnfromfile 14g.\n");
        /* Check if isotope is in the lineinfo file, set isoeq and isodo: */
        isisoline(at->n[ison], at->m[ison], at->isoeq+ison, atisodo,
                  iso->isof, isodo, iso->n_i);
        //transitprint(1, verblevel, "FLAG: getmnfromfile 14h.\n");
        /* Set atm isodo and increase ison counter by 1:                 */
        at->isodo[ison++] = atisodo;

        /* Skip over recently read field, and go to next field: */
        while(*lp!=' ' && *lp!='\0') lp++;
        while(*lp==' ' || *lp=='\t') lp++;
      }
      //transitprint(1, verblevel, "FLAG: getmnfromfile 14i.\n");
      continue;

    /* An isotope is to be taken as proportional to other */
    case 'f':
      //transitprint(1, verblevel, "FLAG: getmnfromfile 15a.\n");
      /* Move pointer 1 character and skip blank and tabs: */
      lp = line+1;
      while(*lp==' ' || *lp=='\t') lp++;
      /* Double size if ipi reached ipa: */
      if(ipi==ipa)
        isoprop=(struct atm_isoprop *)realloc(isoprop,
                                          (ipa<<=1)*sizeof(struct atm_isoprop));

      /* Get the isotope mass and name: */
      isoprop[ipi].m = getds(lp, 0, isoprop[ipi].n, maxeisoname-1);
      /* Skip the field just read and go to next one: */
      lp = nextfield(lp);
      /* Skip an optional equal '=' sign */
      if(*lp=='=' && lp[1]==' ')
        lp = nextfield(lp);

      /* Get factor, which has to be between 0 and 1 */
      isoprop[ipi].f = strtod(lp, NULL);
      if(isoprop[ipi].f < 0)
        transiterror(TERR_CRITICAL,
                     "Abundance ratio has to be positive in atmosphere "
                     "file '%s' in line: %s", atmfilename, line);
      lp = nextfield(lp);

      /* Get name of reference molecule: */
      i = 0;
      while(*lp)
        isoprop[ipi].t[i++] = *lp++;
      isoprop[ipi].t[i] = '\0';

      /* Check if the isotope is in the lineinfo file: */
      isoprop[ipi].eq = -1;
      /* Set eq, isodo and isolineinatm.  Check that mass in atm. file 
         match line info mass:                                          */
      isisoline(isoprop[ipi].n, isoprop[ipi].m, &isoprop[ipi].eq, factor,
                iso->isof, isodo, iso->n_i);

      /* Increase index and go for the next line: */
      ipi++;
      //transitprint(1, verblevel, "FLAG: getmnfromfile 15b.\n");
      continue;

    /* End of keyword variables: */
    default:   
      break;
    }
    break;
  }

  //transitprint(1, verblevel, "FLAG: getmnfromfile 20.\n");
  transitprint(3, verblevel,
               "Read all keywords in atmosphere file without problems.\n");

  /* Check that there was at least one isotope defined and re-allocate 
     array sizes to their final size:                                */
  if(!ison)
    transiterror(TERR_SERIOUS,
                 "No isotopes were found in atmosphere file, make sure to "
                 "specify them in a line starting with the letter 'i'. "
                 "First non-comment line read:\n%s\n", line);
  /* Set position of beginning of data: */
  at->begpos = ftell(fp)-strlen(line)-1;

  /* Ignored non-matched atmosphere-file isotopes: */
  fonly = (struct fonly *)calloc(nfonly, sizeof(struct fonly));
  /* Re-allocate isoprop array to its definitive size: */
  at->ipa = ipa = ipi;
  isoprop = (struct atm_isoprop *)realloc(isoprop,
                                          ipa*sizeof(struct atm_isoprop));
  /* Atmosphere arrays contain both, molecules and proportional isotopes: */
  nmb = at->n_aiso = ison + ipa;
  at->isodo = (enum isodo *)realloc(at->isodo, nmb*sizeof(enum isodo)      );
  at->isoeq = (int *)       realloc(at->isoeq, nmb*sizeof(int)             );
  at->m     = (PREC_ZREC *) realloc(at->m,     nmb*sizeof(PREC_ZREC)       );
  at->n     = (char **)     realloc(at->n,     nmb*sizeof(char *)          );
  at->n[0]  = (char *)      realloc(at->n[0],  nmb*maxeisoname*sizeof(char));
  for(i=1; i<nmb; i++)
    at->n[i] = at->n[0] + i*maxeisoname;

  /* Initialize values for the factorized elements: */
  for(i=ison; i<nmb; i++){
    strncpy(at->n[i], isoprop[i-ison].n, maxeisoname-1);
    at->n[i][maxeisoname-1] = '\0';
    at->isoeq[i]            = i-ison;  /* isotope index in isoprop */
    at->m[i]                = isoprop[i-ison].m;
    at->isodo[i]            = factor;
  }

  /* Resolve what to do with those isotopes that appear in the
     line transition database, but not in the atmosphere file. Get
     the number of non-ignored isotopes in atm_data without linelist: */
  at->n_niso = checknonmatch(tr, at, isodo);

  /* Re-alocate iso to add non-ignored non-linedb isotopes:             */
  nmb = iso->n_i + at->n_niso;  /* Total number of non-ignored isotopes */
  iso->isodo = (enum isodo *)realloc(iso->isodo, nmb*sizeof(enum isodo));
  iso->isof  = (prop_isof  *)realloc(iso->isof,  nmb*sizeof(prop_isof));
  iso->isof[iso->n_i].n = (char *)realloc(iso->isof[iso->n_i].n,
                                       (nmb-iso->n_i)*maxeisoname*sizeof(char));
  for(i=1; i<nmb-iso->n_i; i++)
    iso->isof[iso->n_i+i].n = iso->isof[iso->n_i].n + i*maxeisoname;

  /* Add non-linedb non-ignored atmosphere-file isotopes into iso struct: */
  nmb = iso->n_i;      /* Number of linedb isotopes                */
  ipi = 0;             /* Number of ignored isotopes not in linedb */
  int lineignore = 0;  /* Number of ignored isotopes in linedb     */

  /* Add to nmb the non-ignored non-linedb isotopes: */
  for(i=0; i<ison+ipa; i++) /* Over all atm isotopes */
    /* If the isotope is not associated to the linedb isotopes, then
       associate it. Note that isotopes in linedb that are ignored will
       be associated (see comments for Line \label{isodbassoc}). Hence
       they won't be detected in this IF.   */

    if(at->isoeq[i] == -1){ /* Not in linedb */
      /* If not ignored, add it to iso:      */
      if(at->isodo[i] != ignore){
        at->isoeq[i]     = nmb;
        iso->isodo[nmb]  = at->isodo[i];
        iso->isof[nmb].m = at->m[i];
        strcpy(iso->isof[nmb++].n, at->n[i]); /* Increase counter */
      }
      /* If ignored and non-linedb, must be a referenced molecule
         from another isotope:                                    */
      else{
        /* Set isoeq on fonly struct:             */
        at->isoeq[i] = ipi;
        /* Store in fonly and increase ipi count: */
        strcpy(fonly[ipi++].n, at->n[i]);
      }
    }

    /* Count linedb ignored isotopes: */
    else if(at->isodo[i] == ignore)
      lineignore++;

    /* If it is a factor non-linedb isotope: */
    else if(at->isodo[i] == factor && isoprop[at->isoeq[i]].eq == -1){
      /* FINDME: How can isodo be factor and ignore at the same time? */
      if(at->isodo[i]==ignore)
        transiterror(TERR_CRITICAL,
                     "Trying to ignore a factor isotope, that is not "
                     "posible.\n");
      isoprop[at->isoeq[i]].eq = nmb;
      iso->isodo[nmb]          = at->isodo[i];
      iso->isof[nmb].m         = at->m[i];
      strcpy(iso->isof[nmb++].n, at->n[i]);  /* Increase counter */
    }

  /* Reduce the array to get rid of non-linedb ignored isotopes: */
  iso->n_e = nmb;  /* Set number of extended isotopes */
  iso->isodo = (enum isodo *)realloc(iso->isodo, nmb*sizeof(enum isodo));
  iso->isof  = (prop_isof  *)realloc(iso->isof,  nmb*sizeof(prop_isof ));
  iso->isof[iso->n_i].n = (char *)realloc(iso->isof[iso->n_i].n,
                                          nmb*maxeisoname*sizeof(char));
  for(i=1; i<nmb-iso->n_i; i++)
    iso->isof[iso->n_i+i].n = iso->isof[iso->n_i].n + i*maxeisoname;

  /* Calculate cumulative fraction of 'other'-factor isotopes: */
  double cumulother = 0;
  for(i=0; i<at->n_aiso; i++)
    if(at->isodo[i]==factor){
      int feq = at->isoeq[i];
      if(strcasecmp(isoprop[feq].t, "other")==0)
        cumulother += isoprop[feq].f;
    }
  /* It doesn't make sense for cumulother to be anything different from
     unity (except round-off error): you want to associate the
     remainder of the atmosphere to some isotopic properties.            */
  if(cumulother!=0 && (int)(cumulother*ROUNDOFF+0.5)!=(int)(ROUNDOFF+0.5))
    transiterror(TERR_SERIOUS,
                 "Sum of fractional abundance of isotopes proportional to "
                 "'other' (%g) must add to 1.0.\n", cumulother);

  transitprint(1, verblevel, " FLAG: %i+%i+%i, %i+%i.\n", nfonly, lineignore,
                             nmb-lineignore, ison, ipa);
  transitASSERT(nmb+nfonly != ison+ipa,
                "Number of ignored-nonline elements (%i), plus the "
                "number of ignored-line elements (%i), plus the "
                "number of nonignored (%i), doesn't match the "
                "number of elements found in fields 'i' (%i) and 'f' (%i) "
                "of the atmosphere file '%s'.\n", nfonly, lineignore, 
                nmb-lineignore, ison, ipa, atmfilename);

  transitASSERT(nmb != iso->n_e,
                "Problem in file %s, line %i. Assertion failed: %i != %i!.\n",
                __FILE__, __LINE__, nmb, iso->n_e);

  /* Free unused array, store factor info in at structure and return line
     where T,P start */
  free(isolineinatm);
  at->isoprop = isoprop;

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
readatmfile(FILE *fp,            /* Atmospheric file  */
            struct transit *tr,  /* transit struct    */
            struct atm_data *at, /* Atmosphere struct */
            prop_samp *rads,     /* Radius sampling   */
            int nrad){           /* Number of allocated radii, note that
                                    is not returned updated */
  /* Find abundance related quantities for each radius */
  int lines = at->begline;
  PREC_NREC r = 0; /* Radius index (number of radii being read) */
  PREC_RES tmp;
  char rc;         /* File reading output */
  float allowq = 1-tr->allowrq;
  double sumq;
  char line[maxline], *lp, *lp2;
  prop_isov *isov = at->isov;
  int *isoeq      = at->isoeq;
  struct isotopes *iso = tr->ds.iso;
  enum isodo *isodo    = at->isodo;
  int i, neiso=iso->n_e;

  fseek(fp, at->begpos, SEEK_SET); /* Go to position where data begins */
  while(1){
    /* Reallocate if necessary: */
    if(r==nrad){
      nrad <<= 1;
      rads->v   = (PREC_ATM *)realloc(rads->v,   nrad*sizeof(PREC_ATM));
      at->atm.t = (PREC_ATM *)realloc(at->atm.t, nrad*sizeof(PREC_ATM));
      at->atm.p = (PREC_ATM *)realloc(at->atm.p, nrad*sizeof(PREC_ATM));
      at->mm    = (double   *)realloc(at->mm,    nrad*sizeof(double));
      for(i=0; i<neiso; i++){
        isov[i].d = (PREC_ATM *)realloc(isov[i].d, nrad*sizeof(PREC_ATM));
        isov[i].q = (PREC_ATM *)realloc(isov[i].q, nrad*sizeof(PREC_ATM));
        isov[i].n = nrad;
      }
    }

    /* Skip comments and read next line: */
    while((rc=fgetupto_err(lp=line, maxline, fp, &atmerr, atmfilename, lines++))
          =='#' || rc=='\n');
    /* If it is end of file, stop loop: */
    if(!rc)
      break;

    /* Read and store radius, pressure, and temperature from file: */
    tmp = rads->v[r] = strtod(lp, &lp2) + zerorad; /* Radius       */
    checkposvalue(tmp, 1, lines);     /* Check value is a positive */
    if(lp==lp2) 
      invalidfield(line, lines, 1, "radius");
    tmp = at->atm.p[r] = strtod(lp2, &lp);         /* Pressure     */
    checkposvalue(tmp, 2, lines); 
    if(lp==lp2)
      invalidfield(line, lines, 2, "pressure");
    tmp = at->atm.t[r] = strtod(lp, &lp2);         /* Temperature  */
    checkposvalue(tmp, 3, lines);
    if(lp==lp2)
      invalidfield(line, lines, 3, "temperature");

    /* Variables to be used by factor (except ieq which is general): */
    int ieq, /* Isotope's reference index */
        feq; /* Isoprop's reference index */
    double ref;
    _Bool otherfct[neiso]; /* is this isotope proportional to 'other' */
    memset(otherfct, 0, sizeof(otherfct));

    /* Read abundances for each isotope. Don't process factor isotopes
        because they might be proportional to a fixed element which is
        set below. */
    for(i=0; i<at->n_aiso; i++){
      ieq = isoeq[i];
      switch(isodo[i]){
      /* If fixed, ask for the abundance: */
      case fixed:
        if(!r){ /* Ask only in first radius iteration */
          isov[ieq].q[0] = askforposd("%s abundance for isotope %s: ",
                                    at->mass?"Mass":"Number", iso->isof[ieq].n);
          if(isov[ieq].q[0] >= 1){
            fprintf(stderr, "Abundance must be less than one. Try Again.\n");
            i--;
          }
        }
        else
          isov[ieq].q[r] = isov[ieq].q[0];
        break;
      case factor:
        feq = ieq;
        ieq = isoprop[feq].eq;
        /* Process 'other'-factor isotopes at the end: */
        if(strcasecmp(isoprop[feq].t, "other")==0){
          otherfct[ieq] = 1;
          continue;
        }
        /* Find q of the reference value:              */
        ref = findfactq(isoprop[feq].t, iso->isof, isov, neiso, r);
        isov[ieq].q[r] = isoprop[feq].f*ref;
        break;
      default:
        transiterror(TERR_CRITICAL,
                     "Trying to read isotope in readatmfile() which is "
                     "not 'fixed', 'atmfile', 'ignored', nor 'factor'.\n");
        exit(EXIT_FAILURE);
        break;
      case atmfile:
      case ignore:
        transitASSERT(ieq<0 || 
                      (isodo[i]==ignore && ieq>=nfonly) || 
                      (isodo[i]!=ignore && ieq>=iso->n_e),
                      "Assertion failed in file %s, line %i: %i!=[0,%i]. "
                      "fonly: %i\n", __FILE__, __LINE__, isoeq[i], 
                      isodo[i]==ignore?nfonly:iso->n_e-1, isodo[i]==ignore);
        /* Read the abundance of the isotope:        */
        if(isodo[i]==ignore)  /* factor-only isotope */
          tmp = fonly[ieq].q   = strtod(lp2, &lp);
        else
          tmp = isov[ieq].q[r] = strtod(lp2, &lp);
        checkposvalue(tmp, i+4, lines);
        if(lp==lp2)
          invalidfield(line, lines, 4+i, "isotope abundance");
        lp2 = lp;
        break;
      }
    }

    /* Remainder of the sum of abundances:       */
    ref = 1 - addq(isov, iso->isodo, otherfct, neiso, r);
    /* Set 'other'-factor isotopes' abundances:  */
    for(i=0; i<at->n_aiso; i++)
      if(isodo[i]==factor){
        feq = isoeq[i];
        ieq = isoprop[feq].eq;
        if(otherfct[ieq])
          isov[ieq].q[r] = isoprop[feq].f*ref;
      }

    /* Calculate mean molecular mass and check whether abundances add up
       to one (within roundoff error): */
    sumq = checkaddmm(at->mm+r, r, isov, iso->isof, neiso, at->mass,
                      iso->isodo);
    if((int)(sumq*ROUNDOFF+0.5)<(int)(allowq*ROUNDOFF+0.5))
      transiterror(TERR_WARNING,
                   "In radius %g (%i: %g in file), abundances "
                   "don't add up to 1: %.9g\n",
                   at->rads.v[r], r, at->rads.v[r]-zerorad, sumq);

    /* Calculate densities using ideal gas law: */
    for(i=0; i<neiso; i++)
      isov[i].d[r] = stateeqnford(at->mass, isov[i].q[r], at->mm[r],
                                  iso->isof[i].m, at->atm.p[r]*at->atm.pfct,
                                  at->atm.t[r]*at->atm.tfct);
    r++;
  }

  /* Re-allocate arrays to final size (nrad):  */
  rads->n = nrad = r;
  rads->v   = (PREC_ATM *)realloc(rads->v,   nrad*sizeof(PREC_ATM));
  at->atm.t = (PREC_ATM *)realloc(at->atm.t, nrad*sizeof(PREC_ATM));
  at->atm.p = (PREC_ATM *)realloc(at->atm.p, nrad*sizeof(PREC_ATM));
  at->mm    = (double   *)realloc(at->mm,    nrad*sizeof(double));
  for(i=0; i<neiso; i++){
    isov[i].d = (PREC_ATM *)realloc(isov[i].d, nrad*sizeof(PREC_ATM));
    isov[i].q = (PREC_ATM *)realloc(isov[i].q, nrad*sizeof(PREC_ATM));
    isov[i].n = nrad;
  }

  /* Free arrays that were used only to get the factorizing elements: */
  free(fonly);
  nfonly = 0;

  return nrad;
}


