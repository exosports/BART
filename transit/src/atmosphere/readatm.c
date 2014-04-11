/*
 * readatm.c - Read atmospheric info main file. Component of the
 *               Transit program.
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

/* List of functions:
int    getatm(struct transit *tr) 

double checkaddmm(double *mm, PREC_NREC r, prop_mol *molec,
                  struct molecules *mol, int n, _Bool mass)

void   telldefaults(struct isotopes *iso, struct atm_data *at)

void   saveonept_arr(FILE *out, struct onept *onept)

int    restonept_arr(FILE *in, struct onept *onept)

int    freemem_atmosphere(struct atm_data *at, long *pi)

void   freemem_onept(struct onept *o)                                         */

//#include "at_common.c"
#include <readatm.h>

/* \fcnfh
   Initialize ds.at (atm_data).  Set abundance mass and allowrq parameters.
   Check existence, open, and set pointer to atmosphere file.
   Get keyword variables from atm file (list of isotopes among others).
   Get temperature and isotopes abundances per radius from atm file.

   Return:
     0 on succes, elses
    -1 if no atmospheric info file was specified and no defaults are allowed
    -2 if default handling mode does not exist
    -3 if something really bad happened
    -4 if sum of abundances add up to more than 1
    -5 if requested isotope is an ignored one                */
int
getatm(struct transit *tr){
  int nmol,     /* Number of molecules      */
      nrad,     /* Number of radius samples */
      i;        /* for-loop index           */
  struct transithint *th = tr->ds.th; /* Get transithint                  */
  struct onept *onept = &th->onept;   /* Get onept from transit hint      */
  static struct atm_data at;       /* Declare atm_data structure       */
  static struct molecules mol;
  prop_samp *rads = &at.rads;      /* Radius sample                    */
  memset(&at,  0, sizeof(struct atm_data));  /* Set atm_data  mem to 0  */
  memset(&mol, 0, sizeof(struct molecules)); /* Set molecules mem to 0  */
  tr->ds.at  = &at;                 /* Set transit's atm_data structure */
  tr->ds.mol = &mol;

  PREC_ZREC *f_remainder;

  /* How the atmospheric parameters are input: */
  enum {invalid,       /* Unvalid input                   */
        given,         /* Given from one point            */
        ask,           /* Interactive from standard input */
        file}          /* From file                       */
        inp = invalid; /* Default */

  /* 'fp' and 'newiso' are initialized to avoid compiler output: */
  FILE *fp = NULL;

  /* Pass atmospheric flags into transit struct: */
  transitacceptflag(tr->fl, th->fl, TRU_ATMBITS); /* See transit.h */
  at.mass     = th->mass;  /* bool, abundance by mass (1) or by number (0) */
  tr->allowrq = th->allowrq;

  /* If --onept parameter was given, force a given one point: */
  if(onept->one){
    tr->f_atm = NULL;
    inp  = given;
    nrad = rads->n = 1;
    rads->fct      = 1;
  }
  /* If atmosphere filename was not set, use 1 point default: */
  else if(th->f_atm==NULL || strcmp(th->f_atm, "-")==0){
    tr->f_atm = th->f_atm;
    nrad = rads->n = 1;
    rads->fct      = 1;

    /* Set how to handle the default case: */
    switch(tr->fl & TRU_ATM1PBITS){
    case TRU_ATMNODEF:    /* throw error message:             */
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                   "getatm():: No atmospheric file specified.\n");
      return -1;
      break;
    case TRU_ATMHARDC1P:  /* Use hard-coded values:           */
      sethcdef(tr, &at, rads);
      return 0;
      break;
    case TRU_ATMASK1P:    /* Get values from standard input:  */
      inp = ask;
      break;
    default:              /* Undefined, throw error:          */
      transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
                   "getatm():: Unexistent default handling mode (0x%x)"
                   "requested.\n", tr->fl&TRU_ATM1PBITS);
      return -2;
      break;
    }
  }
  /* Else, atmospheric file was specified:           */
  else{
    /* Check that the file exists and can be opened. Set name and pointer
       in transit structure:                                              */
    if((tr->fp_atm=verbfileopen(th->f_atm, "Atmospheric info ")) == NULL)
      exit(EXIT_FAILURE);
    atmfilename = tr->f_atm = th->f_atm; /* Set file name  */
    /* FINDME: seems that there is no need to set th.f_atm */
    inp = file;
    fp  = tr->fp_atm; /* Pointer to file */
    transitprint(1, verblevel, "Reading atmosphere file: '%s'.\n", atmfilename);

    /* nrad will be the amount of allocated radii so far: */
    nrad      = 8;
    rads->n   = 0;
    rads->fct = 1.0;
  }

  /* Throw error if inp value is invalid: */
  transitassert(inp==invalid, "This shouldn't have happened, 'inp' variable "
                              "was not initialized in file %s, line %i.\n",
                              __FILE__, __LINE__);  /* See transit.h */

  /* Initialize atmosphere temperature-pressure arrays:    */
  at.atm.tfct = 1; /* Default temperature units are cgs */
  at.atm.pfct = 1; /* Default pressure    units are cgs */
  rads->v  = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
  at.atm.t = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
  at.atm.p = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
  rads->v[0] = 1.0;

  /* Set p, t, and newiso: */
  switch(inp){
  case given:
    //newiso = at.n_niso = onept->ne; /* Number of new isotopes, FINDME */
    at.atm.t[0] = onept->t;
    at.atm.p[0] = onept->p;
    break;
  case ask:
    /* FINDME: askonenpt ?? */
    askonenpt(onept, &at, -1);
    break;
  case file:
    /* T, P, and number of extra isotopes will be read from the file.
       Line beginning with 'i' will contain fields with mass and name,
       'newiso' will contain the number of currently allocated extra
       isotopes.                                                     */
    break;

  default: /* Just to avoid compiler output, it should never happen: */
    transiterror(TERR_CRITICAL, "This is not happening in file %s (%li)\n",
                 __FILE__, __LINE__);
    break;
  }

  /* Remainder molecules' abundance fraction: */
  f_remainder = (PREC_ZREC *)calloc(1, sizeof(PREC_ZREC));
  /* Read keyword-variables from file:        */
  if((i=getmnfromfile(fp, &at, tr, f_remainder))<1){
    /* On success getmnfromfile returns at.begline, not an error code */
    transiterror(TERR_SERIOUS, "getmnfromfile() returned error code %i\n", i);
    exit(EXIT_FAILURE);
  }
  nmol = at.n_aiso;

  /* Allocate mass and radius array of the molecules: */
  mol.nmol   = nmol;
  mol.mass   = (PREC_ZREC *)calloc(nmol, sizeof(PREC_ZREC));
  mol.radius = (PREC_ZREC *)calloc(nmol, sizeof(PREC_ZREC));
  mol.molec  = (prop_mol  *)calloc(nmol, sizeof(prop_mol));
  /* Get values from 'molecules.dat' file: */
  getmass(&at, &mol);

  /* Allocate arrays for the mean molecular mass, density, and abundance: */
  at.molec        = (prop_mol *)calloc(nmol, sizeof(prop_mol));
  at.mm           = (double   *)calloc(nrad, sizeof(double));
  for(i=0; i<nmol; i++){
    at.molec[i].d = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
    at.molec[i].q = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
    at.molec[i].n = nrad;
  }

  /* Get isotopic abundances: */
  nrad = readatmfile(fp, tr, &at, rads, nrad, f_remainder);
  transitprint(1, verblevel, "Done.\n\n");
  fclose(fp);

  /* Set required values in 'rads' structure: */
  rads->i = rads->v[0];
  rads->f = rads->v[rads->n-1];
  rads->o = 1;
  rads->d = 0;

  /* 'tauiso' is the isotope index which tau() will calculate afterwards.
      Check that it is between boundaries. */
  //tr->tauiso = 0;
  //if(th->tauiso >= 0  &&  th->tauiso < iso->n_i){
  //  if(iso->isodo[th->tauiso] == ignore  &&  tr->ds.ex->periso){
  //    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
  //                 "Selected isotope to compute tau (#%i: %s) was actually "
  //                 "ignored according to the atmospheric info\n",
  //                 th->tauiso, tr->ds.iso->isof[th->tauiso]);
  //    return -5;
  //  }
  //  else
  //    tr->tauiso = th->tauiso;
  //}

  /* Return succes and set progress indicator */
  tr->pi |= TRPI_GETATM;
  return 0;
}


/* \fcnfh
   Compute the mean molecular mass, check that sum of abundances
   is no larger than 1.

   Return: Sum of abundances.                                      */
double
checkaddmm(double *mm,          /* Mean molecular mass stored      */
           PREC_NREC r,         /* Radius position                 */
           prop_mol *molec,     /* Molecule info                   */
           struct molecules *mol, /* Molecules */
           int n,               /* Number of molecules             */
           _Bool mass){          /* Mass abundances?                */
  double sumq;  /* Fractional abundance sum */
  int i;        /* for index                */

  if(r >= molec[0].n)
    transiterror(TERR_CRITICAL,
                 "In file %s (line %li) a radius beyond the allocated "
                 "has been requested.", __FILE__, __LINE__);

  /* Compute the mean molecular mass: */
  sumq = *mm = 0;
  for(i=0; i<n; i++){
    if(mass)
      *mm += (molec[i].q[r])/(mol->mass[i]);
    else
      *mm += (molec[i].q[r])*(mol->mass[i]);
    sumq += molec[i].q[r];
  }
  if(mass)
    *mm = 1.0/(*mm);

  /* Check that sum of proportional abundances make sense: */
  if(sumq>1.001){
    transiterror(TERR_SERIOUS, "Sum of abundances of isotopes adds up to "
                               "more than 1: %g\n", sumq);
  }

  return sumq;
}


/* \fcnfh
   Print out default values when only one radius is being selected   */
void
telldefaults(struct isotopes *iso,
             struct atm_data *at){
  int i;
  transitprint(1, verblevel,
               "You are using one point atmospheric conditions:\n"
               " Temperature:         %g K\n"
               " Pressure:            %g dyne/cm2\n"
               " Mean molecular mass: %g AMU\n", 
               at->atm.t[0]*at->atm.tfct, at->atm.p[0]*at->atm.pfct, at->mm[0]);
  /* Densities for all isotopes: */
  for(i=0; i<iso->n_i; i++)
      transitprint(1, verblevel, " %-8s: density %8g g/cm3\n", 
                                   iso->isof[i].n,
                                   at->molec[iso->imol[i]].d[0]);
}


/* \fcnfh
   Save onept structure's array */
void
saveonept_arr(FILE *out,
              struct onept *onept){
  if(onept->ne <= 0)
    return;
  fwrite(onept->n[0],                 1, maxeisoname*onept->ne, out);
  fwrite(onept->m,    sizeof(PREC_ZREC), onept->ne,             out);
}


/* \fcnfh
   Restores onept structure's arrays
   @returns 0 on success, else:
           -1 if not all the expected information is read
           -2 if info read is wrong
           -3 if cannot allocate memory
            1 if information read was suspicious           */
int
restonept_arr(FILE *in,
              struct onept *onept){
  if(onept->ne <= 0)
    return 0;

  if((onept->n    = (char     **)calloc(onept->ne, sizeof(char    *)))==NULL ||
     (onept->n[0] = (char      *)calloc(onept->ne, sizeof(char     )))==NULL ||
     (onept->m    = (PREC_ZREC *)calloc(onept->ne, sizeof(PREC_ZREC)))==NULL  )
    return -3;

  int i;
  for(i=1; i<onept->ne; i++)
    onept->n[i] = onept->n[0]+i*maxeisoname;

  if(fread(onept->n[0], 1, maxeisoname*onept->ne, in) != maxeisoname*onept->ne)
    return -1;

  if(fread(onept->m, sizeof(PREC_ZREC), onept->ne, in) != onept->ne)
    return -1;

  return 0;
}


/* \fcnfh
   Free memory from atmosphere structure
   Return: 0 on success                           */
int
freemem_atmosphere(struct atm_data *at,
                   long *pi){
  /* Free structures:                             */
  free_samp(&at->rads);
  for(int i=0; i<at->n_aiso; i++)
    free_mol(at->molec+i);
  free_atm(&at->atm);

  /* Free arrays:                                 */
  free(at->mm);
  free(at->info);

  /* Clear progress indicator and return success: */
  *pi &= ~(TRPI_GETATM);
  return 0;
}


/* \fcnfh
   Free memory from onept structure
   Return: 0 on success                          */
void
freemem_onept(struct onept *o){
  free(o->q);
  if(o->nm){
    free(o->n[0]);
    free(o->n);
    free(o->m);
  }
}
