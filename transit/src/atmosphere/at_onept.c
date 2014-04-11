/*
 * at_onept.c - Take one point atmospheric info. Component of the
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

//#include "at_common.c"
#include <readatm.h>

/* \fcnfh
   ask from standard input for temperature, pressure and number of
   extraisotopes, it stores both in 'onept' and in 'at' at radius 'rad':  */
void
askonenpt(struct onept *onept,
          struct atm_data *at,
          int rad){  /* Radius at which to store values,
                        negative if this a one point run */
  _Bool one = 0;

  /* If this is a onepoint model: */
  if (rad<0){
    one = 1;
    rad = 0;
  }

  onept->nm = onept->nq = 0;

  /* Ask for number of extra isotopes if this is the zeroth radius only: */
  if(!rad)
    onept->ne = at->n_aiso = askforposl(" Number of extra isotopes for which "
                         "only abundance and molecular mass has to be given: ");

  /* Ask for pressure: */
  if(one)
    onept->p = askforposd(" Atmospheric pressure, cgs units (1e6cgs=1atm): ");
  else
    onept->p = askforposd(" Atmsopheric pressure for radius %i%s: ",
                          rad, rad?"":". cgs units (1e6cgs=1atm)");
  at->atm.p[rad] = onept->p;

  /* Ask for temperature: */
  if(one)
    onept->t = askforposd(" Atmospheric temperature, Kelvin degrees: ");
  else
    onept->t = askforposd(" Atmsopheric temperature for radius %i%s: ",
                          rad, rad?"":". Kelvin degrees");
  at->atm.t[rad] = onept->t;
}

/* \fcnfh
   Ask mass and name for new isotopes                                 */
void
askonemn(struct onept *onept, /* Input values are going to be stored here */
         prop_isof *isof,     /* Also stored here for use by the program  */
         int n,               /* Number of new isotopes                   */
         int nf){             /* Number of isotopes with full info        */

  int i;
  char rc;

  /* Allocate onept variables: */
  onept->m    = (PREC_ZREC *)calloc(n,             sizeof(PREC_ZREC));
  onept->n    = (char **)    calloc(n,             sizeof(char *));
  onept->n[0] = (char *)     calloc(n*maxeisoname, sizeof(char));
 
  /* For each isotope ask for mass and name: */
  for(i=0; i<n; i++){
    onept->n[i] = onept->n[0] + i*maxeisoname;

    while(1){
      fprintf(stderr, " Mass and name of extra isotope #%i (Order mandatory, "
                      "e.g. 12.011Carbon):\n  ", i+1);
      onept->m[i] = isof[i+nf].m = readds(stdin, &rc, onept->n[i],
                                          maxeisoname-1);
      strcpy(isof[i+nf].n, onept->n[i]);
      if(rc == 'q'){
        transitprint(0, verblevel, "User interrupt!\n");
        exit(EXIT_SUCCESS);
      }
      if(onept->m[i] <= 0)
        fprintf(stderr, " Invalid value %g, must be positive\n", onept->m[i]);
      else if(!rc)
        break;
      fprintf(stderr, "Try again!\n");
    }
  }
}



/* \fcnfh
   Set hard coded values  */
void
sethcdef(struct transit *tr,
         struct atm_data *at,
         prop_samp *rads){

  int i;
  int nmb = 0;
  //Hard coded values.
  //'hc\_t' temperature.
  //'hc\_meanmass' is both the mean mass and the mass of the only
  //isotope.
  //'hc\_abund' is the abundance of the first isotope, the rest are
  //just set at zero in the hardcode modality.
  PREC_ZREC hc_t       = 1350;
  PREC_ATM hc_pres     = 1.0e3;
  PREC_ATM hc_abund    = 6.5e-4;
  PREC_ATM hc_meanmass = 2.3;

  int nrad = rads->n;
  rads->v = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
  at->atm.tfct = 1;
  at->atm.pfct = 1;
  at->atm.t = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
  at->atm.p = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
  at->mm    = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
  rads->v[0] = 1.0;

  //at->n_niso = 0;
  struct isotopes *iso = tr->ds.iso;
  //nmb=iso->n_e = iso->n_i + at->n_niso;
  //at->isov = (prop_isov *)calloc(nmb, sizeof(prop_isov));

  at->atm.t[0] = hc_t;
  at->atm.p[0] = hc_pres;
  at->mm[0]    = hc_meanmass;
  for(i=0; i<nmb; i++){
    //at->isov[i].d = (PREC_ATM *)calloc(1, sizeof(PREC_ATM));
    //at->isov[i].d[0] = stateeqnford(at->mass, hc_abund, hc_meanmass,
    //                                hc_meanmass, hc_pres, hc_t);
    hc_abund = 0;
  }
  telldefaults(iso, at);
}
