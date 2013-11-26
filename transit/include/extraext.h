/*
 * extraext.h   - Computes extra contributors to extinction. Component
 *                of the Transit program. 
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
 */


/* \fcnfh
   Compute scatering contribution to extinction
*/
static __inline__ void
computeextscat(double *e, 
               long n, 
               struct extscat *sc,
               double *rad,
               double trad,
               double *temp,
               double tcft,
               double wn){
  long i;

  for(i=0; i<n; i++)
    e[i] = 0;
}


/* \fcnfh
   Compute cloud contribution to extinction
*/
static __inline__ void
computeextcloud(double *e, 
               long n,
               struct extcloud *cl,
               prop_samp *rad,
               double *temp,
               double tcft,
               double wn){
  long i;
  double *r = rad->v;
  double rfct = rad->fct;
  double rini = cl->rini*cl->rfct;
  double rfin = cl->rfin*cl->rfct;

  /* If there are no clouds, set array to zero: */
  if(rini==0){
    memset(e, 0, n*sizeof(double));
    return;
  }

  if(rad->d == 0)
    transiterror(TERR_SERIOUS,
                 "Radius needs to be equispaced for clouds prescription.\n");
  double slp = cl->maxe / (rfin - rini);

  /* Find radius index right below the cloud top layer: */
  for(i=n-1; i>=0; i--){
    if(r[i]*rfct <= rini)
      break;
    e[i] = 0; /* Zero extinction in this range          */
  }

  /* Set cloud extinction between cloud top and maxe:   */
  for(; i>=0; i--){
    if(r[i] * rfct <= rfin)
      break;
    e[i] = slp * (r[i]*rfct - rini);
  }

  /* Keep constant extinction maxe until the bottom:    */
  for(; i>=0; i--)
      e[i] = cl->maxe;
}
