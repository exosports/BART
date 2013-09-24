/*
 * idxrefraction.c - Finds the index of refraction. Component of the 
 *                   Transit program.
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

#include <transit.h>

/* \fcnfh
   Calculates the index of refraction. Right now it only gives 1 to all
   levels.

   @returns 0 on success
 */
int
idxrefrac(struct transit *tr)
{
  static struct idxref st_idx;
  long r;

  transitcheckcalled(tr->pi,"idxrefrac",1,
		     "makeradsample",TRPI_MAKERAD
		     );

  tr->ds.ir = &st_idx;
  prop_atm *atm = &tr->atm;

  PREC_ATM nustp = 0; 		/* TD: Allow for ray bending. Tau2 has
				   to be enabled as well */
  PREC_ATM rho;

  //allocate space and initialize
  st_idx.n=(PREC_RES *)calloc(tr->rads.n,sizeof(PREC_RES));

  for(r=0;r<tr->rads.n;r++){
    rho = stateeqnford(1, 1.0, atm->mm[r], 0, atm->p[r], atm->t[r]);
    st_idx.n[r] = 1 + rho*nustp/(LO*AMU*atm->mm[r]);
  }

  //set progress indicator and return success
  tr->pi|=TRPI_IDXREFRAC;
  return 0;
}


/* \fcnfh
   Fees memory for index of refraction structure

   @returns 0 on success;
 */
int
freemem_idexrefrac(struct idxref *ir,
		   long *pi)
{
  //free arrays
  free(ir->n);

  //clear progress indicator and return success
  *pi&=~(TRPI_IDXREFRAC|TRPI_TAU);
  return 0;
}


/* \fcnfh
   Restore hints structure, the structure needs to have been allocated
   before

   @returns 0 on success
            -1 if not all the expected information is read
	    -2 if info read is wrong
	    -3 if cannot allocate memory
	    1 if information read was suspicious
*/
int
restidxref(FILE *in,
	   PREC_NREC nrad,
	   struct idxref *ir)
{

  if(nrad<0)
    return -2;
  if(nrad>1000000)
    return 1;
  if((ir->n=(PREC_RES *)calloc(nrad,sizeof(PREC_RES)))==NULL)
    return -3;
  if(nrad==0)
    return 0;
  if(fread(ir->n,sizeof(PREC_RES),nrad,in)!=nrad)
    return -1;

  return 0;
}


/* \fcnfh
   Saves index of refraction structure
*/
void
saveidxref(FILE *out,
	   PREC_NREC nrad,
	   struct idxref *ir)
{
  if(nrad>0)
    fwrite(ir->n,sizeof(PREC_RES),nrad,out);
}
