/*
 * tau.c
 * tau.txc - Finds the optical depth. Component of the Transit program.
 *
 * Copyright (C) 2003 Patricio Rojo (pato@astro.cornell.edu)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
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
   Calculates optical depth as a function of radii for a spherical
   symmetric planet

   @returns 0 on success
 */
int
tau(struct transit *tr)
{
  /* Different structures */
  static struct optdepth tau;
  tr->ds.tau=&tau;
  prop_isov *isov=tr->isov;
  prop_isof *isof=tr->isof;
  prop_samp *rad=&tr->rads;
  prop_samp *wn=&tr->wns;
  prop_samp *ip=&tr->ips;
  prop_atm *atm=&tr->atm;

  //index and final values
  long ii,ri,wi;
  long wnn,rnn;
  /* 'bb' is the impact parameter, while 'n' is the index of refraction,
     'w' is the wavenumber, 't' is the tau as function of impact
     parameter, 'r' is the radius */
  PREC_RES b,*bb=ip->v;
  PREC_RES *n=tr->ds.ir->n;
  PREC_RES *w=wn->v;
  PREC_RES *t,*r;
  PREC_RES r0;

  transitcheckcalled(tr->pi,"tau",2,
		     "idxrefrac",TRPI_IDXREFRAC,
		     "extwn",TRPI_EXTWN,
		     );

  //set tau structures' value
  if(tr->ds.th.na&TRH_TOOMUCH&&tr->ds.th.toomuch>0)
    tau.toomuch=tr->ds.th.toomuch;
  tau.first=(long *)calloc(wn->n,sizeof(long));
  tau.t=(PREC_RES **)calloc(wn->n,sizeof(PREC_RES *));
  tau.t[0]=(PREC_RES *)calloc(wn->n*ip->n,sizeof(PREC_RES));
  for(i=1;i<wn->n;i++)
    tau.t[i]=tau.t[0]+i*ip->n;

  //final index or number of elements
  wnn=wn->n;
  inn=ip->n-1;
  rnn=rad->n;

  //for each wavenumber
  for(wi=0;wi<wnn;wi++){
    t=tau.t[wi];

    //For each resultant impact parameter
    for(ii=inn;ii>=0;ii--){
      b=bb[ii];

      //Look for closest approach radius
      /* TD from here complete iteration and loops */
      r0=b;
      r0=b/interp(r0,r,n,rnn);
    }
  }

  tr->pi|=TRPI_TAU;
  return 0;
}
