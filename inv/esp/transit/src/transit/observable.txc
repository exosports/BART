/*
 * observable.c
 * observable.txc - Finds the optical depth. Component of the Transit program.
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
   Obtains the quantity that is observable, but before being convolved
   to telescope resolution

   @returns 0 on success
 */
int
modulation(struct transit *tr)	/* Main structure */
{
  long w;
  prop_samp *ip=&tr->ips;
  prop_samp *wn=&tr->wns;

  transitcheckcalled(tr->pi,"modulation",4,
		     "tau",TRPI_TAU,
		     "makeipsample",TRPI_MAKEIP,
		     "makewnsample",TRPI_MAKEWN,
		     "setgeom",TRPI_GEOMETRY
		     );

  //output and geometry variables.
  PREC_RES *out=tr->outpret=(PREC_RES *)calloc(wn->n,sizeof(PREC_RES));
  struct geometry *sg=tr->ds.sg;
  struct optdepth *tau=tr->ds.tau;

  PREC_RES *t;

  setgeom(sg,HUGE_VAL,&tr->pi);

  gsl_interp_accel *acc=gsl_interp_accel_alloc();
  for(w=0;w<wn->n;w++){
    out[w]=tr->sol->obsperwn(tau->t[w],tau->first[w],tau->toomuch,
			     ip->v,ip->n,sg,acc);
  }
  gsl_interp_accel_free(acc);

  return 0;
}
