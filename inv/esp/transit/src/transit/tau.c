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
   Calculates optical depth as a function of radii for a espherical
   symmetric planet

   @returns 0 on success
 */
int
tau(struct transit *tr)
{
  static struct optdepth tau;
  tr->ds.tau=&tau;
  prop_isov *isov=tr->isov;
  prop_isof *isof=tr->isof;
  prop_samp *rad=&tr->rads;
  prop_samp *wn=&tr->wns;
  prop_samp *ip=&tr->ips;
  prop_atm *atm=&tr->atm;

  long i;
  PREC_RES b;

  for(i=ip->n-1;i>=0;i--){
    b=ip->v[i];
  }

  tr->pi|=TRPI_TAU;
  return 0;
}
