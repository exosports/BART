/*
 * extraext.h   - Computes extra contributors to extinction. Component
 *                of the Transit program. 
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


/* \fcnfh
   Compute scatering contribution to extinction
*/
static __inline__ void
computeextscat(double *e, 
	       long n, 
	       struct extscat *sc,
	       double *temp,
	       double tcft,
	       double wn)
{
  long i;

  for(i=0;i<n;i++)
    e[i]=0;
}


/* \fcnfh
   Compute cloud contribution to extinction
*/
static __inline__ void
computeextcloud(double *e, 
	       long n,
	       struct extcloud *cl,
	       double *temp,
	       double tcft,
	       double wn)
{
  long i;

  for(i=0;i<n;i++)
    e[i]=0;
}
