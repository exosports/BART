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
  prop_samp *rad=&tr->rads;
  prop_samp *wn=&tr->wns;
  prop_samp *ip=&tr->ips;
  PREC_RES ***e=tr->ds.ex->e;

  //index, initial and final values
  long ii,wi;
  long inn,rnn,wnn;
  //'bb' is the impact parameter, while 'n' is the index of refraction,
  //'w' is the wavenumber, 't' is the tau as function of impact
  //parameter, 'r' is the radius
  PREC_RES *bb=ip->v;
  PREC_RES *n=tr->ds.ir->n;
  PREC_RES *r=rad->v;
  PREC_RES *t;

  transitcheckcalled(tr->pi,"tau",2,
		     "idxrefrac",TRPI_IDXREFRAC,
		     "extwn",TRPI_EXTWN
		     );

  transitacceptflag(tr->fl,tr->ds.th->fl,TRU_TAUBITS);

  //set tau structures' value
  tau.toomuch=50;
  if(tr->ds.th->na&TRH_TOOMUCH&&tr->ds.th->toomuch>0)
    tau.toomuch=tr->ds.th->toomuch;
  tau.iso=0;
  if(tr->ds.th->na&TRH_TAUISO&&tr->ds.th->tauiso>=0&&tr->ds.th->tauiso<tr->n_i)
    tau.iso=tr->ds.th->tauiso;
  tau.first=(long *)calloc(wn->n,sizeof(long));
  tau.t=(PREC_RES **)calloc(wn->n,sizeof(PREC_RES *));
  tau.t[0]=(PREC_RES *)calloc(wn->n*ip->n,sizeof(PREC_RES));
  for(ii=1;ii<wn->n;ii++)
    tau.t[ii]=tau.t[0]+ii*ip->n;

  //final index or number of elements
  wnn=wn->n;
  inn=ip->n;
  rnn=rad->n;

  //to temporarily store a per radius info, and the ratio of the ip and
  //rad units
  PREC_RES dt[wnn];
  double riw=ip->fct/rad->fct;

  //Need at least three radius to calculate a spline interpolation.
  if(inn<3)
    transiterror(TERR_SERIOUS,
		 "tau(): At least three impact parameters points are required!.\n"
		 );

  PREC_RES *out=tr->outpret=(PREC_RES *)calloc(wnn,sizeof(PREC_RES));
  struct geometry *sg=tr->ds.sg;

  gsl_interp_accel *acc=gsl_interp_accel_alloc();
  //for each wavenumber
  for(wi=0;wi<wnn;wi++){
    t=tau.t[wi];

    //For each resultant impact parameter
    for(ii=inn-1;ii>=0;ii--){
      if((t[ii]=rad->fct*tr->sol->tauperb(bb[ii]*riw,r,n,e,inn
					  ,tau.iso,wi,dt,acc))
	 >tau.toomuch){
	tau.first[wi]=ii;
	break;
      }
    }

    out[wi]=tr->sol->perwn(t,bb,inn,sg,acc);
  }
  gsl_interp_accel_free(acc);

  tr->pi|=TRPI_TAU;
  return 0;
}
