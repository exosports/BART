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
            -1 if impact parameter sampling is not equispaced
 */
int
modulation(struct transit *tr)	/* Main structure */
{
  transitcheckcalled(tr->pi,"modulation",4,
		     "tau",TRPI_TAU,
		     "makeipsample",TRPI_MAKEIP,
		     "makewnsample",TRPI_MAKEWN
		     );

  //initial variables and check that impact parameters was a monospaced
  //array. Stop otherwise.
  long w;
  prop_samp *ip=&tr->ips;
  prop_samp *wn=&tr->wns;
  transit_ray_solution *sol=tr->sol;
  if(ip->d<=0&&sol->monoip){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "To compute %s modulation, the impact parameter has to\n"
		 "be an equispaced array\n"
		 ,sol->name);
    return -1;
  }

  //output and geometry variables.
  PREC_RES *out=tr->outpret=(PREC_RES *)calloc(wn->n,sizeof(PREC_RES));
  struct geometry *sg=tr->ds.sg;
  struct optdepth *tau=tr->ds.tau;

  //set time to the user hinted default
  setgeom(sg,HUGE_VAL,&tr->pi);


  transitprint(1,verblevel,
	       "Integrating for each wavelength. For the current range,\n"
	       "expect %li dots below...\n"
	       ,wn->n/512);
  //integrate for each wavelength
  gsl_interp_accel *acc=gsl_interp_accel_alloc();
  for(w=0;w<wn->n;w++){
    out[w]=sol->obsperwn(tau->t[w],tau->first[w],tau->toomuch,
			 ip,sg,acc);

    if((w&0x1f)==0x1ff)
      transitdot(1,verblevel);
  }
  gsl_interp_accel_free(acc);


  //set progress indicator, and print output
  tr->pi&=TRPI_MODULATION;
  printmod(tr);  
  return 0;
}



/* \fcnfh
   Printout for modulation as function of wavelength
*/
void
printmod(struct transit *tr)
{
  FILE *out=stdout;
  int rn;

  //open file
  if(tr->f_out&&tr->f_out[0]!='-')
    out=fopen(tr->f_out,"w");

  transitprint(1,verblevel,
	       "\nPrinting in-eclipse/out-eclipse ratio for requested\n"
	       "conditions in '%s'\n"
	       ,tr->f_out?tr->f_out:"standard output");

  //print!
  fprintf(out,
	  "#wavenumber[cm-1]\twavelength[nm]\tmodulation\n");
  for(rn=0;rn<tr->wns.n;rn++)
    fprintf(out,"%12.6f%14.6f%17.7g\n"
	    ,tr->wns.fct*tr->wns.v[rn],WNU_O_WLU/tr->wns.v[rn]/tr->wns.fct,
	    tr->outpret[rn]);

  exit(EXIT_SUCCESS);
}