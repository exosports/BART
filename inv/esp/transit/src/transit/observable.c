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
  static struct outputray st_out;
  tr->ds.out=&st_out;

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
  PREC_RES *out=st_out.o=(PREC_RES *)calloc(wn->n,sizeof(PREC_RES));
  struct geometry *sg=tr->ds.sg;
  struct optdepth *tau=tr->ds.tau;

  //set time to the user hinted default
  setgeom(sg,HUGE_VAL,&tr->pi);

  //integrate for each wavelength
  transitprint(1,verblevel,
	       "Integrating for each wavelength. Expect %li dots below...\n"
	       ,wn->n/512);

#ifdef _USE_GSL
  gsl_interp_accel *acc=gsl_interp_accel_alloc();
#endif
  for(w=0;w<wn->n;w++){
    out[w]=sol->obsperwn(tau->t[w],tau->first[w],tau->toomuch,
			 ip,sg,1,acc);

    if((w&0x1ff)==0x1ff)
      transitdot(1,verblevel);
  }
#ifdef _USE_GSL
  gsl_interp_accel_free(acc);
#endif
  transitprint(1,verblevel,"\n");

  //frees no longer needed memory.
  freemem_idexrefrac(tr->ds.ir,&tr->pi);
  freemem_extinction(tr->ds.ex,&tr->pi);
  freemem_tau(tr->ds.tau,&tr->pi);

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
  FILE *outf=stdout;
  struct outputray *outray=tr->ds.out;
  int rn;

  //open file
  if(tr->f_out&&tr->f_out[0]!='-')
    outf=fopen(tr->f_out,"w");

  transitprint(1,verblevel,
	       "\nPrinting in-eclipse/out-eclipse ratio for requested\n"
	       "conditions in '%s'\n"
	       ,tr->f_out?tr->f_out:"standard output");

  //print!
  fprintf(outf,
	  "#wavenumber[cm-1]\twavelength[nm]\tmodulation\n");
  for(rn=0;rn<tr->wns.n;rn++)
    fprintf(outf,"%12.6f%14.6f%17.7g\n"
	    ,tr->wns.fct*tr->wns.v[rn],WNU_O_WLU/tr->wns.v[rn]/tr->wns.fct,
	    outray->o[rn]);

  fclose(outf);

  exit(EXIT_SUCCESS);
}


/*\fcnfh
  Frees memory allocated in struct outputray

  @returns 0 on success
*/
int
freemem_outputray(struct outputray *out,
		  long *pi)
{
  //free arrays
  free(out->o);

  //clear PI and return
  *pi&=~(TRPI_MODULATION);
  return 0;
}
