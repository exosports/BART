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

  transitprint(1,verblevel,
	       "Calculating optical depth at various radius...\n");

#ifdef _USE_GSL
  gsl_interp_accel *acc=gsl_interp_accel_alloc();
#endif
  //for each wavenumber
  for(wi=0;wi<wnn;wi++){
    t=tau.t[wi];

    //For each resultant impact parameter
    for(ii=inn-1;ii>=0;ii--){
      if((t[ii]=rad->fct*tr->sol->tauperb(bb[ii]*riw,r,n,e,inn
					  ,tr->tauiso,wi,dt,acc))
	 >tau.toomuch){
	tau.first[wi]=ii;
	break;
      }
    }

  }
#ifdef _USE_GSL
  gsl_interp_accel_free(acc);
#endif

  transitprint(1,verblevel,
	       "Optical depth calculated up to %g[cm-1]\n"
	       ,tr->ds.tau->toomuch);

  //Print lowest impact parameter before optical gets too big
  if(tr->f_toomuch)
    printtoomuch(tr->f_toomuch,tr->ds.tau, &tr->wns);


  //Set progress indicator and output tau if requested, otherwise return
  //success.
  tr->pi|=TRPI_TAU;
  if(tr->fl&TRU_OUTTAU)
    printtau(tr);
  return 0;
}


/* \fcnfh
   Print lowest impact parameter before optical depth gets too big
*/
void
printtoomuch(char *file, 	/* Filename to save to, a '-' is
				   standard output */
	     struct optdepth *tau, /* Tau information */
	     prop_samp *wn)	/* Wavenumber sampling */
{
  long w;
  FILE *out=stdout;

  //open file
  if(file&&file[0]!='-')
    out=fopen(file,"w");
  if(!out)
    transiterror(TERR_WARNING,
		 "Cannot open '%s' for writing maximum depth before too much\n"
		 "optical depth.\n"
		 ,out==stdout?"STDOUT":file);

  transitprint(1,verblevel,
	       "\nPrinting in '%s' maximum depth before optical depth got\n"
	       "larger than %g and therefore impact parameter was not\n"
	       "calculated for deeper layers.\n"
	       ,file,tau->toomuch);

  fprintf(out,"#Wavelength  Maximum_calculated_depth\n");
  for(w=0;w<wn->n;w++)
    fprintf(out,"%12.3g%12.3g\n",wn->v[w]*wn->fct,tau->first[w]);
	       
}


/* \fcnfh
   Printout for optical depth at a requested radius
*/
void
printtau(struct transit *tr)
{
  int rn;
  FILE *out=stdout;
  prop_samp *rads=&tr->ips;
  PREC_RES **t=tr->ds.tau->t;
  long *frst=tr->ds.tau->first;
  PREC_RES toomuch=tr->ds.tau->toomuch;

  transitcheckcalled(tr->pi,"printtau",1,
		     "tau",TRPI_TAU);

  //open file
  if(tr->f_out&&tr->f_out[0]!='-')
    out=fopen(tr->f_out,"w");
  if(!out)
    transiterror(TERR_WARNING,
		 "Cannot open '%s' for writing optical depth.\n"
		 ,out==stdout?"STDOUT":tr->f_out);

  long rad=
    askforposl("Radius at which you want to print the optical depth(%li - %li): "
	       ,1,rads->n)-1;
  if(rad>rads->n){
    fprintf(stderr,"Value out of range, try again\n");
    printtau(tr);
  }

  transitprint(1,verblevel,
	       "\nPrinting optical depth for radius %li (at %gcm) in '%s'\n"
	       "Optical depth calculated up to %g[cm-1]\n"
	       ,rad+1,rads->fct*rads->v[rad]
	       ,tr->f_out?tr->f_out:"standard output",toomuch);

  //print!
  fprintf(out,
	  "#wavenumber[cm-1]\twavelength[nm]\toptical depth[cm-1]\n");
  for(rn=0;rn<tr->wns.n;rn++)
    fprintf(out,"%12.6f%14.6f%17.7g\n"
	    ,tr->wns.fct*tr->wns.v[rn],WNU_O_WLU/tr->wns.v[rn]/tr->wns.fct,
	    rad<frst[rn]?toomuch:t[rn][rad]);

  exit(EXIT_SUCCESS);
}


/* \fcnfh
   Frees tau structure

   @returns 0 on success
*/
int
freemem_tau(struct optdepth *tau,
	    long *pi)
{
  //frees arrays
  free(tau->t[0]);
  free(tau->t);
  free(tau->first);

  //clear indicator and return success
  *pi&=!(TRPI_TAU);
  return 0;

}
