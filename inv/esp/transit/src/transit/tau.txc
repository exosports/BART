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
   symmetric planet.

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
  struct extinction *ex=tr->ds.ex;
  PREC_RES **e=ex->e[tr->tauiso];
  PREC_RES (*fcn)()=tr->sol->tauperb;


  //index, initial and final values
  long ii,wi,ri;
  long inn,rnn,wnn;
  int rn;
  //'bb' is the impact parameter, while 'n' is the index of refraction,
  //'w' is the wavenumber, 't' is the tau as function of impact
  //parameter, 'r' is the radius
  PREC_RES *bb=ip->v;
  PREC_RES *n=tr->ds.ir->n;
  PREC_RES *r=rad->v;
  PREC_RES *t;
  PREC_ATM *temp=tr->atm.t,tfct=tr->atm.tfct;

  transitcheckcalled(tr->pi,"tau",2,
		     "idxrefrac",TRPI_IDXREFRAC,
		     "extwn",TRPI_EXTWN
		     );

  transitacceptflag(tr->fl,tr->ds.th->fl,TRU_TAUBITS);

  //set tau structures' value
  const int taulevel=tr->taulevel=tr->ds.th->taulevel;
  tau.toomuch=50;
  if(tr->ds.th->toomuch>0)
    tau.toomuch=tr->ds.th->toomuch;
  tau.last=(long *)calloc(wn->n,sizeof(long));
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
  double rfct=rad->fct;
  double riw=ip->fct/rfct;

  //Need at least three radius to calculate a spline interpolation.
  if(inn<4)
    transiterror(TERR_SERIOUS,
		 "tau(): At least four impact parameters points are required!.\n"
		 " (three for spline and one for the analitical part)"
		 );

  transitprint(1,verblevel,
	       "Calculating optical depth at various radius...\n");

  if(ex->periso)
    transitprint(2,verblevel,
		 " Note that I'm computing only for isotope '%s', others"
		 "were ignored\n"
		 ,tr->ds.iso->isof[tr->tauiso].n);

  //to store reordered extinction
  PREC_RES er[rnn];
  _Bool *comp=ex->computed;
  int lastr=rnn-1;
  int wnextout=(int)(wnn/10.0);

  //for each wavenumber
  for(wi=0;wi<wnn;wi++){
    t=tau.t[wi];

    //print output every 10\% that is ready
    if(wi>wnextout){
      transitprint(2,verblevel,
		   "%li%%\r"
		   ,100*wi/wnn);
      wnextout+=(int)(wnn/10.0);
    }

    //Put the extinction values in a new array, the values may be
    //temporarily overwritten by (fnc)(), but they should come back as
    //they went in.
    for(ri=0;ri<rnn;ri++)
      er[ri]=e[ri][wi];

    //For each resultant impact parameter
    for(ri=0;ri<inn;ri++){

      //Compute extinction at new radius if required
      if(bb[ri]*ip->fct<r[lastr]*rfct){

	transitprint(3,verblevel,
		     "Last Tau(bb=%9.3g, wn=%9.3g): %10.4g\n"
		     ,bb[ri-1],t[ri-1]);
	//while the extinction at a radius bigger than the impact
	//parameter is not computed.. go for it
	do{
	  if(!comp[--lastr]){
	    //compute a new extinction at given radius printing error if
	    //something happen
	    if((rn=computeextradius(lastr,r[lastr]*rfct,temp[lastr]*tfct,ex))!=0)
	      transiterror(TERR_CRITICAL,
			   "computeextradius() return error code %i while\n"
			   "computing radius #%i: %g\n"
			   ,rn,r[lastr]*rfct);
	    //otherwise, update the value of the extinction at the right place.
	    else
	      er[lastr]=e[lastr][wi];
	  }
	}while(bb[ri]*ip->fct<r[lastr]*rfct);
      }

      if( (t[ri] = rfct * 
	   fcn(bb[ri]*riw,r+lastr,n+lastr,er+lastr,rnn-lastr,taulevel))
	  > tau.toomuch){
	tau.last[wi]=ri;
	break;
      }
    }

  }

  transitprint(1,verblevel,
	       " DONE\nOptical depth calculated up to %g\n"
	       ,tr->ds.tau->toomuch);

  //Print lowest impact parameter before optical gets too big
  if(tr->f_toomuch)
    printtoomuch(tr->f_toomuch,tr->ds.tau, &tr->wns, &tr->rads);

  //free memory that is no longer needed.
  freemem_lineinfotrans(tr->ds.li,&tr->pi);
  freemem_localextinction();

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
	     prop_samp *wn,	/* Wavenumber sampling */
	     prop_samp *rad)	/* Radius sampling */
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
	       "\nPrinting in '%s'\n"
	       " maximum depth before optical depth got larger than %g, and\n"
	       " therefore impact parameter was not calculated for deeper layers.\n\n"
	       ,file,tau->toomuch);

  fprintf(out,"#Wavelength  Maximum_calculated_depth\n");
  for(w=0;w<wn->n;w++)
    fprintf(out,"%-14.10g%16.12g\n",wn->v[w]*wn->fct,
	    rad->v[tau->last[w]]*rad->fct);

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
  long *last=tr->ds.tau->last;
  PREC_RES toomuch=tr->ds.tau->toomuch;

  transitcheckcalled(tr->pi,"printtau",1,
		     "tau",TRPI_TAU);
  tr->ot=tr->ds.th->ot;

  //open file
  if(tr->f_out&&tr->f_out[0]!='-')
    out=fopen(tr->f_out,"w");
  if(!out)
    transiterror(TERR_WARNING,
		 "Cannot open '%s' for writing optical depth.\n"
		 ,out==stdout?"STDOUT":tr->f_out);

  long rad;
  if(tr->ot<0){
    rad=askforposl("Radius at which you want to print the optical depth(%li - %li): "
		   ,1,rads->n)-1;
    if(rad>rads->n){
      fprintf(stderr,"Value out of range, try again\n");
      printtau(tr);
    }
  }
  else
    rad=tr->ot;

  transitprint(1,verblevel,
	       "\nPrinting optical depth for radius %li (at %gcm) in\n"
	       " '%s'\n"
	       ,rad+1,rads->fct*rads->v[rad]
	       ,tr->f_out?tr->f_out:"standard output");

  transitprint(2,verblevel,
	       "Optical depth calculated up to %g[cm-1]\n"
	       ,toomuch);

  //print!
  fprintf(out,
	  "#wavenumber[cm-1]\twavelength[nm]\toptical depth[cm-1]\n");
  for(rn=0;rn<tr->wns.n;rn++)
    fprintf(out,"%12.6f%14.6f%17.7g\n"
	    ,tr->wns.fct*tr->wns.v[rn],1/tr->wavs.fct/tr->wns.v[rn]/tr->wns.fct,
	    rad>last[rn]?toomuch:t[rn][rad]);

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
  free(tau->last);

  //clear indicator and return success
  *pi&=!(TRPI_TAU);
  return 0;

}
