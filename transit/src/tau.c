/*
 * tau.c   - Finds the optical depth. Component of the Transit program.
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
#include <extraext.h>

#define CIA_DOFLOAT  2
#define CIA_RADFIRST 1

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
  int rn;
  //'bb' is the impact parameter, while 'n' is the index of refraction,
  //'w' is the wavenumber, 't' is the tau as function of impact
  //parameter, 'r' is the radius
  PREC_RES *bb=ip->v;
  PREC_RES *n=tr->ds.ir->n;
  PREC_RES *r=rad->v;
  PREC_RES *tau_wn;
  PREC_ATM *temp=tr->atm.t,tfct=tr->atm.tfct;

  transitcheckcalled(tr->pi,"tau",2,
		     "idxrefrac",TRPI_IDXREFRAC,
		     "extwn",TRPI_EXTWN
		     );

  transitacceptflag(tr->fl,tr->ds.th->fl,TRU_TAUBITS);

  //number of elements
  long int wnn=wn->n;
  long int inn=ip->n;
  long int rnn=rad->n;
  double wfct=wn->fct;

  //set tau structures' value
  struct transithint *trh=tr->ds.th;
  tr->save.ext=trh->save.ext;
  const double blowex=tr->blowex=trh->blowex;
  const int taulevel=tr->taulevel=trh->taulevel;
  tau.toomuch=50;
  if(tr->ds.th->toomuch>0)
    tau.toomuch=trh->toomuch;
  tau.last=(long *)calloc(wnn,sizeof(long));
  tau.t=(PREC_RES **)calloc(wnn,sizeof(PREC_RES *));
  tau.t[0]=(PREC_RES *)calloc(wnn*ip->n,sizeof(PREC_RES));
  for(ii=1;ii<wnn;ii++)
    tau.t[ii]=tau.t[0]+ii*ip->n;

  //set cloud structure
  static struct extcloud cl;
  cl.maxe=tr->ds.th->cl.maxe;
  cl.rini=tr->ds.th->cl.rini;
  cl.rfin=tr->ds.th->cl.rfin;
  if(tr->ds.th->cl.rfct==0)
    cl.rfct=rad->fct;
  else
    cl.rfct=tr->ds.th->cl.rfct;
  tr->ds.cl=&cl;


  _Bool *comp=ex->computed;
  //Restoring savefile if given
  if(tr->save.ext)
    restfile_extinct(tr->save.ext, e, comp, rnn, wnn);


  //compute extinction at the outermost layer
  if(!comp[rnn-1]){
    transitprint(1,verblevel,
		 "Computing extinction in the outtermost layer\n");
    if((rn=computeextradius(rnn-1,
			    tr->atm.t[rnn-1]*tr->atm.tfct, ex))!=0)
      transiterror(TERR_CRITICAL,
		   "computeexradius()returned error code %i\n"
		   ,rn);
  }


  //to temporarily store a per radius info, and the ratio of the ip and
  //rad units
  double rfct=rad->fct;
  double riw=ip->fct/rfct;

  //Need at least three radius to calculate a spline interpolation.
  if(inn<4)
    transiterror(TERR_SERIOUS,
		 "tau(): At least four impact parameters points are\n"
		 " required! (three for spline and one for the analitical\n"
		 " part)"
		 );

  transitprint(1,verblevel,
	       "Calculating optical depth at various radius...\n");

  if(ex->periso)
    transitprint(2,verblevel,
		 " Note that I'm computing only for isotope '%s', others "
		 "were ignored\n"
		 ,tr->ds.iso->isof[tr->tauiso].n);

  //to store reordered extinction
  PREC_RES er[rnn];
  int lastr=rnn-1;
  int wnextout=(long)(wnn/10.0);
  //Following are extinction from scattering and from clouds
  double e_s[rnn], e_c[rnn];
  PREC_CIA **e_cia=tr->ds.cia->e;
  struct extscat *sc=tr->ds.sc;

  //for each wavenumber
  for(wi=0;wi<wnn;wi++){
    tau_wn=tau.t[wi];

    //print output every 10\% that is ready
    if(wi>wnextout){
      transitprint(2,verblevel,
		   "%i%%\r"
		   ,(int)(100*(float)wi/wnn+0.5));
      wnextout+=(long)(wnn/10.0);
    }

    //Calculate extinction coming from scattering, clouds, and CIA for
    //each level
    computeextscat(e_s, rnn, sc, rad->v, rad->fct, 
		   temp, tfct, wn->v[wi]*wfct);
    computeextcloud(e_c, rnn, &cl, rad,
		    temp, tfct, wn->v[wi]*wfct);

    //Put the extinction values in a new array, the values may be
    //temporarily overwritten by (fnc)(), but they should come back as
    //they went in.
    for(ri=0;ri<rnn;ri++)
      er[ri]=e[ri][wi]*blowex+e_s[ri]+e_c[ri]+e_cia[wi][ri];

    if( wi==300 )
      rn=56;

    //For each resultant impact parameter
    for(ri=0;ri<inn;ri++){

      //Compute extinction at new radius if required
      if(bb[ri]*ip->fct<r[lastr]*rfct){

	if(ri)
	  transitprint(3,verblevel,
		       "Last Tau(bb=%9.4g, wn=%9.4g): %10.4g\n"
		       ,bb[ri-1],wn->v[wi],tau_wn[ri-1]);
	//while the extinction at a radius bigger than the impact
	//parameter is not computed.. go for it
	do{
	  if(!comp[--lastr]){
	    //compute a new extinction at given radius printing error if
	    //something happen
	    transitprint(2,verblevel,"Radius %i: %.9g[cm]... "
			 ,lastr+1, r[lastr]);
	    if((rn=computeextradius(lastr,temp[lastr]*tfct,ex))!=0)
	      transiterror(TERR_CRITICAL,
			   "computeextradius() return error code %i while\n"
			   "computing radius #%i: %g\n"
			   ,rn,r[lastr]*rfct);
	    //update the value of the extinction at the right
	    //place.
	    er[lastr] = e[lastr][wi]*blowex + e_s[lastr]
	      + e_c[lastr] + e_cia[wi][lastr];
	  }
	}while(bb[ri]*ip->fct<r[lastr]*rfct);
      }

      if( (tau_wn[ri] = rfct * 
	   fcn(bb[ri]*riw,r+lastr,n+lastr,er+lastr,rnn-lastr,taulevel))
	  > tau.toomuch){
	tau.last[wi]=ri;
	if (ri<3){
	  transitprint(1,verblevel,
		       "WARNING: At wavenumber %g (cm-1), the critical TAU"
		       " value (%g)\n"
		       " was exceeded with tau=%g at the impact parameter"
		       " level %li (%g km), this \n"
		       " should have happened in a deeper layer (check"
		       " IP sampling or ATM file)\n"
		       , wn->v[wi], tau.toomuch, tau_wn[ri]
		       , ri, bb[ri]*rfct/1e5);
	}
	break;
      }
      transitDEBUG(22,verblevel,
		   "Tau(lambda %li=%9.07g, r=%9.4g) : %g  (toomuch: %g)\n"
		   ,wi,wn->v[wi], r[ri], tau_wn[ri],tau.toomuch);
    }

    if(ri==inn){
      transitprint(1,verblevel,
		   "WARNING: At wavenumber %g (cm-1), the bottom of the atmosphere\n"
		   " was reached before obtaining the critical TAU value of %g.\n"
		   " Maximum TAU reached: %g\n",
		   wn->v[wi], tau.toomuch, tau_wn[ri]);
      tau.last[wi] = ri-1;
    }

  }

  transitprint(1,verblevel,
	       " done\nOptical depth calculated up to %g\n"
	       ,tr->ds.tau->toomuch);

  //print detailed output if appropiate
  if(tr->ds.det->tau.n)
    detailout(&tr->wns,&tr->ips,&tr->ds.det->tau,tau.t,0);
  if(tr->ds.det->ext.n)
    detailout(&tr->wns,&tr->rads,&tr->ds.det->ext,e,CIA_RADFIRST);
  if(tr->ds.det->cia.n)
    detailout(&tr->wns,&tr->rads,&tr->ds.det->cia,(double **)e_cia,
	      CIA_DOFLOAT);

  if(tr->save.ext)
    savefile_extinct(tr->save.ext, e, comp, rnn, wnn);

  //Print lowest impact parameter before optical gets too big
  if(tr->f_toomuch)
    printtoomuch(tr->f_toomuch,tr->ds.tau, &tr->wns, &tr->ips);

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


void
outdebtauex(char *name,
	    PREC_RES **e,
	    prop_samp *ip,
	    PREC_RES **t,
	    long rn,
	    long w)
{
  FILE *fp=fopen(name,"w");

  long j;
  for(j=0;j<rn;j++)
    fprintf(fp,"%-15.10g%-15.10g\t%-15.10g\n",ip->v[j],t[w][rn-j-1],e[j][w]);

  fclose(fp);
}

void
outdebex(char *name,
	 PREC_RES **e,
	 PREC_RES *r,
	 long rn,
	 long wi,
	 long wf)
{
  FILE *fp=fopen(name,"w");

  long i,j;
  for(j=0;j<rn;j++){
    fprintf(fp,"%-15.10g\t",r[j]);
    for(i=wi;i<=wf;i++)
      fprintf(fp,"%-15.10g\t",e[j][i]);
    fprintf(fp,"\n");
  }

  fclose(fp);
}
	 

void
outdebtau(char *name,
       prop_samp *ip,
       PREC_RES **t,
       long wi,
       long wf)
{
  FILE *fp=fopen(name,"w");

  long i,j;
  for(j=0;j<ip->n;j++){
    fprintf(fp,"%-15.10g\t",ip->v[j]);
    for(i=wi;i<=wf;i++)
      fprintf(fp,"%-15.10g\t",t[i][j]);
    fprintf(fp,"\n");
  }

  fclose(fp);

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
		 "Cannot open '%s' for writing maximum depth before too\n"
		 " much optical depth.\n"
		 ,out==stdout?"STDOUT":file);

  transitprint(1,verblevel,
	       "\nPrinting in '%s'\n"
	       " maximum depth before optical depth got larger than %g, and\n"
	       " therefore impact parameter was not calculated for deeper\n"
	       " layers.\n\n"
	       ,file,tau->toomuch);

  fprintf(out,"#Wvn[cm-1]  Maximum_calculated_depth\n");
  for(w=0;w<wn->n;w++)
    fprintf(out,"%-14.10g%16.12g\n",wn->v[w]*wn->fct,
	    rad->v[tau->last[w]]*rad->fct);

  fclose(out);
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
    rad=askforposl("Radius at which you want to print the optical "
		   "depth(%li - %li): "
		   ,1,rads->n)-1;
    if(rad>rads->n){
      fprintf(stderr,"Value out of range, try again\n");
      printtau(tr);
    }
  }
  else
    rad=tr->ot;

  transitprint(1,verblevel,
	       "\nPrinting in '%s'\n"
	       " optical depth for radius %li (at %gcm)\n"
	       ,tr->f_out?tr->f_out:"standard output"
	       ,rad+1,rads->fct*rads->v[rad]);

  transitprint(2,verblevel,
	       "Optical depth calculated up to %g[cm-1]\n"
	       ,toomuch);

  //print!
  fprintf(out,
	  "#wavenumber[cm-1]\twavelength[nm]\toptical depth[cm-1]\n");
  for(rn=0;rn<tr->wns.n;rn++)
    fprintf(out,"%12.6f%14.6f%17.7g\n"
	    ,tr->wns.fct*tr->wns.v[rn]
	    ,1/tr->wavs.fct/tr->wns.v[rn]/tr->wns.fct
	    ,rad>last[rn]?toomuch:t[rn][rad]);

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

/* \fcnfh
   Output detailed optical depth info

   @returns 0 on success
 */
int 
detailout(prop_samp *wn,
	     prop_samp *rad,
	     struct detailfld *det,
	     PREC_RES **arr,
	     short flag)
{
  long i,u,d,m;
  _Bool radfirst= (_Bool)(flag&CIA_RADFIRST);
  _Bool dofloat= (_Bool)(flag&CIA_DOFLOAT);

  long idx[det->n];
  double val;
  float **arrf=(float **)arr;
  FILE *out=fopen(det->file,"w");
  if(!out)
    transiterror(TERR_SERIOUS,
		 "Cannot open '%s' for writing fine detail\n"
		 ,det->file);

  transitprint(1,verblevel,
	       "\nPrinting in '%s'\n"
	       " fine detail of %s at selected wavenumbers\n"
	       ,det->file,det->name);

  fprintf(out,"#Radius-w=>    ");
  for(i=0;i<det->n;i++){
    val=det->ref[i];
    u=wn->n-1;
    if(val==wn->v[u])
      d=u;
    else{
      d=0;
      while(u-d>1){
	m=(u+d)/2;
	if(wn->v[m]>val)
	  u=m;
	else
	  d=m;
      }
    }

    idx[i]=d;
    fprintf(out,"%-15.8g",wn->v[idx[i]]);
  }
  fprintf(out,"\n");

  if(radfirst){
    for(m=0;m<rad->n;m++){
      fprintf(out,"%-15.7g",rad->v[m]);
      for(i=0;i<det->n;i++){
	if(dofloat)
	  val=arrf[m][idx[i]];
	else
	  val=arr[m][idx[i]];
	fprintf(out,"%-15.7g",val);
      }
      fprintf(out,"\n");
    }
  }
  else{
    for(m=0;m<rad->n;m++){
      fprintf(out,"%-15.7g",rad->v[m]);
      for(i=0;i<det->n;i++){
	if(dofloat)
	  val=arrf[idx[i]][m];
	else
	  val=arr[idx[i]][m];
	fprintf(out,"%-15.7g",val);
      }
      fprintf(out,"\n");
    }
  }

  fclose(out);

  return 0;
}

