/*
 * cia.c   - Computes collision-indeced absorption. Component of the
 *           Transit program.
 *
 * Copyright (C) 2004-2006 Patricio Rojo (pato@astro.cornell.edu)
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


/* \fcnfh
   Read CIA info from given tabulated files

   @returns 0 on success
*/
int
interpolatecia(struct transit *tr)
{
  FILE *fp;
  char *file,*colname1=0,*colname2=0;
  double **a,*t=NULL,*wn;

  static struct cia st_cia;
  tr->ds.cia=&st_cia;
  int npairs=tr->ds.cia->n=tr->ds.th->ncia, p;
  long nwn,nt=0;
  char rc;
  char *lp,*lpa;
  int maxline=300,n;
  long lines,i;
  char line[maxline+1];
  double *tmpw = malloc(tr->wns.n  * sizeof(double));
  double *tmpt = malloc(tr->rads.n * sizeof(double));
  //  double tmpw[tr->wns.n],tmpt[tr->rads.n];
  prop_atm *atm=&tr->atm;
  double *densiso1,*densiso2,amagat2;
  struct isotopes *iso=tr->ds.iso;

  transitcheckcalled(tr->pi,"interpolatecia",2,
		     "makewnsample",TRPI_MAKEWN,
		     "makeradsample",TRPI_MAKERAD
		     );

 //invert temperature becase it increases inwards.
  for(i=0;i<tr->rads.n;i++)
    tmpt[i]=atm->tfct*atm->t[tr->rads.n-i-1];
  for(i=0;i<tr->wns.n;i++)
    tmpw[i]=tr->wns.fct *tr->wns.v[i];


  //Return succes if there is no files to read, or give some info
  if(!npairs){
    st_cia.e=(PREC_CIA **)calloc(tr->wns.n,sizeof(PREC_CIA *));
    st_cia.e[0]=(PREC_CIA *)calloc(tr->wns.n*tr->rads.n,sizeof(PREC_CIA));
    for(p=1;p<tr->wns.n;p++)
      st_cia.e[p]=st_cia.e[0]+p*tr->rads.n;
    return 0;
  }
  transitprint(1,verblevel,
	       "\nComputing CIA opacities for %i database%s...\n"
	       ,npairs,npairs>1?"s":"");

  //prepare and read each datafile
  st_cia.file=(char **)calloc(npairs,sizeof(char *));

  double **e;
  e=(double **)calloc(tr->wns.n,sizeof(double *));
  e[0]=(double *)calloc(tr->wns.n*tr->rads.n,sizeof(double));
  for(p=1;p<tr->wns.n;p++)
    e[p]=e[0]+p*tr->rads.n;
  memset(e[0],0,tr->wns.n*tr->rads.n*sizeof(double));

  st_cia.e=(PREC_CIA **)calloc(tr->wns.n,sizeof(PREC_CIA *));
  st_cia.e[0]=(PREC_CIA *)calloc(tr->wns.n*tr->rads.n,sizeof(PREC_CIA));
  for(p=1;p<tr->wns.n;p++)
    st_cia.e[p]=st_cia.e[0]+p*tr->rads.n;

  for(p=0;p<npairs;p++){
    file = st_cia.file[p] = strdup(tr->ds.th->ciafile[p]);

    if((fp=fopen(file,"r"))==NULL)
      transiterror(TERR_SERIOUS,
		   "Cannot read file '%s' for CIA information\n"
		   ,file);
    lines=0;
    lpa=0;
    //following loop will stop right after reading the 't' line
    //(ignoring blanks and comments)
    while(1){
      //Skip comments,blanks and read next line
      while((rc=fgetupto_err(lp=line,maxline,fp,&ciaerr,file,lines++))
	    =='#'||rc=='\n');
      //if it is end of file, stop loop
      if(!rc)
	transiterror(TERR_SERIOUS,
		     "File '%s' finished before any info about CIA\n"
		     " opacities\n"
		     ,file);

      switch(rc){
      case 'i':
	while(isblank(*++lp));
	long ni=countfields(lp,',');
	if(ni!=2)
	  transiterror(TERR_SERIOUS,
		       "Wrong line %i in CIA file '%s', if it begins with\n"
		       " a 'i' then it should have the comma-separated\n"
		       " fields with the names of the isotopes in\n"
		       " collision. Rest of line:\n"
		       " %s\n"
		       ,lines,file,lp);
	ni=0;
	while(lp[ni++]!=',');
	colname1=(char *)calloc(ni,sizeof(char));
	strncpy(colname1,lp,ni-1);
	lp+=ni-1;
	ni=0;
	while(isblank(*lp++));
	while(lp[ni++]);
	if((--ni)<=0)
	  transiterror(TERR_SERIOUS,
		       "Inexistent second isotope name in CIA file '%s'\n"
		       ,file);
	colname2=(char *)calloc(ni+1,sizeof(char));
	strncpy(colname2,lp,ni);
	continue;

      case 't':
	while(isblank(*++lp));
	nt=countfields(lp,0);
	if(!nt)
	  transiterror(TERR_SERIOUS,
		       "Wrong line %i in CIA file '%s', if it begins with\n"
		       " a 't' then it should have the blank-separated\n"
		       " fields with the temperatures. Rest of line:\n"
		       " %s\n"
		       ,lines,file,lp);
	t=(double *)calloc(nt,sizeof(double));
	n=0;
	lpa=lp;
	while(n<nt){
	  while(isblank(*lpa++));
	  t[n]=strtod(--lpa,&lp);
	  if(lp==lpa)
	    transiterror(TERR_CRITICAL,
			 "Less fields(%i) than expected(%i) were read\n"
			 "for temperature in the CIA file '%s'\n"
			 ,n,nt,file);
	  if((lp[0]|0x20)=='k') lp++;
	  lpa=lp;
	  n++;
	}
	continue;
      default:
	break;
      }
      break;
    }

    if(!colname2||nt==0)
      transiterror(TERR_SERIOUS,
		   "There is either a temperature 't' or isotope name 'i'\n"
		   " missing in CIA file '%s'\n"
		   ,file);

    //set an initial value for allocated wavenumber fields
    long wa=32;
    wn=(double *)calloc(wa,sizeof(double));
    a=(double **)calloc(wa,sizeof(double *));
    a[0]=(double *)calloc(wa*nt,sizeof(double));
    for(i=1;i<wa;i++)
      a[i]=a[0]+i*nt;
    n=0;
    //this loop will read all the information
    while(1){
      //Skip comments,blanks and read next line
      if (n)
	while((rc=fgetupto_err(lp=line,maxline,fp,&ciaerr,file,lines++))
	      =='#'||rc=='\n');
      //if it is end of file, stop loop
      if(!rc)
	break;

      //enlarge allocation if necessary
      if(n==wa){
	wn=(double *)realloc(wn, sizeof(double)*(wa<<=1));
	a=(double **)realloc(a, sizeof(double)*wa);
	a[0]=(double *)realloc(a[0], sizeof(double)*wa*nt);
	for(i=1;i<wa;i++)
	  a[i]=a[0]+i*nt;
      }

      //store new line, starting with the wavenumber and then looping
      //over the cross-sections
      while(isblank(*lp++));
      wn[n]=strtod(lp-1,&lpa);
      if(lp==lpa+1)
	transiterror(TERR_CRITICAL,
		     "Not even one valid field for the %ith wavenumber\n"
		     " in the CIA file '%s'\n"
		     ,n+1,file);
      i=0;
      while(i<nt){
	a[n][i]=strtod(lpa,&lp);
	if(lp==lpa)
	  transiterror(TERR_CRITICAL,
		       "Less fields(%i) than expected(%i) were read\n"
		       "for the %ith wavenumber in the CIA file '%s'\n"
		       ,i,nt,n+1,file);
	lpa=lp;
	i++;
      }
      n++;
    }

    //diminish allocation if necessary
    if(n<wa){
      wn=(double *)realloc(wn, sizeof(double)*n);
      a[0]=(double *)realloc(a[0], sizeof(double)*n*nt);
      for(i=1;i<n;i++)
	a[i]=a[0]+i*nt;
    }
    nwn=n;

    //Now interpolate the data to the correct sampling
    bicubicinterpolate(e,a,wn,nwn,t,nt,
		       tmpw,tr->wns.n, tmpt,tr->rads.n);


    //TD: Linearly adding CIA when there is more than 1 such
    //database. Check that this is really OK.
    densiso1=densiso2=0;
    long j;
    for(j=0;j<iso->n_e;j++){
      if(strcmp(iso->isof[j].n,colname1)==0)
	densiso1=iso->isov[j].d;
      if(strcmp(iso->isof[j].n,colname2)==0)
	densiso2=iso->isov[j].d;
    }

    if(!densiso1||!densiso2)
      transiterror(TERR_SERIOUS,
		   "One or both of the names of the isotopes in CIA (%s, %s)\n"
		   " file '%' doesn't match any in the atmsopheric database '%s'"
		   ,colname1,colname2,file,tr->f_atm);

    nt=tr->rads.n;
    for(i=0;i<nt;i++){
      amagat2=densiso1[i]*densiso2[i]/RHOSTP/RHOSTP;
      for(j=0;j<tr->wns.n;j++)
	st_cia.e[j][i] += e[j][nt-i-1]*amagat2;
    }

    free(a[0]);
    free(a);
    free(wn);
    free(t);
    free(colname1);
    free(colname2);
    fclose(fp);
  }

  free(e[0]);
  free(e);
  free(tmpw);
  free(tmpt);

  transitprint(1,verblevel,
	       " done\n");

  tr->pi|=TRPI_CIA;
  return 0;
}


/* \fcnfh
   Interpolates 'src' into 'res' according to the new dimensions, first
   interpolates the second dimension and then the first. The result is
   added to whatever it is already existent in 'res'

   @returns 0 on success
*/
int
bicubicinterpolate(double **res, /* target array [t1][t2] */
		   double **src, /* Source array [x1][x2] */
		   double *x1,
		   long nx1,
		   double *x2,
		   long nx2,
		   double *t1,
		   long nt1,
		   double *t2,
		   long nt2)
{
  long i,j;
  double fx1=x1[0],fx2=x2[0],lx1=x1[nx1-1],lx2=x2[nx2-1];
  long lj=nt2,fj=0;
  long li=nt1,fi=0;

#ifndef _USE_GSL
#error We cannot spline interpolate without GSL to obtain CIA opacities
#endif

  if(t1[0]>lx1||t1[nt1-1]<fx1||
     t2[0]>lx2||t2[nt2-1]<fx2)
    return 0;


  while(t1[fi++]<fx1);
  fi--;
  for(i=0;i<li;i++)
    if(t1[i]>lx1)
      li=i;
  while(t2[fj++]<fx2);
  fj--;
  for(j=0;j<lj;j++)
    if(t2[j]>lx2)
      lj=j;

  double **f2 = (double **)malloc(nt2*sizeof(double *));
  f2[0] = (double *)malloc(nt2*nx1*sizeof(double));
  for(i=1 ; i<nt2 ; i++)
    f2[i] = f2[0] + i*nx1;
  gsl_interp_accel *acc;
  gsl_interp       *spl;

  for(i=0 ; i<nx1 ; i++){
    acc = gsl_interp_accel_alloc ();
    spl = gsl_interp_alloc(gsl_interp_cspline,nx2);
    gsl_interp_init(spl,x2,src[i],nx2);
    for(j=fj;j<lj;j++)
      f2[j][i]=gsl_interp_eval(spl,x2,src[i],t2[j],acc);
    gsl_interp_free(spl);
    gsl_interp_accel_free(acc);
  }

  for(j=fj ; j<lj ; j++){
    acc = gsl_interp_accel_alloc ();
    spl = gsl_interp_alloc(gsl_interp_cspline,nx1);
    gsl_interp_init(spl,x1,f2[j],nx1);
    for(i=fi;i<li;i++)
      res[i][j]+=gsl_interp_eval(spl,x1,f2[j],t1[i],acc);
    gsl_interp_free(spl);
    gsl_interp_accel_free(acc);
  }


  free(f2[0]);
  free(f2);

  return 0;
    
}

  
/* \fcnfh
   This function is used when a longer than expected line is found in the
   CIA file
*/
void 
ciaerr(int max,
       char *name,
       int line)
{
  transiterror(TERR_SERIOUS,
	       "Line %i of CIA file '%s' is longer than %i characters...\n"
	       " hard coded values in file '%s' need to be changed\n"
	       ,line,name,max,__FILE__);
}


/* \fcnfh
   frees cia structure 

   @returns 0 on success
*/
int
freemem_cia(struct cia *cia,
	    long *pi)
{
  //free arrays
  free(cia->e[0]);
  free(cia->e);
  for (int i=0 ; i<cia->n ; i++)
    free(cia->file[i]);
  if(cia->n)
    free(cia->file);

  //unset appropiate flags
  *pi&=~(TRPI_CIA);
  return 0;
}
