/*
 * readatm.c - Read atmospheric info main file. Component of the
 *               Transit program.
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

#include "at_common.c"



/* \fcnfh
   frees memory from the onept structure

   @returns 0 on success
*/
void
freemem_onept(struct onept *o)
{

  free(o->q);
  if(o->nm){
    free(o->n[0]);
    free(o->n);
    free(o->m);
  }
}


/* \fcnfh
   frees memory from the atmosphere structure

   @returns 0 on success
*/
int
freemem_atmosphere(struct atm_data *at,
		   long *pi)
{
  //free structures
  free_samp(&at->rads);
  for(int i=0 ; i<at->n_nonignored ; i++)
    free_isov(at->isov+i);
  free_atm(&at->atm);
  free(at->isoprop);


  //free arrays
  free(at->isov);
  free(at->mm);
  free(at->n[0]);
  free(at->n);
  free(at->m);
  free(at->isoeq);
  free(at->isodo);
  free(at->info);


  //clear progress indicator and return success
  *pi&=~(TRPI_GETATM);
  return 0;
}


/* \fcnfh
   saves onept structure's array
*/
void
saveonept_arr(FILE *out,
	      struct onept *onept)
{
  if(onept->ne<=0)
    return;
  fwrite(onept->n[0],1,maxeisoname*onept->ne,out);
  fwrite(onept->m,sizeof(PREC_ZREC),onept->ne,out);
}


/* \fcnfh
   restores onept structure's arrays

   @returns 0 on success
            -1 if not all the expected information is read
	    -2 if info read is wrong
	    -3 if cannot allocate memory
	    1 if information read was suspicious
*/
int
restonept_arr(FILE *in,
	      struct onept *onept)
{
  if(onept->ne<=0)
    return 0;

  if((onept->n=(char **)calloc(onept->ne,sizeof(char *)))==NULL ||
     (onept->n[0]=(char *)calloc(onept->ne,sizeof(char)))==NULL ||
     (onept->m=(PREC_ZREC *)calloc(onept->ne,sizeof(PREC_ZREC)))==NULL)
    return -3;

  int i;
  for(i=1;i<onept->ne;i++)
    onept->n[i]=onept->n[0]+i*maxeisoname;

  if(fread(onept->n[0],1,maxeisoname*onept->ne,in)!=maxeisoname*onept->ne)
    return -1;

  if(fread(onept->m,sizeof(PREC_ZREC),onept->ne,in)!=onept->ne)
    return -1;

  return 0;
}


/* \fcnfh
   Computes mean molecular mass, check that sum of abundances is no
   bigger than 1, and returns it

   @returns sum of abundances
*/
double
checkaddmm(double *mm,		/* mean molecular mass stored here */
	   PREC_NREC r,		/* radius position */
	   prop_isov *isov,	/* variable isotope info */
	   prop_isof *isof,	/* fixed isotope info */
	   int n,		/* number of isotopes in atmfile */
	   _Bool mass,		/* mass abundances? */
	   enum isodo *isodo)	/* whether to ignore isotopes */
{
  double sumq;
  int i;

  if(r >= isov[0].n)
    transiterror(TERR_CRITICAL,
		 "In file %s (line %li) a radius beyond the allocated"
		 " has been requested.", __FILE__, __LINE__);

  //Computes mean molecular mass
  sumq=*mm=0;
  for(i=0;i<n;i++){
    if(isodo[i]!=ignore){
      if(mass)
	*mm+=(isov[i].q[r])/(isof[i].m);
      else
	*mm+=(isov[i].q[r])*(isof[i].m);
      sumq+=isov[i].q[r];
    }
  }
  if(mass)
    *mm=1.0/(*mm);

  //Check that abundances make sense
  if(sumq>1.001){
    transiterror(TERR_SERIOUS,
		 "Sum of abundances of isotopes adds up to more\n"
		 "than 1: %g\n"
		 ,sumq);
  }

  return sumq;
}



/* \fcnfh
   Tell defaults when only one radius is being selected
 */
void
telldefaults(struct isotopes *iso, 
	     struct atm_data *at)
{
  int i;

  transitprint(1,verblevel,
	       "You are using one point atmospheric conditions:\n"
	       " Temperature:         %g K\n"
	       " Pressure:            %g dyne/cm2\n"
	       " Mean molecular mass: %g AMU\n"
	       ,at->atm.t[0]*at->atm.tfct,at->atm.p[0]*at->atm.pfct,at->mm[0]);
  //Densities for all isotopes
  for(i=0;i<iso->n_e;i++)
    if(iso->isodo[i]!=ignore)
      transitprint(1,verblevel,
		   " %-8s: density %8g g/cm3\n"
		   ,iso->isof[i].n,at->isov[i].d[0]);
}


/*\fcnfh
   getatm: Get atmospheric parameters from a given file. It has to get
   all the extra isotopes data that was not in line database.

   @returns 0 On success
            -1 If no atmospheric info file was specified, and no
               defaults are allowed
	    -2 Default handling mode inexistent
	    -3 Something really bad happened!
	    -4 sum of abundances add up to more than 1
	    -5 requested isotope is an ignored one
*/
int getatm(struct transit *tr) /* Containing filename of atmosphere
				       info file */
{
  //'nmb' auxiliar variable for a number.
  int nmb,nrad,i;
  double sumq;
  //'at' and 'hints' are the atmospheric and hint's structure,
  //respectively
  struct transithint *th=tr->ds.th;
  struct onept *onept=&th->onept;
  static struct atm_data st_at;
  prop_samp *rads=&st_at.rads;
  memset(&st_at,0,sizeof(struct atm_data));
  tr->ds.at=&st_at;
  enum {invalid,given,ask,file} inp=invalid;
  float allowq;
  //'fp' and 'newiso' are initialized to avoid compiler output
  FILE *fp=NULL;
  int newiso=0;
  struct isotopes *iso=tr->ds.iso;

  //Pass atmospheric flags into info struct
  transitacceptflag(tr->fl,th->fl,TRU_ATMBITS);
  st_at.mass=th->mass;
  allowq = 1 - (tr->allowrq=th->allowrq);

  //If --onept parameter was given, then force a given one point.
  if(onept->one){
    tr->f_atm=NULL;
    inp=given;
    nrad=rads->n=1;
    rads->fct=1;
  }
  //otherwise if filename is not given or is "-", then use 1 point
  //default behavior
  else if(th->f_atm==NULL||strcmp(th->f_atm,"-")==0){
    tr->f_atm=th->f_atm;
    nrad=rads->n=1;
    rads->fct=1;

    //See how does the user wants the default to be handled
    switch(tr->fl&TRU_ATM1PBITS){
    case TRU_ATMNODEF:		/* Wants an error message */
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "getatm():: No atmospheric file specified\n");
      return -1;
      break;
    case TRU_ATMHARDC1P:	/* wants to use hard-coded values */
      sethcdef(tr,&st_at,rads);
      return 0;
      break;
    case TRU_ATMASK1P:		/* Wants to get values from standard
				   input */
      inp=ask;
      break;
    default:			/* Oops */
      transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
		   "getatm():: Unexistent default handling mode (0x%x)\n"
		   "requested\n"
		   ,tr->fl&TRU_ATM1PBITS);
      return -2;
      break;
    }
  }
  //Otherwise if info is to be taken from file.
  else{
    //If a filename was given check its existence and that it can be
    //opened.
    if((tr->fp_atm=verbfileopen(th->f_atm,"Atmospheric info ")) == NULL)
      exit(EXIT_FAILURE);
    atmfilename=tr->f_atm=th->f_atm;
    inp=file;
    fp=tr->fp_atm;
    transitprint(1,verblevel,
		 "\nReading atmosphere file '%s'...\n"
		 ,atmfilename);

    //nrad willl be the amount of allocated radius so far.
    nrad=8;
    rads->n=0;
    rads->fct=1.0;
  }

  transitassert(inp==invalid,
		"This shouldn't have happened, 'inp' variable is\n"
		"not initialized in file %s, line %i\n"
		,__FILE__,__LINE__);

  //initialize TP arrays, by default temp and pressure are given in cgs.
  rads->v=(PREC_ATM *)calloc(nrad,sizeof(PREC_ATM));
  st_at.atm.tfct=1;
  st_at.atm.pfct=1;
  st_at.atm.t= (PREC_ATM *)calloc(nrad,sizeof(PREC_ATM));
  st_at.atm.p= (PREC_ATM *)calloc(nrad,sizeof(PREC_ATM));
  rads->v[0]=1.0;

  //set p,t and neiso
  switch(inp){
  case given:
    newiso=st_at.n_niso=onept->ne;
    st_at.atm.t[0]=onept->t;
    st_at.atm.p[0]=onept->p;
    break;
  case ask:
    askonenpt(onept,&st_at,-1);
    newiso=st_at.n_niso;
    break;
  case file:
    //T and P will be read from the file, not knowing in advance how
    //many will be.\par
    //Same for the number of extra isotopes, first non-comment lines
    //begginning with i will contain fields with mass and name, 'newiso'
    //will contain the number of currently allocated extra isotopes,
    //'nmb' the total number of allocated isotopes.
    newiso=iso->n_i;
    break;
  default:			/* Just to avoid compiler output, it
				   should never happen */
    transiterror(TERR_CRITICAL,
		 "This is not happening in file %s (%li)\n"
		 ,__FILE__,__LINE__);
    break;
  }

  //Set values and allocate accordingly
  nmb = iso->n_e = iso->n_i+newiso;
  iso->isof = (prop_isof *)realloc(iso->isof,nmb*sizeof(prop_isof));
  iso->isof[iso->n_i].n = (char *)calloc(newiso*maxeisoname,sizeof(char));
  for(i=1;i<newiso;i++)
    iso->isof[iso->n_i+i].n=iso->isof[iso->n_i].n+i*maxeisoname;

  //set mass and names for the new guys
  switch(inp){
  case given:
    //if there were given an exact number of masses and names for the
    //extra isotopes.
    if(onept->nm==newiso){
      for(i=0;i<newiso;i++){
	strcpy(iso->isof[i+iso->n_i].n,onept->n[i]);
	iso->isof[i+iso->n_i].m=onept->m[i];
      }
    }
    else if(onept->nm>0){
      free(onept->m);
      free(onept->n[0]);
      free(onept->n);
    }
  case ask:
    //if there are new isotopes and they were not given exactly. Note
    //that no ignoring is possible if using interactive input.
    if(newiso>0&&onept->nm!=newiso)
      askonemn(onept,iso->isof,newiso, iso->n_i);
    st_at.n=onept->n;
    st_at.m=onept->m;
    break;
  case file:
    if((i=getmnfromfile(fp,&st_at,tr,nmb))<1){
      transiterror(TERR_SERIOUS,
		   "getmnfromfile() returned error code %i\n"
		   ,i);
      exit(EXIT_FAILURE);
    }
    break;
  default:
    break;
  }

  //Allocate isotope information for depth dependent information
  st_at.n_nonignored=nmb=iso->n_e;
  st_at.isov=(prop_isov *)calloc(nmb,sizeof(prop_isov));
  iso->isov=(prop_isov *)realloc(iso->isov,nmb*sizeof(prop_isov));
  st_at.mm=(double *)calloc(nrad,sizeof(double));
  for(i=0;i<nmb;i++){
    st_at.isov[i].d=(PREC_ATM *)calloc(nrad,sizeof(PREC_ATM));
    st_at.isov[i].q=(PREC_ATM *)calloc(nrad,sizeof(PREC_ATM));
    st_at.isov[i].n=nrad;
  }
  iso->isodo=(enum isodo *)realloc(iso->isodo,iso->n_e*sizeof(enum isodo));

  //get isotope abundance
  switch(inp){
  case given:
    if(onept->nq==nmb)
      for(i=0;i<nmb;i++){
	st_at.isov[i].q[0]=onept->q[i];
	iso->isodo[i]=fixed;
      }
    else if(onept->nq>0)
      free(onept->q);
  case ask:
    if(onept->nq!=nmb){
      onept->q=(PREC_ATM *)calloc(nmb,sizeof(PREC_ATM));
      for(i=0;i<nmb;i++){
	st_at.isov[i].q[0]=onept->q[i]
	  =askforposd(" %s abundance for isotope %s: "
		      ,st_at.mass?"Mass":"Number",iso->isof[i].n);
	if(onept->q[i]>=1){
	  transitprint(1,verblevel,
		       " Abundance given is greater than 1, hence ignoring this\n"
		       "isotope from now on\n");
	  iso->isodo[i]=ignore;
	}
	else
	  iso->isodo[i]=fixed;
      }
    }

    if((sumq=checkaddmm(st_at.mm,0,st_at.isov,iso->isof,nmb,st_at.mass,iso->isodo))
       <allowq)
      transiterror(TERR_WARNING,
		   "Abundances don't add up to 1: %.9g under one\n"
		   "point conditions\n"
		   ,sumq);
    //Calculate densities
    for(i=0;i<nmb;i++)
      st_at.isov[i].d[0]=stateeqnford(st_at.mass,
				      st_at.isov[i].q[0],
				      st_at.mm[0],
				      iso->isof[i].m,
				      st_at.atm.p[0]*st_at.atm.pfct,
				      st_at.atm.t[0]*st_at.atm.tfct);
    break;
  case file:
    //All the extra isotopes are going to have density from file (cannot
    //be ignored or fixed). FALSE, they can be FACTOR!
    //    for(i=iso->n_i;i<iso->n_e;i++)
    //      iso->isodo[i]=atmfile;
    nrad=readatmfile(fp,tr,&st_at,rads,nrad);
    transitprint(1,verblevel,
		 " done\n");
    fclose(fp);
    break;
  default:
    break;
  }

  //set required values in 'rads' structure
  rads->i=rads->v[0];
  rads->f=rads->v[rads->n-1];
  rads->o=1;
  rads->d=0;

  //If one point is being used, then let user know values
  if(nrad==1)
    telldefaults(iso,&st_at);

  transitASSERT(iso->n_e<iso->n_i,
		"Uyuyuyy!, number of isotopes after extension (%i)\n"
		"is smaller than number of isotopes before that (%i)\n"
		,iso->n_e,iso->n_i);


  //'tauiso' is the isotope index which tau() will calculate
  //afterwards. Just check that it is between boundaries.
  tr->tauiso=0;
  if(th->tauiso>=0&&
     th->tauiso<iso->n_i){
    if(iso->isodo[th->tauiso]==ignore&&tr->ds.ex->periso){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "Selected isotope to compute tau (#%i: %s) was actually\n"
		   "ignored according to the atmospheric info\n"
		   ,th->tauiso,tr->ds.iso->isof[th->tauiso]);
      return -5;
    }
    else
      tr->tauiso=th->tauiso;
  }

  //Return succes and set progress indicator
  tr->pi|=TRPI_GETATM;
  return 0;
}
