/*
 * at_file.c
 * at_file.txc - Read atmospheric info from file. Component of the
 *               Transit program.
 *
 * Copyright (C) 2004 Patricio Rojo (pato@astro.cornell.edu)
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


#include <readatm.h>

static double zerorad=0;

char *atmfilename;
_Bool *isolineinatm;
static struct atm_isoprop *isoprop;

struct fonly {
  char n[maxeisoname];
  PREC_ATM *q;
} *fonly;
static int nfonly=0;



/* \fcnfh
   Add all the abundance from different isotopes except those that are
   going to split the remainder.

   @returns total abundance that is going to be always between 0 and 1
                  inclusive
*/
static inline double
addq(prop_isov *isov,		/* Variable isotope info (abundance
				   among other) */
     enum isodo *isodo,		/* Just add those that are ignored */
     _Bool *other,		/* Don't add those that are going to
				   take care of the remainder to unity
				*/ 
     int n,			/* Number of isotopes */
     PREC_NREC r)		/* Radius at which to compute everything
				 */
{
  double res=0;

  while(--n)
    if(isodo[n]!=ignore && !other[n])
      res+=isov[n].q[r];
  if(*isodo!=ignore && !*other)
    res+=isov->q[r];
  
  if(res>1||res<0)
    transiterror(TERR_SERIOUS,
		 "Without processing 'other' molecules, abundance\n"
		 "addition(%g) is either bigger than 1 or negative!\n"
		 ,res);

  return res;
}


/* \fcnfh
   Finds the abundance of the molecule called 'iso' among the arrays
   given by 'isov' and 'isof' for the radius at level 'r'

   @returns abundance of reference level
*/
static inline double
findfactq(char *iso, 		/* reference isotope looked for */
	  prop_isof *isof, 	/* fixed isotope info (name among
				   others) */
	  prop_isov *isov, 	/* variable isotope info like abundance
				 */
	  int n,		/* Number of isotopes */
	  PREC_NREC r)		/* Radius index from which get the
				   abundance */
{

  while (--n)
    if(strcasecmp(iso,isof[n].n)==0)
      return isov[n].q[r];
  if(strcasecmp(iso,isof->n)==0)
    return isov->q[r];

  transiterror(TERR_SERIOUS,
	       "Isotope you want to reference(%s) was not found among those that\n"
	       "abundance was given.\n"
	       ,iso);

  return -1;
}


/* \fcnfh
   Check that val is positive
*/
static inline void 
checkposvalue(PREC_RES val,	/* value to check */
	      int field,	/* field where it was read */
	      long line)	/* line from which it was read */
{
  if(val<0)
    transiterror(TERR_SERIOUS,
		 "While reading the %ith field in line %li of atmosphere\n"
		 "file %s, a negative value was found (%g)\n"
		 ,field,line-1,atmfilename,val);
}


/* \fcnfh
   Ask what to do with the isotopes that don't have a match in the just
   read atmosphere file

   return number of extra isotopes
*/
static int
checknonmatch(struct transit *tr, /* info about existent isotopes */
	      struct atm_data *at, /* info about just read isotopes */
	      enum isodo *isodo) /* what the user will want with each
				    isotope */ 
{
  int i,j,rn,ison=at->n_aiso;
  struct isotopes *iso=tr->ds.iso;

  //for each of the non-matched isotopes, ask if they want to be
  //matched, ignored or given a fixed value.
  for(i=0;i<iso->n_i;i++){
    if(isodo[i]==unclear){
      //If want to match then ask with what, if ignored or fixed it will
      //be dealt with later
      isodo[i]=askforposl("\nIsotope %s is not in atmosphere file, what do you\n"
			  "want to do? (1:match to some isotope, 2:ignore"
			  " this\nisotope, 3:give a fixed abundance value) "
			  ,iso->isof[i].n);
      if(isodo[i]>3){
	fprintf(stderr,"Invalid value, Try again!:\n");
	isodo[i--]=0;
	continue;
      }
      if(isodo[i]==atmfile){
	while(1){
	  rn=0;
	  for(j=0;j<ison;j++)
	    if(at->isoeq[j]==-1&&at->isodo[j]!=ignore){
	      rn=1;
	      fprintf(stderr,"  %2i: %s (%gAMU)\n", j+1,at->n[j],at->m[j]);
	    }
	  if(!rn){
	    fprintf(stderr,
		    "Sorry but there are no isotopes to match with, so try\n"
		    "another option\n");
	    isodo[i--]=0;
	    continue;
	  }
	  j=askforposd("Select a isotope number from the above list to\n"
		       "match %s with: ",iso->isof[i].n);
	  if(at->isoeq[j-1]==-1){
	    at->isoeq[j-1]=i;
	    break;
	  }
	  else
	    fprintf(stderr,"Invalid value, Try again\n");
	}
      }
    }
  }

  rn=0;
  for(j=0 ; j<ison ; j++)
    if( at->isoeq[j]==-1 && at->isodo[j]!=ignore )
      rn++;

  return rn;
}



/* \fcnfh
   This function is called if a line of 'file' was longer than
   'max' characters 
*/
static void 
atmerr(int max,			/* Maxiumum length of an accepted line
				 */ 
       char *file,		/* File from which we were reading */
       int line)		/* Line who was being read */
{
  transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
	       "Line %i of file '%s' has more than %i characters,\n"
	       "that is not allowed\n"
	       ,file,max);
  exit(EXIT_FAILURE);
}


/* \fcnfh
   print out an error, it is called by readdatarng if one of the field
   with transition info is invalid

*/
static void
invalidfield(char *line,	/* Contents of the line */
	     int nmb,		/* File number */
	     int fld,		/* field with the error */
	     char *fldn)	/* Name of the field */
{
  transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
	       "Line %i of file '%s': Field %i (%s) has\n"
	       " not a valid value:\n%s\n"
	       ,nmb,atmfilename,fld,fldn,line);
  exit(EXIT_FAILURE);
}


/* \fcnfh
   Check whether isotope 'name' found in atmosphere file do
   correspond to an isotope in the lineinfofile for which
   extinction is going to be calculated
*/
static void
isisoline(char *name,		/* Isotope's name */
	  PREC_ATM mass,	/* Isotope's mass */
	  int *isoeq,		/* Isotope's lineinfo position to be stored */
	  enum isodo atisodo,	/* Action for current isotope */
	  prop_isof *isof,	/* Info from lineinfo file */
	  enum isodo *isodo,	/* Action to be taken storage */
	  PREC_NREC nliso)	/* Number of isotopes in the lineinfo */
{
  PREC_NREC i;

  for(i=0;i<nliso;i++){
    if(strcasecmp(name,isof[i].n)==0){
      *isoeq=i;
      if(isolineinatm[i])
	transiterror(TERR_SERIOUS,
		     "Isotope %s has been defined more than once in the\n"
		     "atmosphere file.\n"
		     ,name);
      isolineinatm[i]=1;
      if(isodo[i]!=ignore)
	isodo[i]=atisodo;
      if(isodo[i]==ignore)
	transitprint(2,verblevel,
		     "Ignoring isotope %s (%g AMU)\n"
		     ,isof[i].n,isof[i].m);
      else if(isof[i].m!=mass&&mass!=0)
	transiterror(TERR_WARNING,
		     "Mass of isotope %s, is not the same\n"
		     "in the atmosphere file %s(%g) than in\n"
		     "the transition info file(%g)\n"
		     ,name,atmfilename,mass,
		     isof[i].m);
      break;
    }
  }
  if(*isoeq==-1&&atisodo==ignore)
    nfonly++;
}


/* \fcnfh
   get number of isotopes from file and set index 

   @returns number of lines read
*/
int
getmnfromfile(FILE *fp,
	      struct atm_data *at,
	      struct transit *tr,
	      int nmb)
{
  char line[maxline],*lp;
  int ison=0,i;
  struct isotopes *iso=tr->ds.iso;
  prop_isof *isof=iso->isof;
  enum isodo *isodo=iso->isodo;

  //Set variable to handle proportional to isotopes
  int ipi=0,ipa=at->ipa=4;
  isolineinatm=(_Bool *)calloc(iso->n_i,sizeof(_Bool));
  isoprop=(struct atm_isoprop *)calloc(ipa,sizeof(struct atm_isoprop));

  at->begline=0;
  enum isodo atisodo;
  at->isodo=(enum isodo *)calloc(nmb,sizeof(enum isodo));
  at->isoeq=(int *)calloc(nmb,sizeof(int));
  at->m=(PREC_ZREC *)calloc(nmb,sizeof(PREC_ZREC));
  at->n=(char **)calloc(nmb,sizeof(char *));
  at->n[0]=(char *)calloc(nmb*maxeisoname,sizeof(char));
  at->isoeq[0]=-1;
  for(i=1;i<nmb;i++){
    at->n[i]=at->n[0]+i*maxeisoname;
    at->isoeq[i]=-1;
  }

  //while t,p data doesn't start, check for the various modifiers
  while(1){
    switch(fgetupto(line,maxline,fp,&atmerr,atmfilename,at->begline++)){
    case '\n':			//Ignore comments and
    case '#':			//  blank lines
      continue;

    case 0:			//Error if EOF
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "readatminfo:: EOF unexpectedly found at line %i\n"
		   "of file %s while no t,p data points have been read\n"
		   ,at->begline,atmfilename);
      exit(EXIT_FAILURE);
      continue;

    case 'q':			//Whether is mass or number abundance
      lp=line+1;
      while(*lp++==' ');
      lp--;
      switch(*lp|0x20){
      case 'n':
	at->mass=0;
	break;
      case 'm':
	at->mass=1;
	break;
      default:
	transiterror(TERR_SERIOUS,
		     "'q' option in the atmosphere file can only be followed\n"
		     "by 'm' (for abundances by mass) or 'n' (for abundances by\n"
		     "number). '%s' is invalid.\n"
		     ,line);
	break;
      }
      continue;

    case 'z':			//Zero radius value
      zerorad=atof(line+1);
      continue;

    case 'f':			//An isotope is to be taken as
				//proportional to other.
      lp=line+1;
      while(*lp==' '||*lp=='\t') lp++;

      if(ipi==ipa)
	isoprop=(struct atm_isoprop *)realloc(isoprop,(ipa<<=1)*
					      sizeof(struct atm_isoprop));
      isoprop[ipi].m=getds(lp,0,isoprop[ipi].n,maxeisoname-1);
      //skip over recently read field, and go to next field.
      lp=nextfield(lp);
      //skip an optional equal '=' sign
      if(*lp=='=' && lp[1]==' ')
	lp=nextfield(lp);
      //get factor, which has to be between 0 and 1
      isoprop[ipi].f=strtod(lp,NULL);
      if(isoprop[ipi].f<0 )
	transiterror(TERR_CRITICAL,
		     "Abundance ratio has to be positive in atmosphere\n"
		     "file '%s' in line: %s"
		     ,atmfilename,line);
      lp=nextfield(lp);
      //get name of reference and increase index
      i=0;
      while(*lp)
	isoprop[ipi].t[i++]=*lp++;
      isoprop[ipi].t[i]='\0';

      //now check if that isotope is one of the given in the lineinfo
      //file 
      isoprop[ipi].eq=-1;
      isisoline(isoprop[ipi].n,isoprop[ipi].m,&isoprop[ipi].eq,factor,
		isof,isodo,iso->n_i);

      //advance index and go for the next line
      ipi++;
      continue;

    case 'u':			//Change factorization of radius, temp,
				//or press
      switch(line[1]){
      case 'r':
	at->rads.fct=atof(line+2);
	break;
      case 'p':
	at->atm.pfct=atof(line+2);
	break;
      case 't':
	at->atm.tfct=atof(line+2);
	break;
      default:
	transiterror(TERR_SERIOUS,
		     "Invalid unit factor indication in atmosphere file\n");
	exit(EXIT_FAILURE);
      }
      continue;

    case 'n':			//Name or identifier for file data
      storename(at,line+1);
      continue;

    case 'i':			//Isotope information
      lp=line+1;
      while(*lp==' '||*lp=='\t') lp++;

      //for each field
      while(*lp){
	atisodo=atmfile;
	//Allocate if necessary
	if(ison==nmb){
	  nmb<<=1;
	  at->isodo=(enum isodo *)realloc(at->isoeq,nmb*sizeof(enum isodo));
	  at->isoeq=(int *)realloc(at->isoeq,nmb*sizeof(int));
	  at->m=(PREC_ZREC *)realloc(at->m,nmb*sizeof(PREC_ZREC));
	  char *tmp=at->n[0];
	  at->n=(char **)realloc(at->n,nmb*sizeof(char *));
	  at->n[0]=(char *)realloc(tmp,nmb*maxeisoname*sizeof(char));
	  for(i=1;i<nmb;i++)
	    at->n[i]=at->n[0]+i*maxeisoname;
	  for(i=nmb/2;i<nmb;i++)
	    at->isoeq[i]=-1;
	}

	//get mass and name, checking that is correct. First see if this
	//isotope wants to be ignored.
	if(*lp=='!'){
	  lp++;
	  atisodo=ignore;
	}
	at->m[ison]=getds(lp,0,at->n[ison],maxeisoname-1);
	if(at->m[ison]<0||at->n[ison]=='\0'){
	  transiterror(TERR_SERIOUS,
		       "Invalid field in file %s, line %i while reading isotope"
		       " info at:\n%s\n"
		       ,atmfilename,at->begline,lp);
	}

	//now check if that isotope is one of the given in the lineinfo
	//file
	isisoline(at->n[ison],at->m[ison],at->isoeq+ison,atisodo,
		  isof,isodo,iso->n_i);
	at->isodo[ison++]=atisodo;

	//skip over recently read field, and go to next field.
	while(*lp!=' '&&*lp!='\0') lp++;
	while(*lp==' '||*lp=='\t') lp++;
      }
      continue;

    default:			//T,P seems to be starting
      break;
    }
    break;
  }

  //Check if there was at least an isotope identification and allocate new
  //arrays
  if(!ison)
    transiterror(TERR_SERIOUS,
		 "No isotopes were found in atmosphere file, make sure to\n"
		 "specify them in a line starting with the letter 'i'.\n"
		 "First non-comment line read:\n%s\n"
		 ,line);
  at->begpos=ftell(fp)-strlen(line)-1;

  //shorten extra length of arrays
  fonly=(struct fonly *)calloc(nfonly,sizeof(struct fonly));
  at->ipa=ipa=ipi;
  isoprop=(struct atm_isoprop *)realloc(isoprop,ipa*
				    sizeof(struct atm_isoprop));

  //Makes at arrays bigger, so that they can hold the factorized values.
  nmb=at->n_aiso=ison+ipa;
  at->isodo=(enum isodo *)realloc(at->isodo,nmb*sizeof(enum isodo));
  at->isoeq=(int *)realloc(at->isoeq,nmb*sizeof(int));
  at->m=(PREC_ZREC *)realloc(at->m,nmb*sizeof(PREC_ZREC));
  lp=at->n[0];
  at->n=(char **)realloc(at->n,nmb*sizeof(char *));
  at->n[0]=(char *)realloc(lp,nmb*maxeisoname*sizeof(char));
  for(i=1;i<nmb;i++)
    at->n[i]=at->n[0]+i*maxeisoname;

  //initialize values for the factorized elements
  for(i=ison;i<nmb;i++){
    at->isoeq[i]=isoprop[i-ison].eq;
    at->m[i]=isoprop[i-ison].m;
    strncpy(at->n[i],isoprop[i-ison].n,maxeisoname-1);
    at->n[i][maxeisoname-1]='\0';
    at->isodo[i]=factor;
  }

  //Resolve what to do with those isotopes that appear in the transition
  //database, but not in the atmosphere file.
  at->n_niso=checknonmatch(tr,at,isodo);

  //Set full isotope info in the transit structure
  nmb=iso->n_e=iso->n_i+at->n_niso;
  iso->isof=(prop_isof *)realloc(iso->isof,nmb*sizeof(prop_isof));
  iso->isof[iso->n_i].n=(char *)realloc(iso->isof[iso->n_i].n,
				      nmb*maxeisoname*sizeof(char));
  for(i=1;i<at->n_niso;i++)
    iso->isof[iso->n_i+i].n=iso->isof[iso->n_i].n+i*maxeisoname;
  
  nmb=iso->n_i;
  ipi=0;
  for(i=0;i<ison+ipa;i++)
    if(at->isoeq[i]==-1){
      if(at->isodo[i]!=ignore){
	at->isoeq[i]=nmb;
	iso->isof[nmb].m=at->m[i];
	strcpy(iso->isof[nmb++].n,at->n[i]);
      }
      else{
	at->isoeq[i]=ipi;
      	strcpy(fonly[ipi++].n,at->n[i]);
      }
    }
  at->n_nonignored=nmb;

  transitASSERT(nmb+nfonly!=ison+ipa,
		"Oooops number of ignored-nonline elements (%i) plus the\n"
		"number of nonignored (%i), doesn't match the number of elements\n"
		"found in fields 'i'(%i) and 'f'(%i) of the atmosphere\n"
		"file '%s'\n"
		,nfonly,nmb,ison,ipa,atmfilename);

  transitASSERT(nmb!=iso->n_e,
		"Uyuyuyuyu! Problem in file %s, line %i,\n"
		"assertion failed: %i != %i!!\n"
		,__FILE__,__LINE__,nmb,iso->n_e);

  free(isolineinatm);
  at->isoprop=isoprop;

  return at->begline;
}


/* \fcnfh
   Read abundances and pressure for each isotope and radius

   @returns number of radius point
*/
int
readatmfile(FILE *fp,		/* File */
	    struct transit *tr, /* transit info */
	    struct atm_data *at, /* atmosphere info */
	    prop_samp *rads,	/* radius sampling */
	    int nrad)		/* number of allocated radii, note that
				   is not returned updated */
{
  //find abundance related quantities for each radius
  int lines=at->begline;
  PREC_NREC r=0;
  PREC_RES tmp;
  char rc;
  float allowq=1-tr->allowrq;
  double sumq;
  char line[maxline],*lp,*lp2;
  prop_isov *isov=at->isov;
  int *isoeq=at->isoeq;
  struct isotopes *iso=tr->ds.iso;
  enum isodo *isodo=iso->isodo;
  int i,nmb=iso->n_e;

  for(i=0;i<nfonly;i++)
    fonly[i].q=(PREC_ATM *)calloc(nrad,sizeof(PREC_ATM));

  fseek(fp,at->begpos,SEEK_SET);
  while(1){
    //reallocate if necessary
    if(r==nrad){
      nrad<<=1;
      rads->v=(PREC_ATM *)realloc(rads->v,nrad*sizeof(PREC_ATM));
      at->atm.t= (PREC_ATM *)realloc(at->atm.t,nrad*sizeof(PREC_ATM));
      at->atm.p= (PREC_ATM *)realloc(at->atm.p,nrad*sizeof(PREC_ATM));
      at->mm=(double *)realloc(at->mm,nrad*sizeof(double));
      for(i=0;i<nmb;i++){
	isov[i].d=(PREC_ATM *)realloc(isov[i].d,
				      nrad*sizeof(PREC_ATM));
	isov[i].q=(PREC_ATM *)realloc(isov[i].q,
				      nrad*sizeof(PREC_ATM));
      }
      for(i=0;i<nfonly;i++)
	fonly[i].q=(PREC_ATM *)realloc(fonly[i].q,
				       nrad*sizeof(PREC_ATM));
    }

    //Skip comments and read next line
    while((rc=fgetupto(lp=line,maxline,fp,&atmerr,atmfilename,lines++))
	  =='#'||rc=='\n');
    //if it is end of file, stop loop
    if(!rc)
      break;

    tmp=rads->v[r]=strtod(lp,&lp2)+zerorad;
    checkposvalue(tmp,1,lines);
    if(lp==lp2) 
      invalidfield(line, lines, 1, "radius");
    tmp=at->atm.p[r]=strtod(lp2,&lp);
    checkposvalue(tmp,2,lines);
    if(lp==lp2)
      invalidfield(line, lines, 2, "pressure");
    tmp=at->atm.t[r]=strtod(lp,&lp2);
    checkposvalue(tmp,3,lines);
    if(lp==lp2)
      invalidfield(line, lines, 3, "temperature");

    //now read abundances for every isotope
    for(i=0;i<at->n_aiso;i++){
      //stop reading from file if we reach the factorized elements
      if(isodo[i]==factor)
	continue;
      //Read the abundance of the new element. There are two ways:      
      transitASSERT(isoeq[i]<0 || 
		    (isodo[i]==ignore&&isoeq[i]>=nfonly) || 
		    (isodo[i]!=ignore&&isoeq[i]>=iso->n_e),
		    "Assertion failed in file %s, line %i: %i!=[0,%i]. Fonly: %i\n"
		    ,__FILE__, __LINE__, isoeq[i], 
		    isodo[i]==ignore?nfonly:iso->n_e-1, isodo[i]==ignore);
      //If processing one of the factor only elements
      if(isodo[i]==ignore)
	tmp=fonly[isoeq[i]].q[r]=strtod(lp2,&lp);
      //otherwise if this element is going to be considered
      else
	tmp=isov[isoeq[i]].q[r]=strtod(lp2,&lp);
      checkposvalue(tmp, i+4, lines);

      if(lp==lp2)
	invalidfield(line, lines, 4+i, "isotope abundance");
      lp2=lp;
    }
    int factorfrom=i;

    for(i=0;i<iso->n_i;i++)
      if(isodo[i]==fixed){
	if(!r){
	  isov[i].q[0]=askforposd(" %s abundance for isotope %s: "
				  ,at->mass?"Mass":"Number"
				  ,iso->isof[i].n);
	  if(isov[i].q[0]>=1){
	    fprintf(stderr," Abundance for any single isotope has to be"
		    " less than one\n Try Again!\n");
	    i--;
	  }
	}
	else
	  isov[i].q[r]=isov[i].q[0];
      }

    /* TD: this is WRONG, it goes two loops above. */
    //process the factorized elements if there is any
    if(at->n_aiso>factorfrom){
      double ref;
      int fcti;
      _Bool otherfct[nmb];
      memset(otherfct,0,sizeof(otherfct));
      //go for each of such element defining fcti as the index in that
      //realm.
      for(i=factorfrom;i<at->n_aiso;i++){
	fcti=i-factorfrom;
	//don't process yet those that will use whatever abundance is left
	//to complete unity
	if(strcasecmp(isoprop[fcti].t,"other")==0){
	  otherfct[isoeq[i]]=1;
	  continue;
	}
	//find the reference value
	ref=findfactq(isoprop[fcti].n,iso->isof,isov,nmb,r);
	isov[isoeq[i]].q[r]=isoprop[fcti].f*ref;
      }
      //process factorized elements that will take care of the rest of the
      //atmosphere
      ref=1-addq(isov,isodo,otherfct,nmb,r);
      for(i=factorfrom;i<at->n_aiso;i++)
	if(otherfct[isoeq[i]])
	  isov[isoeq[i]].q[r]=isoprop[i-factorfrom].f*ref;
    }


    //calculate mean molecular mass and check whether abundances add up
    //correctly
    if((sumq=checkaddmm(at->mm+r,r,isov,iso->isof,nmb,at->mass,isodo))
       <allowq)
      transiterror(TERR_WARNING,
		   "In radius %g(%i: %g in file), abundances don't add up to 1: %.9g\n"
		   ,at->rads.v[r],r,at->rads.v[r]-zerorad,sumq);


    //Calculate densities
    for(i=0;i<nmb;i++)
      isov[i].d[r]=stateeqnford(at->mass,
				isov[i].q[r],
				at->mm[r],
				iso->isof[i].m,
				at->atm.p[r]*at->atm.pfct,
				at->atm.t[r]*at->atm.tfct);
    r++;
  }
  rads->n=nrad=r;
  rads->v=(PREC_ATM *)realloc(rads->v,nrad*sizeof(PREC_ATM));
  at->atm.t= (PREC_ATM *)realloc(at->atm.t,nrad*sizeof(PREC_ATM));
  at->atm.p= (PREC_ATM *)realloc(at->atm.p,nrad*sizeof(PREC_ATM));
  at->mm=(double *)realloc(at->mm,nrad*sizeof(double));
  for(i=0;i<nmb;i++){
    isov[i].d=(PREC_ATM *)realloc(isov[i].d,
				  nrad*sizeof(PREC_ATM));
    isov[i].q=(PREC_ATM *)realloc(isov[i].q,
				  nrad*sizeof(PREC_ATM));
  }

  //free arrays that were used only to get the factorizing elements
  for(i=0;i<nfonly;i++)
    free(fonly[i].q);
  free(fonly);


  return nrad;
}


/* \fcnfh
   Stores info about the atmopshere file
*/
void
storename(struct atm_data *at,
	  char *line)
{
  while(*line==' '||*line=='\t') line++;

  int len=strlen(line);

  //only store name if it has not been stored before
  if(!at->info){
    at->info=calloc(len+1,sizeof(char));
    strcpy(at->info,line);
  }
}



