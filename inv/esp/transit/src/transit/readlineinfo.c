/*
 * readlineinfo.c
 * readlineinfo.txc - reads line info as returned by
 *                    lineread. Component of Transit program
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
#include <math.h>

static void asciierr(int max, char *file, int line);
static void notyet(int lin, char *file);

static inline void datafileBS(FILE *fp,	/* File pointer */
			      PREC_NREC initial, /* initial index */
			      PREC_NREC final, /* last index */
			      double lookfor, /* target value */
			      PREC_NREC *resultp, /* result index */
			      int reclength); /* Total length of record */

#define checkprepost(pointer,pre,omit,post) do{                           \
   if(pre)                                                                \
     transiterror(TERR_SERIOUS,                                           \
                  "Pre-condition failed on line %i(%s) while reading:\n"  \
		  "%s\nTWII_Ascii format most likely invalid\n"           \
                  ,__LINE__,__FILE__,line);                               \
   while(omit)                                                            \
     pointer++;                                                           \
   if(post)                                                               \
     transiterror(TERR_SERIOUS,                                           \
                  "Pre-condition failed on line %i(%s) while reading:\n"  \
		  "%s\nTWII_Ascii format most likely invalid\n"           \
                  ,__LINE__,__FILE__,line);                               \
                                             }while(0)


/* \fcnfh
   Read binary TWII file already open in 'fp'

   @returns 0 on success
 */
int 
readtwii_bin(FILE *fp, 
	     struct transit *tr,
	     struct lineinfo *li)
{
  //'rn', 'i' and 'db' are auxiliary.
  //'ndb', 'nT' and 'nIso' number of databases, temp per db and iso per
  //database respectively.
  //'isonames' array with isotope names.
  //'iniw' and 'finw' initial and final wavelength of database.
  //'CS' is an auxiliary cross section pointer.
  //'T' and 'Z' auxiliary temperature and partition function pointers.
  //'acumiso' keeps the cumulative number of isotopes per database.
  double iniw,finw;
  int ndb,db;
  int i,rn;
  int nT,nIso;
  PREC_ZREC *T,*Z;
  PREC_CS *CS;
  int acumiso=0;

  //Read datafile name, initial, final wavelength, and
  //number of datrabases.
  fread(&li->twii_ver,sizeof(int),1,fp);
  fread(&li->twii_rev,sizeof(int),1,fp);
  fread(&iniw,sizeof(double),1,fp);
  fread(&finw,sizeof(double),1,fp);
  fread(&ndb,sizeof(int),1,fp);

  //Allocate pointers according to the number of databases
  tr->db=(prop_db *)calloc(ndb,sizeof(prop_db));
  li->db=(prop_dbnoext *)calloc(ndb,sizeof(prop_dbnoext));

  //Read info for each database
  for(i=0;i<ndb;i++){
    //Allocate and get DB's name
    fread(&rn,sizeof(int),1,fp);
    tr->db[i].n=(char *)calloc(rn+1,sizeof(char));
    fread(tr->db[i].n,sizeof(char),rn,fp);
    tr->db[i].n[rn]='\0';
    
    //Get number of temperature and isotopes
    fread(&nT,sizeof(int),1,fp);
    fread(&nIso,sizeof(int),1,fp);
    li->db[i].t=nT;
    tr->db[i].i=nIso;

    //allocate for variable and fixed isotope info as well as for
    //temperature points.
    T=li->db[i].T=(PREC_ZREC *) calloc(nT,sizeof(PREC_ZREC)  );
      
    //read temperature points
    fread(T,sizeof(PREC_ZREC),nT,fp);
    
    //Update acumulated isotope count
    tr->db[i].s=acumiso;
    acumiso+=nIso;
  }
  //read total number of isotopes.
  fread(&tr->n_i,sizeof(int),1,fp);
  transitASSERT(tr->n_i!=acumiso,
		"Given number of isotopes (%i), doesn't equal\n"
		"real total number of isotopes (%i)\n"
		,tr->n_i,acumiso);

  //allocate structure that are going to receive the isotope info.
  li->isov=(prop_isov *)calloc(tr->n_i,sizeof(prop_isov));
  tr->isof=(prop_isof *)calloc(tr->n_i,sizeof(prop_isof));
  tr->isov=(prop_isov *)calloc(tr->n_i,sizeof(prop_isov));
  transitDEBUG(21,verblevel,
	       "Isotopes:%i\n"
	       "databases: %i\n"
	       "position %li\n"
	       ,tr->n_i,ndb,ftell(fp)); 

  //info for each isotope in database
  for(i=0;i<tr->n_i;i++){
    transitDEBUG(21,verblevel,"isotope %i/%i\n",i,tr->n_i);

    //read database index
    fread(&tr->isof[i].d,sizeof(int),1,fp);

    //set auxiliary variables
    db=tr->isof[i].d;
    nT=li->db[db].t;

    transitDEBUG(21,verblevel,
		 "belongs to DB %i\n"
		 "which have %i temperatures %i isotopes\n"
		 "and starts at isotope %i\n"
		 ,tr->isof[i].d,li->db[db].t, tr->db[db].i,tr->db[db].s);
    //if this is first isotope in database then allocate room for
    //partition function and cross section
    if(tr->db[db].s==i){
      nIso=tr->db[db].i;
      Z=li->isov[i].z=(PREC_ZREC *)calloc(nIso*nT,sizeof(PREC_ZREC));
      CS=li->isov[i].c=(PREC_CS *)calloc(nIso*nT,sizeof(PREC_CS));
      transitDEBUG(21,verblevel,
		   "allocating %i * %i = %i spaces of Z and CS\n"
		   " at %p and %p\n"
		   ,nIso,nT,nIso*nT,(void *)Z,(void *)CS);
    }
    //Otherwise, just position the pointer to the appropiate place in
    //the array. In this part, 'nIso' will indicate the first isotope of
    //the database.
    else{
      //note that for a little while nIso is the first isotope of the
      //database 
      nIso=tr->db[db].s;
      transitASSERT(nIso>=i,
		    "readinfo_twii():: Somehow the first isotope of\n"
		    "current database %i, has an index (%i) greater than\n"
		    "current isotope (%i)\n"
		    ,db,nIso,i);
      Z=li->isov[i].z=li->isov[nIso].z+nT*(i-nIso);
      CS=li->isov[i].c=li->isov[nIso].c+nT*(i-nIso);
      transitDEBUG(21,verblevel,
		   "Partition pointer allocated at %p\n"
		   "and CS at %p\n"
		   ,(void *)Z,(void *)CS);
    }

    //read mass
    fread(&tr->isof[i].m,sizeof(PREC_ZREC),1,fp);

    transitDEBUG(21,verblevel,
		 "Mass read: %g * %g = %g\n"
		 "position: %li, size %i\n"
		 ,tr->isof[i].m,AMU,tr->isof[i].m*AMU
		 ,ftell(fp),sizeof(tr->isof[i].m));

    //allocate and read isotope names
    fread(&rn,sizeof(int),1,fp);
    transitDEBUG(21,verblevel,
		 "Name's length: %i\n"
		 "position: %li, size %i\n"
		 ,rn,ftell(fp),sizeof(int));
    tr->isof[i].n=(char *)calloc(rn+1,sizeof(char));
    fread(tr->isof[i].n,sizeof(char),rn,fp);
    tr->isof[i].n[rn]='\0';

    transitDEBUG(21,verblevel,
		 "Name: %s\n"
		 ,tr->isof[i].n);

    //read partition function
    fread(Z,sizeof(PREC_ZREC),nT,fp);

    /* TD: add cross section in info file */
    //read cross section
    for(rn=0;rn<nT;rn++)
      CS[rn]=SIGWATER;
  }

  //update structure values
  tr->n_db=ndb;
  li->endinfo=ftell(fp);

  li->wi=iniw;
  li->wf=finw;

  return 0;
}


/* \fcnfh
   Read an TWII-ASCII formated file from an already open file 'fp'
*/
int 
readtwii_ascii(FILE *fp, 
	       struct transit *tr,
	       struct lineinfo *li)
{
  char rc;
  char line[maxline+1],*lp,*lp2;
  int ndb,db;
  int nIso,nT,acumiso;
  int rn,i;
  prop_isov *isov;
  PREC_ZREC *T;

  //Format of the TWII-ASCII file is the following(names should not
  //contain spaces ('_' are replaced by spaces)
  //\begin{verb}
  //<m-database>
  //<DATABASE1-name> <n1-iso> <nt1-temp>
  //<NAME1> <MASS1> ... <MASSn1> <MASSn1>
  //<TEMP1>    <Z1-1>  ...  <Z1-n1>   <CS1-1>  ...  <CS1-n1>
  //...
  //<TEMPnt1> <Znt1-1> ... <Znt1-n1> <CSnt1-1> ... <CSnt1-n1>
  //<DATABASE2-name> <n2-iso> <nt2-temp>
  //<NAME1> <MASS1> ... <MASSn2> <MASSn2>
  //<MASS1> ... <MASSn2>
  //<TEMP1>    <Z1-1>  ...  <Z1-n2>   <CS1-1>  ...  <CS1-n2>
  //...
  //<TEMPnt2> <Znt2-1> ... <Znt2-n2> <CSnt2-1> ... <CSnt2-n2>
  //....
  //....
  //<DATABASEm-name> <nm-iso> <ntm-temp>
  //<NAME1> <MASS1> ... <MASSnm> <MASSnm>
  //<TEMP1>    <Z1-1>  ...  <Z1-nm>   <CS1-1>  ...  <CS1-nm>
  //...
  //<TEMPntm> <Zntm-1> ... <Zntm-nm> <CSntm-1> ... <CSntm-nm>
  //<CTRWAV1> <ISOID1> <LOWENER1> <LOGGF1>
  //<CTRWAV2> <ISOID2> <LOWENER2> <LOGGF2>
  //...
  //...
  //\end{verb}
  //get Number of databases from first line
  while((rc=fgetupto(line,maxline,fp,&asciierr,tr->f_line,li->asciiline++))
	=='#'||rc=='\n');
  if(!rc) notyet(li->asciiline,tr->f_line);
  ndb=strtol(line,&lp,0);
  checkprepost(lp,errno&ERANGE,*lp==' ',*lp!='\0');
  //Allocate pointers according to the number of databases
  tr->db=(prop_db *)calloc(ndb,sizeof(prop_db));
  li->db=(prop_dbnoext *)calloc(ndb,sizeof(prop_dbnoext));

  //for each database
  for(db=0;db<ndb;db++){
    //get name
    while((rc=fgetupto(lp=line,maxline,fp,&asciierr,tr->f_line,li->asciiline++))
	  =='#'||rc=='\n');
    if(!rc) notyet(li->asciiline,tr->f_line);
    if((tr->db[db].n=readstr_sp(lp,&lp,'_'))==NULL)
      transitallocerror(0);

    //go to next field and get number of temperatures and isotopes
    checkprepost(lp,0,*lp==' '||*lp=='\t',*lp=='\0');
    rn=getnl(2,' ',lp,&nT,&nIso);
    checkprepost(lp,rn!=2,0,0);
    li->db[db].t=nT;
    tr->db[db].i=nIso;

    //Update acumulated isotope count
    tr->db[db].s=tr->n_i;
    tr->n_i+=nIso;

    //allocate for variable and fixed isotope info as well as for
    //temperature points.
    T=li->db[db].T=(PREC_ZREC *) calloc(nT,sizeof(PREC_ZREC)  );

    //allocate structure that are going to receive the isotope
    //info. If it is not first database then just reallocate
    if(!db){
      isov=li->isov=(prop_isov *)calloc(tr->n_i,sizeof(prop_isov));
      tr->isof=(prop_isof *)calloc(tr->n_i,sizeof(prop_isof));
      tr->isov=(prop_isov *)calloc(tr->n_i,sizeof(prop_isov));
    }
    else{
      isov=li->isov=(prop_isov *)realloc(li->isov,tr->n_i*sizeof(prop_isov));
      tr->isof=(prop_isof *)realloc(tr->isof,tr->n_i*sizeof(prop_isof));
      tr->isov=(prop_isov *)realloc(tr->isov,tr->n_i*sizeof(prop_isov));
    }

    //Allocate cross section and temperature for this database. set
    //isov at the beggining of this database
    isov+=(acumiso=tr->db[db].s);
    isov->z=(PREC_ZREC *)calloc(nIso*nT,sizeof(PREC_ZREC));
    isov->c=(PREC_CS *)calloc(nIso*nT,sizeof(PREC_CS));

    //get isotope name and mass
    while((rc=fgetupto(lp2=line,maxline,fp,&asciierr,tr->f_line,li->asciiline++))
	  =='#'||rc=='\n');
    if(!rc) notyet(li->asciiline,tr->f_line);
    //for each isotope
    for(i=0;i<nIso;i++){
      isov[i].z=isov->z+nT*i;
      isov[i].c=isov->c+nT*i;
      //get name
      if((tr->isof[acumiso+i].n=readstr_sp(lp2,&lp,'_'))==NULL)
	transitallocerror(0);
      //get mass and convert to cgs
      tr->isof[acumiso+i].m=strtod(lp,&lp2);

      if(i!=nIso-1)
	checkprepost(lp2,lp==lp2,*lp2==' '||*lp2=='\t',*lp2=='\0');
    }
    //Last isotope has to be followed by an end of string
    checkprepost(lp2,0,*lp2==' '||*lp2=='\t',*lp2!='\0');

    //get for each temperature
    for(rn=0;rn<nT;rn++){
      //Get a line with temperature, partfcn and cross-sect info
      while((rc=fgetupto(lp=line,maxline,fp,&asciierr,tr->f_line,li->asciiline++))
	    =='#'||rc=='\n');
      if(!rc) notyet(li->asciiline,tr->f_line);
      while(*lp==' ')
	lp++;
      //read temperature
      T[rn]=strtod(lp,&lp);
      checkprepost(lp,*lp=='\0',*lp==' '||*lp=='\t',*lp=='\0');

      //for each isotope in database, read partition function
      for(i=0;i<nIso;i++){
	isov[i].z[rn]=strtod(lp,&lp);
	checkprepost(lp,*lp=='\0',*lp==' '||*lp=='\t',*lp=='\0');
      }
      //and cross section
      for(i=0;i<nIso-1;i++){
	isov[i].c[rn]=strtod(lp,&lp);
	checkprepost(lp,*lp=='\0',*lp==' '||*lp=='\t',*lp=='\0');
      }
      //the last field needs a different check at the end
      isov[i].c[rn]=strtod(lp,&lp);
      checkprepost(lp,0,*lp==' '||*lp=='\t',*lp!='\0');
    }
  }

  //store number of database and beginning of data
  tr->n_db=ndb;
  li->endinfo=ftell(fp);

  //We don't know beforehand what is the range in ascii storage, it has
  //to be looked for
  getinifinasctwii(&li->wi,&li->wf,fp, tr->f_line);

  return 0;
}


/* \fcnfh
   Find initial and final wavelength from file pointed by 'fp'. who is
   at the end of the non-isotope info.

   @returns 0 on success
*/
int
getinifinasctwii(double *ini,	/* where initial value would be stored */
		 double *fin,	/* where initial value would be stored */
		 FILE *fp,	/* file pointer */
		 char *file)	/* File name of file being read */
{
  char rc;
  char line[maxline+1],*lp,*lplast;

  //get first line
  while((rc=fgetupto(lp=line,maxline,fp,&asciierr,file,0))
	=='#'||rc=='\n');
  //there should be at least one isotope!
  if(!rc){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "readlineinfo:: There was no transition info in file '%s',\n"
		 "only general isotope info.\n"
		 ,file);
    exit(EXIT_FAILURE);
  }

  *ini=strtod(line,&lp);
  if(line==lp){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "readlineinfo:: First central wavelength of transitions in\n"
		 "file '%s' is not valid in line:\n%s\n"
		 ,file,line);
    exit(EXIT_FAILURE);
  }

  fseek(fp,0,SEEK_END);
  if(ftell(fp)<maxline){
    transiterror(TERR_WARNING,
		 "readlineinfo:: weird, TWII-Ascii file has less than %i bytes.\n"
		 "That looks improbable\n"
		 ,maxline);
    fseek(fp,0,SEEK_SET);
  }
  else
    fseek(fp,1-maxline,SEEK_END);
  fread(lp=line,sizeof(char),maxline-1,fp);
  line[maxline-1]='\0';

  lplast=NULL;
  while(*lp){
    lp++;
    if(*lp!='\0'&&lp[-1]=='\n')
      lplast=lp;
  }
  
  if(!lplast){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		"Last line in '%s' is longer than %i bytes\n"
		,file,maxline);
    exit(EXIT_FAILURE);
  }

  *fin=strtod(lplast,&lp);
  if(line==lp){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "readlineinfo:: Last central wavelength of transitions in\n"
		 "file '%s' is not valid in line:\n%s\n"
		 ,file,lplast);
    exit(EXIT_FAILURE);
  }

  return 0;
}

/* TD: data info in a structure */

/* \fcnfh
  Checks to see if the hinted parameters regarding wavelength range are
  acceptable to the program. it also sets transit->ds.li.wavs structure.

  @returns 0        if all sugested values were accepted
	   positive if something was changed
	   0x20     if sug final wavelength was changed
	   0x10     if sug initial wavelength was changed
           negative if something bad happen
	   -1       if sug initial wavelength is larger than
	            sug. final
           -2       if sug initial is larger than largest allowed final
	   -3       if sug final is shorter than smallest allowed initial.
	   -4       Margin value too big.
*/
int checkrange(struct transit *tr, /* General parameters and
					 hints */ 
	       struct lineinfo *li) /* Values returned by
					  readinfo\_twii */ 
{
  //'res' is the return value.
  //'margin' margin value.
  //'th' is the hint structure.
  //'iniw' and 'finw' are pointers to hinted initial and final
  //wavelengths.
  //'dbini' and 'dbfin' are the maximum and minimum of the database.
  int res=0;
  PREC_RES margin;
  struct transithint *th=tr->ds.th;
  prop_samp *msamp=&li->wavs;
  prop_samp *hsamp=&th->wavs;
  PREC_LNDATA dbini=li->wi,dbfin=li->wf;
  double extra;

  //initialize modified hints
  msamp->n=-1;
  msamp->d=-1;
  msamp->v=NULL;

  //First check that the margin value is reasonable. i.e. whether its
  //leaves a non-zero range if applied to the whole line dataset.
  if(th->na&TRH_WM){
    if(2*th->m > (li->wf - li->wi)){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "Margin value (%g) is too big for this dataset whose\n"
		   "range is %g to %g nanometers\n"
		   ,th->m,li->wi,li->wf);
      return -4;
    }
    transitaccepthint(margin=tr->m, th->m, th->na, TRH_WM);
  }
  else
    margin=tr->m=0.0;
 
  //If an TWII ascii file, then we'll be using limits twice the margin
  //from the minimum or maximum values in the dataset. This is because,
  //as opposed to binary archiving, I'm not guaranteed that I have all
  //the lines between 'dbini' and 'dbfin', instead they are the minimum
  //and maximum wavelength of the given transitions.
  if(li->asciiline){
    if(margin==0.0)
      transiterror(TERR_WARNING,
		   "Wavelength margin to be used is zero in a TWII-ASCII\n"
		   " file. Hence, there will be no points to the left or\n"
		   " right of the extreme central wavelengths.\n"
		   );
    extra=2*margin;
  }
  else
    extra=0;
  
  transitDEBUG(21,verblevel,
	       "hinted initial %g, final %g\n"
	       "Databse max %g and min %g\n"
	       ,hsamp->i,hsamp->f,li->wi,li->wf);

  //If final wavelength was not hinted then default it to zero
  if(!(th->na&TRH_WF)){
    hsamp->f=0;
    th->na|=TRH_WF;
    transiterror(TERR_WARNING,
		 "Setting hinted upper wavelength limit before\n"
		 "extraction as %g. It was not user-hinted.\n"
		 ,hsamp->f);
  }

  //Check final wavelength. If it is 0, modify it. Do not return
  //special value, it is assumed that user knows what he is doing.
  //\linelabel{finalwav}
  if(hsamp->f<=0){
      msamp->f=dbfin+extra;
  }
  //otherwise
  else{
    transitDEBUG(20,verblevel,
		 "dbini: %g   margin:%g  sampf:%g\n"
		 ,dbini,margin,hsamp->f);

    //check that is not below the minimum value
    if(dbini+margin>hsamp->f){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "Considering margin, final wavelength (%g) is\n"
		   " smaller than first allowed value in database\n"
		   " (%g = %g + %g)\n"
		   ,hsamp->f,dbini+margin,dbini,margin);
      return -3;
    }
    //warn if it is above maximum value with information
    if(hsamp->f>dbfin-margin)
      transiterror(TERR_WARNING,
		   "Final wavelength is above the maximum informative\n"
		   "value in database\n"
		   );
    //set initial wavelength to be extracted
    msamp->f=hsamp->f;
  }

  if(!(th->na&TRH_WI)){
    hsamp->i=0;
    th->na|=TRH_WI;
    transiterror(TERR_WARNING,
		 "Setting hinted lower wavelength limit before\n"
		 "extraction as %g. It was not user-hinted.\n"
		 ,hsamp->i);
  }
  //Check initial wavelength. See final wavelength treatment above for
  //more detailed comments about the code (\lin{finalwav})
  if(hsamp->i<=0)
    msamp->i=dbini-extra;
  else{
    transitDEBUG(20,verblevel,
		 "dbfin: %g   margin:%g  sampi:%g\n"
		 ,dbfin,margin,hsamp->i);
    if(dbfin-margin<hsamp->i){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "Initial wavelength (%g) is larger than larger\n"
		   "allowed value in database (%g = %g + %g),\n"
		   "note that current margin was considered\n"
		   ,hsamp->i,dbfin-margin,dbfin,margin);
      return -2;
    }
    if(hsamp->i<dbini+margin)
      transiterror(TERR_WARNING,
		   "Final wavelength is above the maximum informative\n"
		   "value in database\n"
		   ); 
    msamp->i=hsamp->i;
  }

  //Check that we still have a range considering margin
  if(2*margin>msamp->f-msamp->i){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Usable final (%g) has to be larger than usable\n"
		 "initial wavelength (%g). Note that those values\n"
		 "could have been modified according to the\n"
		 "database range (%g - %g) and margin (%g)\n"
		 ,msamp->i+margin,msamp->f-margin,dbini,dbfin,margin);
    return -1;
  }

  //set progress indicator and return status
  tr->pi|=TRPI_CHKRNG;
  return res;
}



/* \fcnfh
   This function is called if a line of 'file' was longer than 'max'
   characters
*/
static void 
asciierr(int max,		/* Maxiumum length of an accepted line
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
   It outputs error. Used when EOF is found before expected
*/
static void
notyet(int lin, char *file)
{
  transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
	       "readlineinfo:: EOF unexpectedly found at line %i in\n"
	       "ascii-TWII linedb info file '%s'\n"
	       ,lin,file);
  exit(EXIT_FAILURE);
}



/* \fcnfh
  readinfo\_twii: Read infofile as returned by lineread.

  @todo    checks on allocation errors
  @returns 1 on success
           -1 unavailable file
	   -2 Filename not hinted
	   -3 TWII format not valid (missing magic bytes)
	   -4 Improper TWII-ASCII input
*/
int readinfo_twii(struct transit *tr,
		  struct lineinfo *li)
{
  //'fp' file pointer of info file.
  int rn;
  FILE *fp;
  union {char sig[2];short s[2];} sign={{(char)(('T'<<4)|'W'),
					 (char)(('I'<<4)|'I')}};
  char line[maxline+1];
  //Auxiliary hint structure pointer
  struct transithint *th=tr->ds.th;

  li->asciiline=0;

  //Open info file and transport name into transit from hint (setting
  //corresponding flag) if the file opened succesfully
  if(!(th->na&TRH_FL)){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "File information name needs to be hinted...\n"
		 " it was not\n"
		 );
    return -2;
  }
  if((rn=fileexistopen(th->f_line,&fp))!=1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Line info file '%s' is not available.\n"
		 " fileexistopen() error code %i.\n"
		 ,th->f_line,rn);
    return -1;
  }
  tr->fp_line=fp;
  transitaccepthint(tr->f_line,th->f_line,th->na,TRH_FL);

  //Read first two bytes, they should be either `(T<<4|W)(I<<4|I)' or
  //'\#T'. They are stored as integer, so this check also serves to check
  //whether the machine were the data file and the one this program is
  //being run have the same endian order. If the first two are '\#T',
  //then the first line also start as '\#TWII-ascii'
  fread(sign.s+1,sizeof(short),1,fp);
  //is it a binary TWII
  if(sign.s[0]!=sign.s[1]){
    //does it look like being a Ascii TWII?, if so check it.
    rn=strncasecmp((char *)(sign.s+1),"#T",2);
    if(!rn){
      strcpy(line,"#T");
      fread(line+2,sizeof(char),9,fp);
      rn=strncasecmp(line,"#TWII-ascii",11);
    }
    //If it wasn't any valid TWII, then error and exit
    if(rn){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "The file '%s' has not a valid TWII\n"
		   " format. It may also be because the machine were it\n"
		   " was created have different endian order, which is\n"
		   " incompatible"
		   ,tr->f_line);
      return -3;
    }
    li->asciiline=1;
    //ignore the rest of the first line
    fgetupto(line,maxline,fp,&asciierr,tr->f_line,1);
  }

  //if ascii format of TWII read ascii file
  if(li->asciiline){
    if((rn=readtwii_ascii(fp,tr,li))!=0){
      transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
		   "readtwii_ascii() return error code %i\n"
		   ,rn);
      return -5;
    }
  }
  //If binary storage is used.
  else
    if((rn=readtwii_bin(fp,tr,li))!=0){
      transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
		   "readtwii_bin() return error code %i\n"
		   ,rn);
      return -6;
    }

  //close file, set progres indicator and return success.
  fclose(fp);
  tr->pi|=TRPI_READINFO;
  return 1;
}



/* \fcnfh
   print out an error, it is called by readdatarng if one of the field
   with transition info is invalid

   @return -5 always
*/
int
invalidfield(char *line,	/* Contents of the line */
	     char *file,	/* File name */
	     int nmb,		/* File number */
	     int fld,		/* field with the error */
	     char *fldn)	/* Name of the field */
{
  transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
	       "Line %i of file '%s': Field %i (%s) has\n"
	       " not a valid value:\n%s\n"
	       ,nmb,file,fld,fldn,line);
  return -5;
}



/*
  readdatarng: Read a wavelength range from datafile as returned by
  lineread. This function doesn't check for correct boundaries of
  data. You should run checkrange() before running this

  @returns Number of records read on success
           -1 on unexpected EOF
	   -2 file non-seekable
	   -3 on non-integer number of structure records
	   -4 First field is not valid while looking for starting
              point. 
	   -5 One of the fields contained an invalid flaoating point.
*/
int readdatarng(struct transit *tr, /* General parameters and
				       hints */ 
		struct lineinfo *li) /* Values returned by
					readinfo\_twii */ 
{
  //'fp' datafile file pointer.
  //'alloc' number of allocated line\_transition structures.
  //'res' is an auxiliar pointer to tr->lt, where line transition
  //information is kept.
  //'iniw' and 'finw' are auxiliary to keep chosen extractable range.
  //'wltemp' a temporal variable to store wavelength.
  /*
    'dbini' is the first wavelength in database, required for marks.
    PREC_LNDATA dbini=li->wi;
  */
  FILE *fp;
  int rn;
  long i,j;
  long offs;
  PREC_NREC alloc=8;
  struct line_transition *res;
  PREC_LNDATA iniw=li->wavs.i,finw=li->wavs.f;
  PREC_LNDATA wltmp;
  PREC_NREC nfields;
  char line[maxline+1],rc,*lp,*lp2;

  //Open data file 
  if((rn=fileexistopen(tr->f_line,&fp))!=1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "Data file '%s' is not available.\n"
		   " fileexistopen() error code %i.\n"
		   ,tr->f_line,rn);
    return -1;
  }

  //Allocate initial room for line transition's structure
  if((res=(struct line_transition *)
      calloc(alloc,sizeof(struct line_transition)))==NULL)
    transiterror(TERR_CRITICAL|TERR_ALLOC,
		 "Cannot allocate memory for linetran\n"
		 "structure array of length %i, in function readdatarng()\n"
		 ,alloc);

  //find starting point in datafile.\par
  //if it is TWII-ascii file then do a sequential search
  if(li->asciiline){
    //go to where readinfo\_twii left off and skip all comments, since
    //then, 'li->asciiline' won't increase
    fseek(fp,li->endinfo,SEEK_SET);
    while((rc=fgetupto(lp=line,maxline,fp,&asciierr,tr->f_line,li->asciiline++))
	  =='#'||rc=='\n');
    li->endinfo=ftell(fp)-strlen(line)-1;

    //'offs' is the number of lines since the first transition.
    offs=0;
    while(!feof(fp)){
      wltmp=strtod(lp,&lp2);
      //if the first field is not double
      if(lp==lp2){
	transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		     "First field of line %i in file '%s' is not a valid\n"
		     " floating point value:\n%s\n"
		     ,li->asciiline+offs,tr->f_line);
	return -4;
      }

      //if this line's central wavelength is bigger or equal to the
      //requested initial value, then stop loop: we found starting
      //point.
      if(wltmp>=iniw)
	break;

      //skip following comments and increase 'offs'
      while((rc=fgetupto(lp=line,maxline,fp,&asciierr,tr->f_line,
			 li->asciiline+offs++))
	    =='#'||rc=='\n');
    }
    fseek(fp,li->endinfo,SEEK_SET);
  }
  //if is a binary file then we can look through binary search
  else{
    //Check seekability.
    if(fseek(fp,0,SEEK_END)){
      transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
		   "File '%s' was not seekable when trying to go to the end\n"
		   ,tr->f_line);
      return -2;
    }
    //Find number of fields, checking that there are an integer number
    //of them
    offs=li->endinfo;
    j=ftell(fp);
    rn=sizeof(struct line_transition);
    nfields=((j-offs)/rn);
    if(nfields*rn+offs!=j){
      transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
		   "Data file does not have an integer number of records\n"
		   "Initial byte %i, final %i, record size %i\n"
		   ,offs,j,rn);
      return -3;
    }

    //Do the binary search
    datafileBS(fp, offs, nfields, iniw, &j, rn);
    transitDEBUG(21,verblevel,"Beginning found at position %li ",j);
    //check whether we need to start reading from records further back
    //because of repetition of wavelength
    if (j){
      do{
	fseek(fp,offs+--j*sizeof(struct line_transition),SEEK_SET);
	fread(&wltmp,sizeof(PREC_LNDATA),1,fp);
      }while(wltmp>=iniw);
      j++;
    }
    //seek file to starting point.
    transitDEBUG(21,verblevel,"and then slide to %li\n",j);
    fseek(fp,offs+j*sizeof(struct line_transition),SEEK_SET);
  }

  /* Main loop to read all the data. */
  //read the data in the main loop
  i=0;
  while(1){
    //enlarge structure if allocated spaced got filled up
    if(i==alloc)
      if((res=(struct line_transition *)
	  realloc(res,sizeof(struct line_transition)*(alloc<<=1)))
	 ==NULL)
	transiterror(TERR_CRITICAL|TERR_ALLOC,
		     "Cannot enlarge memory allocation area\n"
		     "for array of linetran structure to length %i in\n"
		     "function readdatarng\n",alloc);

    //if ascii, after skipping comments read the 4 fields: center, isoid,
    //lowE, loggf
    if(li->asciiline){
      while((rn=fgetupto(lp=line,maxline,fp,&asciierr,
			 tr->f_line,li->asciiline+offs++))
	    =='#'||rc=='\n');
      //if it is not end of file, read the records.
      if(rn){
	res[i].wl=strtod(lp,&lp2);
	if(lp==lp2) 
	  return invalidfield(line,tr->f_line,li->asciiline+offs
			      ,1,"central wavelength");
	res[i].isoid=strtol(lp2,&lp,0);
	if(lp==lp2)
	  return invalidfield(line,tr->f_line,li->asciiline+offs
			      ,2,"isotope ID");
	res[i].elow=strtod(lp,&lp2);
	if(lp==lp2)
	  return invalidfield(line,tr->f_line,li->asciiline+offs
			      ,3,"lower energy level");
	res[i].gf=strtod(lp2,&lp);
	if(lp==lp2)
	  return invalidfield(line,tr->f_line,li->asciiline+offs
			      ,4,"log(gf)");
      }
    }
    //If binary read a data structure
    else{
      rn=fread(res+i,sizeof(struct line_transition),1,fp);
      transitDEBUG(22,verblevel,"Wavelength:%.8f iso:%i\n",res[i].wl,
		   res[i].isoid);
    }

    //Warn if EOF was found. This should be OK if you want to read
    //until the very last line.
    if(!rn){
      transiterror(TERR_WARNING,
		   "End-of-file in datafile '%s'.\n"
		   "Last wavelength read (%f) was in record %i.\n"
		   "If you are reading the whole range, you can safely ignore\n"
		   "this warning.\n\n"
		   ,tr->f_line,res[i-1].wl,i);
      break;
    }

    /* TD: Skip isotope at read time */

    //End if we passed maximum required wavelength
    if(res[i].wl>finw)
      break;

    i++;
  }
  transitDEBUG(21,verblevel,
	       "Number of lines just read: %li\n"
	       ,i);

  //Realloc final strucuture array to its required size and set the
  //pointer in the transit structure.
  if((tr->lt=
      (struct line_transition *)realloc(res,sizeof(struct line_transition)*i))
     ==NULL){
    transiterror(TERR_CRITICAL|TERR_ALLOC,
		 "Cannot diminish memory allocation area\n"
		 "for array of linetran structure to length %i in\n"
		 "function readdatarng\n",alloc);
    return -1;
  }  

  //close file, set progress indicator and exit the number of read line.
  fclose(fp);
  tr->pi|=TRPI_READDATA;
  return i;
}

/* \fcnfh
  datafileBS: Perform a binary search in file pointed by 'fp'(FILE *)
  between 'initial'(PREC\_NREC) and 'final'(PREC\_NREC) looking for
  'lookfor'(PREC\_NREC) at the first item of a record, result is stored
  in 'resultp'(PREC\_NREC *). Records are of length 'reclength'(int) and
  the first item of each of them is of type PREC\_BREC.
*/
static inline void 
datafileBS(FILE *fp,		/* File pointer */
	   long offs,	        /* initial position of data in twii 
				   file */
	   PREC_NREC nfields,	/* last position */
	   double lookfor,	/* target value */
	   PREC_NREC *resultp,	/* result index */
	   int reclength)	/* Total length of record */
{
  PREC_LNDATA temp;
  const int trglength=sizeof(PREC_LNDATA);
  PREC_NREC ini=0,fin=nfields-1;

  transitDEBUG(21,verblevel,
	       "BS: Start looking from %li in %li fields for %f\n"
	       ,offs,nfields,lookfor);
  do{
    *(resultp)=(fin+ini)/2;
    fseek(fp,offs+reclength*(*resultp),SEEK_SET);
    fread(&temp,trglength,1,fp);
    transitDEBUG(21,verblevel,"BS: found wl %f at position %li\n"
		 ,temp,(*resultp));
    if(lookfor>temp)
      ini=*(resultp);
    else
      fin=*(resultp);
  }while (fin-ini>1);
  *resultp=ini;
}



/* TD: Accept isotope list */
/* TD: Error return codes */
/* TD: check that strtok doesn't overwrite first string */
/* \fcnfh
   readlineinfo: It only function is to call all the functions regarding
   reading TWII information.

   @returns 0 on success
*/
int readlineinfo(struct transit *transit) /* General parameters and
					     hints  */
{
  long rn;
  struct transithint *th=transit->ds.th;
  static struct lineinfo st_li;
  memset(&st_li,0,sizeof(struct lineinfo));
  transit->ds.li=&st_li;

  //Try to read hinted info file
  transitprint(1,verblevel, "Reading info file '%s'... "
	       ,th->f_line);
  if((rn= readinfo_twii(transit,&st_li))!=1)
    transiterror(TERR_SERIOUS,
		 "readinfo_twii() returned an error code %i!\n"
		 ,rn);
  transitprint(1,verblevel, "done%c\n",'.');

  //Check the remainder (margin and range) of the hinted values
  //related to line database reading.
  if((rn=checkrange(transit,&st_li))<0)
    transiterror(TERR_SERIOUS,
		 "checkrange() returned error code %i!\n"
		 ,rn);
  //output status so far if the verbose level is enough
  if(rn>0&&verblevel>1)
    transiterror(TERR_WARNING,
		 "checkrange() modified the suggested parameters,\n"
		 "it returned code 0x%x\n\n"
		 ,rn);
  transitprint(2,verblevel,
	       "   After checking limits, the wavelength range to be\n"
	       "  used is %g to %g. Including a margin of %g,\n"
	       "  the range to be extracted is %g to %g\n"
	       ,transit->ds.li->wavs.i,transit->ds.li->wavs.f
	       ,transit->m,st_li.wi,st_li.wf);

  //read data file
  transitprint(1,verblevel, "Reading data... ");
  if((transit->n_l=readdatarng(transit,&st_li))<1)
    transiterror(TERR_SERIOUS,
		 "readdatarng() returned an error code %li\n"
		 ,transit->n_l);
  transitprint(1,verblevel, "done%c\n",'.');
  transitDEBUG(21,verblevel,
	       "Record 97: wav: %.10g\n"
	       ,transit->lt[97].wl);

  //Status so far
  transitprint(2,verblevel,
	       "OK, So far:\n"
	       " * I read %li records from the datafile.\n"
	       " * The wavelength range read was %.8g to %.8g.\n"
	       " * Current margin is %.4g.\n"
	       " * Usable range is thus %.8g to %.8g.\n"
	       , transit->n_l
	       , st_li.wavs.i,st_li.wavs.f,transit->m
	       ,st_li.wavs.i+transit->m,st_li.wavs.f-transit->m);



#ifndef NODEBUG_TRANSIT
  rn=97; //Some random number to test
  struct line_transition *lp=transit->lt+rn;
  transitDEBUG(21,verblevel,
	       " * And the record %li has the following info\n"
	       "Wavelength: %.10g\n"
	       "Lower Energy Level: %.10g\n"
	       "Log(gf): %.10g\n"
	       "Isotope: %i\n"
	       ,rn,lp->wl,lp->elow, lp->gf
	       ,(lp->isoid));
#endif

  transitDEBUG(21,verblevel,
	       "database min and max: %.10g(%.10g) and %.10g(%.10g)\n"
	       ,st_li.wi,transit->ds.li->wi
	       ,st_li.wf,transit->ds.li->wf);

  return 0;
}


#ifdef DBGREADLINEINFO

int main(int argc, char **argv)
{
  struct transit tr;
  struct transithint th;
  struct lineinfo *li;
  int i,ti1,nbins,ans;
  struct line_transition *lp,*lines;

  tr.ds.th=&th;
  th.na=0;

  verblevel=20;

  th.m=0.001;
  th.na|=TRH_WM;
  char defile_line[]="./res/lineread.inf";
  th.f_line=(char *)calloc(strlen(defile_line)+1,sizeof(char));
  strcpy(th.f_line,defile_line);

  th.na|=TRH_FL;

  nbins=20;
  Pprintf(2,"Number of bins[%i]?: ",nbins);
  if(Pgeti(0,&ti1,6)>0)
    nbins=ti1;

  if((i=readlineinfo(&tr))!=0)
    transiterror(TERR_CRITICAL,
		 "Error code: %i\n"
		 ,i);
  transitDEBUG(20,verblevel,
	       "range: %.10g to %.10g\n"
	       ,tr.ds.li->wi,tr.ds.li->wf);
  li=tr.ds.li;
  lines=tr.lt;

  ti1=(int)(log10(tr.n_l)+1);

  printf("Done reading the file.\n\n"
	 "dbread_pands() test results:\n");
  printf("Chosen wavelength range was from %.10g to %.2f [nm]\n"
	 " %*li lines read\n"
	 " Choosing %i equal-sized bins, the result is\n"
	 ,li->wi,li->wf,ti1,tr.n_l,nbins);

  long qb[tr.n_i];
  float szb=(li->wf-li->wi)/nbins;
  double endb;
 
  lp=lines;
  if(!nbins)
    Pprintf(1,"  hmmm, you chose 0 bins!\n");
  for(i=0;i<nbins;i++){
    memset(qb,0,sizeof(*qb)*4);
    endb=li->wi+(i+1)*szb;
    //    PPprintf(1,2,"KK %g %f\n",lp->wl,endb);
    while((lp->wl)<endb&&lp-lines<tr.n_l){
      qb[lp->isoid]++;
      lp++;
    }

    Pprintf(1," %*i = %i + %i + %i + %i lines shorter than %.3f\n"
	   ,ti1,qb[0]+qb[1]+qb[2]+qb[3],qb[0],qb[1],qb[2],qb[3],endb);
  }


  Pprintf(1,"\nWanna know the value of a single record?\n"
	  "If so, write record number (range 0 - %i), else "
	    "press ^C: "
	  ,tr.n_l-1);

  while(Pgeti(0,&ans,(int)(log10(tr.n_l))+1)>=0){
    if(ans<tr.n_l&&ans>=0){
      lp=lines+ans;
      Pprintf(1,"Wavelength: %.10g\n"
	      ,lp->wl);
      Pprintf(1,"Lower Energy Level: %.10g\nLog(gf): %.10g\n"
	      , lp->elow, lp->gf);
      printf("Isotope Name: %s\n"
	     ,tr.isof[lp->isoid].n);
    }
    else
      Pprintf(1,"\nInvalid record number, so ...");

    Pprintf(1,"\nWanna know the value of another single record?\n"
	    "If so, write the record number (range 0 - %i), else just "
	    "press ^C: "
	    ,tr.n_l-1);
  }
  
}
#endif
