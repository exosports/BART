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

/* TD: data info in a structure */

/*
  datafileBS: Perform a binary search in file pointed by 'fp'(FILE *)
  between 'initial'(PREC\_NREC) and 'final'(PREC\_NREC) looking for
  'lookfor'(PREC\_NREC) at the first item of a record, result is stored
  in 'resultp'(PREC\_NREC *). Records are of length 'reclength'(int) and
  the first item of each of them is of type PREC\_BREC.
*/
static inline void datafileBS(FILE *fp, 	   /* File pointer */
			      PREC_NREC initial,  /* initial index */
			      PREC_NREC final,	   /* last index */
			      double lookfor,     /* target value */
			      PREC_NREC *resultp, /* result index */
			      int reclength);	   /* Total length of record */

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

  //initialize modified hints
  msamp->n=-1;
  msamp->d=-1;
  msamp->v=NULL;

  //First check that the margin value is reasonable. i.e. whether its
  //leaves a non-zero range if applied to the whole line dataset.
  if(2*th->m > (li->wf - li->wi)){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Margin value (%g) is too big for this dataset whose\n"
		 "range is %g to %g nanometers\n"
		 ,th->m,li->wi,li->wf);
    return -4;
  }
  if(th->na&TRH_WM)
    transitaccepthint(margin=tr->m, th->m, th->na, TRH_WM);
  else
    margin=tr->m=0.0;

  transitDEBUG(21,verblevel,
	       "hinted initial %g, final %g\n"
	       "Databse max %g and min %g\n"
	       ,hsamp->i,hsamp->f,li->wi,li->wf);

  if(!(th->na&TRH_WF)){
    hsamp->f=0;
    th->na|=TRH_WF;
    transiterror(TERR_WARNING,
		 "Setting hinted upper wavelength limit before\n"
		 "extraction as %g. It was not user-hinted.\n"
		 ,hsamp->f);
  }
  //Check final wavelength. If it is 0 then use maximun possible less
  //the margin. Do not return special value, it is assumed that user
  //knows what he is doing.
  if(hsamp->f<=0)
    msamp->f=li->wf;
  //otherwise
  else{
    transitDEBUG(20,verblevel,
		 "dbini: %g   margin:%g  sampf:%g\n"
		 ,dbini,margin,hsamp->f);

    //check that is not below the minimum value
    if(dbini+margin>hsamp->f){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "Final wavelength (%g) is smaller than first\n"
		   "allowed value in database (%g = %g + %g),\n"
		   "considering margin\n"
		   ,hsamp->f,dbini+margin,dbini,margin);
      return -3;
    }
    //check and change if it is above maximum allowed value
    if(hsamp->f>dbfin){
      msamp->f=dbfin;
      res|=0x20;
    }
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
  //more detailed comments about the code
  if(hsamp->i<=0)
    msamp->i=li->wi;
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
    if(hsamp->i<dbini){
      msamp->i=dbini;
      res|=0x10;
    }
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
*/
int readinfo_twii(struct transit *tr,
		  struct lineinfo *li)
{
  //'rn', 'i' and 'db' are auxiliary.
  //'ndb', 'nT' and 'nIso' number of databases, temp per db and iso per
  //database respectively.
  //'isonames' array with isotope names.
  //'dnames' array with databases names.
  //'fp' file pointer of info file.
  //'iniw' and 'finw' initial and final wavelength of database.
  //'CS' is an auxiliary cross section pointer.
  //'T' and 'Z' auxiliary temperature and partition function pointers.
  //'acumiso' keeps the cumulative number of isotopes per database.
  int rn,i,db,ndb,nT,nIso;
  FILE *fp;
  double iniw,finw;
  PREC_CS *CS;
  PREC_ZREC *T, *Z;
  int acumiso=0;
  char rc;
  union {char sig[2];short s[2];} sign={{(char)(('T'<<4)|'W'),
					 (char)(('I'<<4)|'I')}};
  int asciiline=0;
  int maxline=200;
  char line[maxline];

  //Auxiliary structure pointers
  struct transithint *th=tr->ds.th;
  struct iso_noext *in=tr->ds.in;

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
    rn=strncmp((char *)&sign.s,"#T",2);
    if(!rn){
      strcpy(line,"#T");
      fread(line+2,sizeof(char),9,fp);
      rn=strncmp(line,"#TWII-ascii",11);
    }
    //If it wasn't any valid TWII, then error and exit
    if(rn){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "The file '%s' has not a valid TWII format. It may\n"
		   "also be because the machine were it was created have\n"
		   "different endian order, which is incompatible"
		   ,tr->f_line);
      return -3;
    }
    asciiline=1;
    //ignore the rest of the first line
    fgetupto(line,maxline,fp,&asciierr,tr->f_line,1);
  }

  if(asciiline){
    while((rc=fgetupto(line,maxline,fp,&asciierr,tr->f_line,asciiline++))
	  =='#');
    if(!rc) notyet(asciiline,tr->f_line);

    /* TD from here!: read ascii-twii input */


  }
  else{
    //Read datafile name, initial, final wavelength, and
    //number of datrabases.
    fread(&li->twii_ver,sizeof(int),1,fp);
    fread(&li->twii_rev,sizeof(int),1,fp);
    fread(&iniw,sizeof(double),1,fp);
    fread(&finw,sizeof(double),1,fp);
    fread(&ndb,sizeof(int),1,fp);

    //Allocate pointers according to the number of databases
    tr->db=(prop_db *)calloc(ndb,sizeof(prop_db));
    in->db=(prop_dbnoext *)calloc(ndb,sizeof(prop_dbnoext));

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
      in->db[i].t=nT;
      tr->db[i].i=nIso;

      //allocate for variable and fixed isotop info as well as for
      //temperature points.
      T=in->db[i].T=(PREC_ZREC *) calloc(nT,sizeof(PREC_ZREC)  );
      
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
  }

  //allocate structure that are going to receive the isotope info.
  in->isov=(prop_isov *)calloc(tr->n_i,sizeof(prop_isov));
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
    nT=in->db[db].t;

    transitDEBUG(21,verblevel,
		 "belongs to DB %i\n"
		 "which have %i temperatures %i isotopes\n"
		 "and starts at isotope %i\n"
		 ,tr->isof[i].d,in->db[db].t, tr->db[db].i,tr->db[db].s);
    //if this is first isotope in database then allocate room for
    //partition function and cross section
    if(tr->db[db].s==i){
      nIso=tr->db[db].i;
      Z=in->isov[i].z=(PREC_ZREC *)calloc(nIso*nT,sizeof(PREC_ZREC));
      CS=in->isov[i].c=(PREC_CS *)calloc(nIso*nT,sizeof(PREC_CS));
      transitDEBUG(21,verblevel,
		   "allocating %i * %i = %i spaces of Z and CS\n"
		   " at %p and %p\n"
		   ,nIso,nT,nIso*nT,(void *)Z,(void *)CS);
    }
    //Otherwise, just position the pointer to the appropiate place in
    //the array. In this part, 'nIso' will indicate the first isotope of
    //the database.
    else{
      nIso=tr->db[db].s;
      transitASSERT(nIso>=i,
		    "readinfo_twii():: Somehow the first isotope of\n"
		    "current database %i, has an index (%i) greater than\n"
		    "current isotope (%i)\n"
		    ,db,nIso,i);
      Z=in->isov[i].z=in->isov[nIso].z+nT*(i-nIso);
      CS=in->isov[i].c=in->isov[nIso].c+nT*(i-nIso);
      transitDEBUG(21,verblevel,
		   "Partition pointer allocated at %p\n"
		   "and CS at %p\n"
		   ,(void *)Z,(void *)CS);
    }

    //read mass
    fread(&tr->isof[i].m,sizeof(PREC_ZREC),1,fp);
    tr->isof[i].m*=AMU;

    transitDEBUG(21,verblevel,
		 "Mass read: %g * %g = %g\n"
		 "position: %li, size %i\n"
		 ,tr->isof[i].m/AMU,AMU,tr->isof[i].m
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

  li->endinfo=ftell(fp);

  //update structure values
  li->wi=iniw;
  li->wf=finw;
  tr->n_db=ndb;
  fclose(fp);

  //set progres indicator and return success.
  tr->pi|=TRPI_READINFO;
  return 1;
}


/*
  readdatarng: Read a wavelength range from datafile as returned by
  lineread. This function doesn't check for correct boundaries of
  data. You should run checkrange() before running this

  @returns Number of records read on success
           -1 on unexpected EOF
	   -2 file non-seekable
	   -3 on non-integer number of structure records
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

  /* Finding starting point in datafile 
     i=(int)((iniw-dbini)/li->dwmark);
     transitDEBUG(20,verblevel,
     "iniw: %.6g    dbini: %.6g  dmark:%g\n"
     ,iniw,dbini,li->dwmark);
     datafileBS(fp, li->mark[i], li->mark[i+1], iniw,
  */
  //find starting point in datafile through binary search.
  if(fseek(fp,0,SEEK_END)){
    transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
		 "File '%s' was not seekable when trying to go to the end\n"
		 ,tr->f_line);
    return -2;
  }
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

  datafileBS(fp, offs, nfields, iniw, &j, rn);
  transitDEBUG(20,verblevel,"Beginning found at position %li ",j);
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
  transitDEBUG(20,verblevel,"and then slide to %li\n",j);
  fseek(fp,offs+j*sizeof(struct line_transition),SEEK_SET);

  /* Main loop to read all the data. */
  //read the data
  i=0;
  while(1){
    //enlarge structure if allocated spaced got filled up
    if(i==alloc)
      if((res=
	  (struct line_transition *)realloc(res,sizeof(struct line_transition)*(alloc<<=1)))
	 ==NULL)
	transiterror(TERR_CRITICAL|TERR_ALLOC,
		     "Cannot enlarge memory allocation area\n"
		     "for array of linetran structure to length %i in\n"
		     "function readdatarng\n",alloc);

    //read a data structure and warn if EOF was found. This should be OK
    //if you want to read until the very last line.
    if(fread(res+i,sizeof(struct line_transition),1,fp)==0){
      transiterror(TERR_WARNING,
		   "End-of-file in datafile '%s'.\n"
		   "Last wavelength read (%f) was in record %i.\n"
		   "If you are reading the whole range ignore this warning.\n\n"
		   ,tr->f_line,res[i-1].wl,i);
      break;
    }
    transitDEBUG(22,verblevel,"Wavelength:%.8f iso:%i\n",res[i].wl,
		 res[i].isoid);

    /* TD: Skip isotope at read time */

    //End if we passed maximum required wavelength
    if(res[i].wl>finw)
      break;

    i++;
  }
  transitDEBUG(20,verblevel,
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
  transitDEBUG(22,verblevel,
	       "inside readdata. Record 97: wav: %.10g\n"
	       ,tr->lt[97].wl);

  //close file, set progress indicator and exit the number of read line.
  fclose(fp);
  tr->pi|=TRPI_READDATA;
  return i;
}


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

  transitDEBUG(20,verblevel,
	       "BS: Start looking from %li in %li fields for %f\n"
	       ,offs,nfields,lookfor);
  do{
    *(resultp)=(fin+ini)/2;
    fseek(fp,offs+reclength*(*resultp),SEEK_SET);
    fread(&temp,trglength,1,fp);
    transitDEBUG(20,verblevel,"BS: found wl %f at position %li\n"
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
   reading TWII information.x
   
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
  static struct iso_noext st_in;
  memset(&st_in,0,sizeof(struct iso_noext));
  transit->ds.in=&st_in;

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
	       " After checking limits, the wavelength range to be\n"
	       "used is %g to %g.\n"
	       "Including a margin of %g,\n"
	       "the range to be extracted is %g to %g\n"
	       ,transit->ds.li->wavs.i,transit->ds.li->wavs.f
	       ,transit->m,st_li.wi,st_li.wf);

  //read data file
  transitprint(1,verblevel, "Reading data... ");
  if((transit->n_l=readdatarng(transit,&st_li))<1)
    transiterror(TERR_SERIOUS,
		 "readdatarng() returned an error code %li\n"
		 ,transit->n_l);
  transitprint(1,verblevel, "done%c\n",'.');
  transitDEBUG(20,verblevel,
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
  transitDEBUG(20,verblevel,
	       " * And the record %li has the following info\n"
	       "Wavelength: %.10g\n"
	       "Lower Energy Level: %.10g\n"
	       "Log(gf): %.10g\n"
	       "Isotope: %i\n"
	       ,rn,lp->wl,lp->elow, lp->lgf
	       ,(lp->isoid));
#endif

  transitDEBUG(20,verblevel,
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
	      , lp->elow, lp->lgf);
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
