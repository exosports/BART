/*
 * readlineinfo.c   - reads line info as returned by
 *                    lineread. Component of Transit program
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

#define TLI_WAV_UNITS 1e-4  //Microns as of TLI v4
#define TLI_E_UNITS   1

static double tli_to_microns = TLI_WAV_UNITS/1e-4;

static void notyet(int lin, char *file);
static int invalidfield(char *line, char *file, int nmb,
			int fld, char *fldn);


#define checkprepost(pointer,pre,omit,post) do{                            \
   if(pre)                                                                 \
     transiterror(TERR_SERIOUS,                                            \
                  "Pre-condition failed on line %i(%s)\n while reading:\n" \
		  "%s\n\nTLI_Ascii format most likely invalid\n"           \
                  ,__LINE__,__FILE__,line);                                \
   while(omit)                                                             \
     pointer++;                                                            \
   if(post)                                                                \
     transiterror(TERR_SERIOUS,                                            \
                  "Post-condition failed on line %i(%s)\n while reading:\n"\
		  "%s\n\nTLI_Ascii format most likely invalid\n"           \
                  ,__LINE__,__FILE__,line);                                \
                                             }while(0)


/* \fcnfh
  datafileBS: Perform a binary search in file pointed by 'fp'(FILE *)
  between 'initial'(PREC\_NREC) and 'final'(PREC\_NREC) looking for
  'lookfor'(PREC\_NREC) at the first item of a record, result is stored
  in 'resultp'(PREC\_NREC *). Records are of length 'reclength'(int) and
  the first item of each of them is of type PREC\_BREC.
*/
static inline void 
datafileBS(FILE *fp,		/* File pointer */
	   long offs,	        /* initial position of data in tli 
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
    transitDEBUG(21,verblevel,"BS: found wl %f microns at position %li\n"
		 ,temp*tli_to_microns,(*resultp));
    if(lookfor>temp)
      ini=*(resultp);
    else
      fin=*(resultp);
  }while (fin-ini>1);
  *resultp=ini;
}

/* \fcnfh
   Read binary TLI file already open in 'fp'

   @returns 0 on success
 */
int 
readtli_bin(FILE *fp, 
	     struct transit *tr,
	     struct lineinfo *li)
{
  //'rs' is auxiliary.
  //'ndb', 'nT' and 'nIso' number of databases, temp per db and iso per
  //database respectively.
  //'isonames' array with isotope names.
  //'iniw' and 'finw' initial and final wavelength of database.
  //'CS' is an auxiliary cross section pointer.
  //'T' and 'Z' auxiliary temperature and partition function pointers.
  //'acumiso' keeps the cumulative number of isotopes per database.
  //'correliso' indicates this isotope's correlative number.
  double iniw,finw;
  unsigned short ndb;
  unsigned short rs;
  unsigned short nT,nIso;
  PREC_ZREC *T,*Z;
  PREC_CS *CS;
  int acumiso=0;
  struct isotopes *iso=tr->ds.iso;
  int correliso = 0;

  //Read datafile name, initial, final wavelength, and
  //number of databases.
  fread(&li->tli_ver,sizeof(unsigned short),1,fp);
  fread(&li->lr_ver ,sizeof(unsigned short),1,fp);
  fread(&li->lr_rev ,sizeof(unsigned short),1,fp);

  if(li->tli_ver!=compattliversion)
    transiterror(TERR_SERIOUS,
		 "The version of the TLI file: %i (lineread v%i.%i) is not compatible with\n"
		 "this version of transit, which can only read version %i\n"
		 ,li->tli_ver,li->lr_ver,li->lr_rev,compattliversion);

  fread(&iniw,sizeof(double),1,fp);
  fread(&finw,sizeof(double),1,fp);
  fread(&rs,sizeof(unsigned short),1,fp);
  char undefined_string[rs+1];
  fread(undefined_string,sizeof(char),rs,fp);
  undefined_string[rs]='\0';
  fread(&ndb,sizeof(unsigned short),1,fp);

  //Allocate pointers according to the number of databases
  iso->db=(prop_db *)calloc(ndb,sizeof(prop_db));
  li->db=(prop_dbnoext *)calloc(ndb,sizeof(prop_dbnoext));

  //Allocate isotope info of size 1, it will be reallocated on demand as
  //new isotopes are read.
  li->isov    = (prop_isov *)calloc(1, sizeof(prop_isov));
  iso->isof   = (prop_isof *)calloc(1, sizeof(prop_isof));
  iso->isov   = (prop_isov *)calloc(1, sizeof(prop_isov));
  iso->isodo  = (enum isodo *)calloc(1, sizeof(enum isodo));
  li->isov->z = (double *)calloc(1,sizeof(double));
  li->isov->c = (double *)calloc(1,sizeof(double));


  //Read info for each database
  for(short i=0 ; i<ndb ; i++){
    //Allocate and get DB's name
    fread(&rs, sizeof(unsigned short), 1, fp);
    iso->db[i].n = (char *)calloc(rs+1, sizeof(char));
    fread(iso->db[i].n, sizeof(char), rs, fp);
    iso->db[i].n[rs] = '\0';
    
    //Get number of temperature and isotopes
    fread(&nT,   sizeof(unsigned short), 1, fp);
    fread(&nIso, sizeof(unsigned short), 1, fp);
    li->db[i].t  = nT;
    iso->db[i].i = nIso;

    //Allocate for the different temperature points and read
    T = li->db[i].T = (double *) calloc(nT, sizeof(double)  );
    fread(T, sizeof(double), nT, fp);

    //reallocate memory to account for new isotopes
    li->isov    = (prop_isov *)realloc(li->isov,
				       (correliso+nIso)*sizeof(prop_isov));
    iso->isof   = (prop_isof *)realloc(iso->isof, 
				       (correliso+nIso)*sizeof(prop_isof));
    iso->isov   = (prop_isov *)realloc(iso->isov, 
				       (correliso+nIso)*sizeof(prop_isov));
    iso->isodo  = (enum isodo *)realloc(iso->isodo,
					(correliso+nIso)*sizeof(enum isodo));
    li->isov[correliso].z = (double *)calloc((correliso+nIso)*nT,
					     sizeof(double));
    li->isov[correliso].c = (double *)calloc((correliso+nIso)*nT,
					     sizeof(double));
    transitDEBUG(21,verblevel,
		 "So far... CumIsotopes:%i, at databases: %i, position %li\n"
		 ,correliso+nIso,i,ftell(fp)); 

    //Reading isotopes from this database
    transitDEBUG(23,verblevel,
		 "DB %i: \"%s\"\n"
		 "  it has %i (%i) temperatures, %i (%i) isotopes,\n"
		 "  and starts at cumulative isotope %i\n"
		 ,iso->isof[correliso].d, iso->db[i].n, li->db[i].t,nT, iso->db[i].i,nIso,iso->db[i].s);

    for (unsigned int j=0 ; j<nIso ; j++){
      transitDEBUG(22,verblevel,"isotope %i/%i for DB %i\n",j+1,nIso,i);

      //initialize to be modified by getatm() and store DB number
      iso->isodo[correliso]  = unclear;
      iso->isof[correliso].d = i;

      fread(&rs, sizeof(unsigned short), 1, fp);
      iso->isof[correliso].n = (char *)calloc(rs+1, sizeof(char));
      fread(iso->isof[correliso].n, sizeof(char), rs, fp);
      iso->isof[correliso].n[rs] = '\0';
      transitDEBUG(21, verblevel,
		   "  Name's length: %i,"
		   "position: %li, value: %s\n"
		   ,rs,(long)(ftell(fp)),iso->isof[i].n);
      
      //read mass
      fread(&iso->isof[correliso].m, sizeof(double), 1, fp);

      transitDEBUG(21,verblevel,
		   "  Mass read: %g * %g = %g, "
		   "position: %li, size %i\n"
		   ,iso->isof[i].m,AMU,iso->isof[i].m*AMU
		   ,(long)(ftell(fp)),(int)(sizeof(iso->isof[i].m)));

      //read partition function  and cross section
      Z = li->isov[correliso].z = li->isov[correliso-j].z+nT*j;
      fread(Z, sizeof(double), nT, fp);
      CS = li->isov[correliso].c = li->isov[correliso-j].c+nT*j;
      fread(CS, sizeof(double), nT, fp);
      li->isov[correliso].n = nT;

      transitDEBUG(12, verblevel,
		   "Z(%i/%i):%g %g ... %g\n", j+1, nIso, Z[0], Z[1], Z[nT-1]);

      correliso++;
    }

    //Update acumulated isotope count (right before current DB)
    iso->db[i].s=acumiso;
    acumiso+=nIso;

    fread(&rs, sizeof(unsigned short), 1, fp);
    if (i!=rs)
      transiterror(TERR_SERIOUS,
		   "Problem in TLI file: database correlative number (%i)"
		   " doesn't match information read (%i)\n"
		   "Isotopes read: %i\n"
		   "Last DB #temps: %i\n"
		   "Last DB #iso: %i\n"
		   , i, rs, acumiso,nT,nIso);
  }
  //read total number of isotopes.
  fread(&iso->n_i,sizeof(unsigned short),1,fp);
  if(iso->n_i!=acumiso)
    transiterror(TERR_SERIOUS,
		  "Given number of isotopes (%i), doesn't equal\n"
		  "real total number of isotopes (%i)\n"
		  ,iso->n_i, acumiso);

  //update structure values
  li->ni=iso->n_i;
  li->ndb=ndb;
  iso->n_db=ndb;
  li->endinfo=ftell(fp);

  li->wi=iniw;
  li->wf=finw;

  return 0;
}


/* \fcnfh
   Read an TLI-ASCII formated file from an already open file 'fp'
*/
int 
readtli_ascii(FILE *fp, 
	       struct transit *tr,
	       struct lineinfo *li)
{
  char rc;
  char line[maxline+1],*lp,*lp2;
  int ndb,db;
  unsigned int nIso,nT;
  int acumiso;
  int rn,i;
  prop_isov *isov;
  PREC_ZREC *T;
  struct isotopes *iso=tr->ds.iso;

  //Format of the TLI-ASCII file is the following(names should not
  //contain spaces ('_' are replaced by spaces)
  //\begin{verb}
  //<m-database>
  //<DATABASE1-name> <n1-iso> <nt1-temp>
  //<NAME1> <MASS1> ... <NAMEn1> <MASSn1>
  //<TEMP1>    <Z1-1>  ...  <Z1-n1>   <CS1-1>  ...  <CS1-n1>
  //...
  //<TEMPnt1> <Znt1-1> ... <Znt1-n1> <CSnt1-1> ... <CSnt1-n1>
  //<DATABASE2-name> <n2-iso> <nt2-temp>
  //<NAME1> <MASS1> ... <NAMEn2> <MASSn2>
  //<MASS1> ... <MASSn2>
  //<TEMP1>    <Z1-1>  ...  <Z1-n2>   <CS1-1>  ...  <CS1-n2>
  //...
  //<TEMPnt2> <Znt2-1> ... <Znt2-n2> <CSnt2-1> ... <CSnt2-n2>
  //....
  //....
  //<DATABASEm-name> <nm-iso> <ntm-temp>
  //<NAME1> <MASS1> ... <NAMEnm> <MASSnm>
  //<TEMP1>    <Z1-1>  ...  <Z1-nm>   <CS1-1>  ...  <CS1-nm>
  //...
  //<TEMPntm> <Zntm-1> ... <Zntm-nm> <CSntm-1> ... <CSntm-nm>
  //<CTRWAV1> <ISOID1> <LOWENER1> <LOGGF1>
  //<CTRWAV2> <ISOID2> <LOWENER2> <LOGGF2>
  //...
  //...
  //\end{verb}
  //get Number of databases from first line
  while((rc=fgetupto_err(line,maxline,fp,&linetoolong,tr->f_line,
			 li->asciiline++))=='#'||rc=='\n');
  if(!rc) notyet(li->asciiline,tr->f_line);
  ndb=strtol(line,&lp,0);
  fprintf(stderr,"%i",ndb);
  checkprepost(lp,errno&ERANGE,*lp==' '||*lp=='\t',*lp!='\0');
  //Allocate pointers according to the number of databases
  iso->db=(prop_db *)calloc(ndb,sizeof(prop_db));
  li->db=(prop_dbnoext *)calloc(ndb,sizeof(prop_dbnoext));

  //for each database
  for(db=0;db<ndb;db++){
    //get name
    while((rc=fgetupto_err(lp=line,maxline,fp,&linetoolong,tr->f_line,
			   li->asciiline++)) == '#' || rc == '\n');
    if(!rc) notyet(li->asciiline,tr->f_line);
    if((iso->db[db].n=readstr_sp_alloc(lp,&lp,'_'))==NULL)
      transitallocerror(0);

    //go to next field and get number of temperatures and isotopes
    checkprepost(lp,0,*lp==' '||*lp=='\t',*lp=='\0');
    rn=getnl(2,' ',lp,&nIso,&nT);
    checkprepost(lp,rn!=2,0,nIso==0||nT==0);
    li->db[db].t=nT;
    iso->db[db].i=nIso;

    //Update acumulated isotope count
    iso->db[db].s=iso->n_i;
    iso->n_i+=nIso;

    //allocate for variable and fixed isotope info as well as for
    //temperature points.
    T=li->db[db].T=(PREC_ZREC *) calloc(nT,sizeof(PREC_ZREC)  );

    //allocate structure that are going to receive the isotope
    //info. If it is not first database then just reallocate
    if(!db){
      isov=li->isov=(prop_isov *)calloc(iso->n_i,sizeof(prop_isov));
      iso->isof=(prop_isof *)calloc(iso->n_i,sizeof(prop_isof));
      iso->isov=(prop_isov *)calloc(iso->n_i,sizeof(prop_isov));
      iso->isodo=(enum isodo *)calloc(iso->n_i,sizeof(enum isodo));
    }
    else{
      isov=li->isov=(prop_isov *)realloc(li->isov,iso->n_i*sizeof(prop_isov));
      iso->isof=(prop_isof *)realloc(iso->isof,iso->n_i*sizeof(prop_isof));
      iso->isov=(prop_isov *)realloc(iso->isov,iso->n_i*sizeof(prop_isov));
      iso->isodo=(enum isodo *)realloc(iso->isodo,iso->n_i*sizeof(enum isodo));
    }

    //Allocate cross section and temperature for this database. set
    //isov at the beggining of this database
    isov+=(acumiso=iso->db[db].s);
    isov->z=(PREC_ZREC *)calloc(nIso*nT,sizeof(PREC_ZREC));
    isov->c=(PREC_CS *)calloc(nIso*nT,sizeof(PREC_CS));

    //get isotope name and mass
    while((rc=fgetupto_err(lp2=line,maxline,fp,&linetoolong,tr->f_line,
			   li->asciiline++)) == '#' || rc == '\n');
    if(!rc) notyet(li->asciiline,tr->f_line);
    //for each isotope
    for(i=0;i<nIso;i++){
      isov[i].z=isov->z+nT*i;
      isov[i].c=isov->c+nT*i;
      //get name
      if((iso->isof[acumiso+i].n=readstr_sp_alloc(lp2,&lp,'_'))==NULL)
	transitallocerror(0);
      //get mass and convert to cgs
      iso->isof[acumiso+i].m=strtod(lp,&lp2);

      if(i!=nIso-1)
	checkprepost(lp2,lp==lp2,*lp2==' '||*lp2=='\t',*lp2=='\0');
    }
    //Last isotope has to be followed by an end of string
    checkprepost(lp2,0,*lp2==' '||*lp2=='\t',*lp2!='\0');

    //get for each temperature
    for(rn=0;rn<nT;rn++){
      //Get a line with temperature, partfcn and cross-sect info
      while((rc=fgetupto_err(lp=line,maxline,fp,&linetoolong,tr->f_line,
			     li->asciiline++)) == '#' || rc == '\n');
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
  iso->n_db=ndb;
  li->endinfo=ftell(fp);

  //We don't know beforehand what is the range in ascii storage, it has
  //to be looked for
  getinifinasctli(&li->wi,&li->wf,fp, tr->f_line);

  return 0;
}


/* \fcnfh
   Find initial and final wavelength from file pointed by 'fp'. who is
   at the end of the non-isotope info.

   @returns 0 on success
*/
int
getinifinasctli(double *ini,	/* where initial value would be stored */
		 double *fin,	/* where initial value would be stored */
		 FILE *fp,	/* file pointer */
		 char *file)	/* File name of file being read */
{
  char rc;
  char line[maxline+1],*lp,*lplast;

  //get first line
  while((rc=fgetupto_err(lp=line,maxline,fp,&linetoolong,file,0))
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
		 "readlineinfo:: weird, TLI-Ascii file has less than %i bytes.\n"
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
					  readinfo\_tli */ 
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
  PREC_LNDATA dbini=li->wi*TLI_WAV_UNITS,dbfin=li->wf*TLI_WAV_UNITS;
  double extra;

  //initialize modified hints
  msamp->n=-1;
  msamp->d=-1;
  msamp->v=NULL;
  msamp->fct=1;

  //Check that the factor is positive non-zero
  if(hsamp->fct<0)
    transiterror(TERR_SERIOUS,
		 "User specified wavelength factor is negative (%g)"
		 ,hsamp->fct);
  if(hsamp->fct>0)
    msamp->fct=hsamp->fct;

  double fct=msamp->fct;
  double fct_to_microns = msamp->fct/1e-4;

  //Check that the margin value is reasonable. i.e. whether its
  //leaves a non-zero range if applied to the whole line dataset.
  if(2*th->margin*fct > (dbfin - dbini)){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Margin value (%g microns) is too big for this dataset whose\n"
		 "range is %g to %g microns.\n"
		 "Factor to convert user margin (%g) to centimeters is %g\n"
		 ,th->margin*fct_to_microns
		 ,li->wi*tli_to_microns, li->wf*tli_to_microns
		 ,th->margin, msamp->fct);
    return -4;
  }
  margin=tr->margin=th->margin*msamp->fct;
 
  //If an TLI ascii file, then we'll be using limits twice the margin
  //from the minimum or maximum values in the dataset. This is because,
  //as opposed to binary archiving, I'm not guaranteed that I have all
  //the lines between 'dbini' and 'dbfin', instead they are the minimum
  //and maximum wavelength of the given transitions.
  if(li->asciiline){
    if(margin==0.0)
      transiterror(TERR_WARNING,
		   "Wavelength margin to be used is zero in a TLI-ASCII\n"
		   " file. Hence, there will be no points to the left or\n"
		   " right of the extreme central wavelengths.\n"
		   );
    extra=2*margin;
  }
  else
    extra=0;

  transitDEBUG(21,verblevel,
	       "in cgs: Hinted initial %g, final %g\n"
	       "Databse max %g and min %g\n"
	       ,hsamp->i*fct,hsamp->f*fct,dbini,dbfin);

  //If final wavelength was not hinted correctly then default it to zero
  if(hsamp->f<0){
    hsamp->f=0;
    transiterror(TERR_WARNING,
		 "Setting hinted upper wavelength limit before\n"
		 "extraction as %g. It was not user-hinted.\n"
		 ,hsamp->f*fct);
  }
  //Check final wavelength. If it is 0, modify it. Do not return
  //special value, it is assumed that user knows what he is doing.
  //\linelabel{finalwav}
  if(hsamp->f<=0){
      msamp->f=(dbfin + extra)/fct;
  }
  //otherwise
  else{
    transitDEBUG(20,verblevel,
		 "dbini: %g   margin:%g  sampf:%g\n"
		 ,dbini,margin,hsamp->f);

    //check that is not below the minimum value
    if( dbini + margin > fct*hsamp->f ){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "Considering margin, final wavelength (%g) is\n"
		   " smaller than first allowed value in database\n"
		   " (%g = %g + %g)\n"
		   ,hsamp->f*fct,dbini+margin
		   ,dbini,margin);
      return -3;
    }
    //warn if it is above maximum value with information
    if( (hsamp->f)+margin > dbfin )
      transiterror(TERR_WARNING,
		   "Final requested wavelength (%g microns) is above\n"
		   " the maximum informative\n"
		   "value in database (%g microns)\n"
		   ,hsamp->f*hsamp->fct*1e4,dbfin*hsamp->fct*1e4);
    //set initial wavelength to be extracted
    msamp->f=hsamp->f;
  }

  if(hsamp->i<0){
    hsamp->i=0;
    transiterror(TERR_WARNING,
		 "Setting hinted lower wavelength limit before\n"
		 "extraction as %g cgs. It was not user-hinted.\n"
		 ,hsamp->i*fct);
  }
  //Check initial wavelength. See final wavelength treatment above for
  //more detailed comments about the code (\lin{finalwav})
  if(hsamp->i<=0)
    msamp->i = (dbini - extra)/fct;
  else{
    transitDEBUG(20,verblevel,
		 "dbfin: %g   margin:%g  sampi:%g\n"
		 ,dbfin,margin,fct*hsamp->i);
    if( dbfin < margin+fct*hsamp->i ){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "Initial wavelength (%g) is larger than larger\n"
		   "allowed value in database (%g = %g + %g),\n"
		   "note that current margin was considered\n"
		   ,fct*hsamp->i,dbfin-margin
		   ,dbfin,margin);
      return -2;
    }
    if( fct*hsamp->i-margin < dbini )
      transiterror(TERR_WARNING,
		   "Initial requested wavelength (%g microns) is below\n"
		   " the minimum informative\n"
		   "value in database (%g microns)\n"
		   ,hsamp->i*1e4,dbini*1e4);
    msamp->i=hsamp->i;
  }

  //Check that we still have a range considering margin
  if(2*margin > (msamp->f-msamp->i)*fct){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Usable final (%g) has to be larger than usable\n"
		 "initial wavelength (%g). Note that those values\n"
		 "could have been modified according to the\n"
		 "database range (%g - %g) and margin (%g) -\n"
		 ,fct*msamp->i+margin
		 ,fct*msamp->f-margin
		 ,dbini,dbfin
		 ,margin);
    return -1;
  }

  //set progress indicator and return status
  tr->pi|=TRPI_CHKRNG;
  return res;
}



/* \fcnfh
   It outputs error. Used when EOF is found before expected
*/
static void
notyet(int lin, char *file)
{
  transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
	       "readlineinfo:: EOF unexpectedly found at line %i in\n"
	       "ascii-TLI linedb info file '%s'\n"
	       ,lin,file);
  exit(EXIT_FAILURE);
}



/* \fcnfh
  readinfo\_tli: Read infofile as returned by lineread.

  @todo    checks on allocation errors
  @returns 1 on success
           -1 unavailable file
	   -2 Filename not hinted
	   -3 TLI format not valid (missing magic bytes)
	   -4 Improper TLI-ASCII input
*/
int readinfo_tli(struct transit *tr,
		  struct lineinfo *li)
{
  //'fp' file pointer of info file.
  int rn;
  FILE *fp;

  union {char sig[4];int32_t s[2];} sign=
    {.s={0,((0xff-'T')<<24)|((0xff-'L')<<16)|((0xff-'I')<<8)|(0xff)}};
  char line[maxline+1];
  //Auxiliary hint structure pointer
  struct transithint *th=tr->ds.th;

  //Open info file and transport name into transit from hint (setting
  //corresponding flag) if the file opened succesfully
  if(!th->f_line){
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
  tr->f_line=th->f_line;

  //Read first four bytes, they should be either
  //`(0xff-T)(0xff-L)(0xff-I)(0xff)' or '\#TLI'. They are stored as
  //integer, so this check also serves to check whether the machine were
  //the data file and the one this program is being run have the same
  //endian order. If the first two are '\#TLI', then the first line
  //should also start as '\#TLI-ascii' 
  fread(sign.s,sizeof(int32_t),1,fp);

  //is it a binary TLI? (asciiline=0 indicates binary, asciiline=1
  //indicates ASCII TLI)
  li->asciiline=0;
  transitDEBUG(13,verblevel,
	       "Comparing %i and %i for Magic Number (len: %li)\n"
	       ,sign.s[0],sign.s[1],sizeof(sign.s[0]));
  
  if(sign.s[0]!=sign.s[1]){
    //does it look like being an Ascii TLI?, if so check it.
    rn=strncasecmp(sign.sig,"#TLI",4);
    if(!rn){
      strcpy(line,"#TLI");
      fread(line+4,sizeof(char),6,fp);
      rn=strncasecmp(line,"#TLI-ascii",10);
    }
    //If it wasn't any valid TLI, then error and exit
    if(rn){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		   "The file '%s' has not a valid TLI\n"
		   " format. It may also be because the machine were the\n"
		   " file was created have different endian order, which is\n"
		   " incompatible"
		   ,tr->f_line);
      return -3;
    }
    li->asciiline=1;
    //ignore the rest of the first line
    fgetupto_err(line,maxline,fp,&linetoolong,tr->f_line,1);
  }

  //if ascii format of TLI read ascii file
  if(li->asciiline){
    if((rn=readtli_ascii(fp,tr,li))!=0){
      transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
		   "readtli_ascii() return error code %i\n"
		   ,rn);
      return -5;
    }
  }
  //If binary storage is used.
  else
    if((rn=readtli_bin(fp,tr,li))!=0){
      transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
		   "readtli_bin() return error code %i\n"
		   ,rn);
      return -6;
    }

  transitprint(3,verblevel,
	       "TLI file read from %g to %g microns\n"
	       , li->wi, li->wf);
  

  //Wavelength in a TLI file is always in microns (as of TLI v4) and
  //lower energy is always in cm-1.
  struct line_transition *lt=&li->lt;
  lt->wfct=TLI_WAV_UNITS;
  lt->efct=TLI_E_UNITS;

  //close file, set progres indicator and return success.
  fclose(fp);
  tr->pi|=TRPI_READINFO;
  return 1;
}



/* \fcnfh
   print out an error, it is called by readdatarng if one of the field
   with transition info is invalid

   @returns -5 always
*/
static int
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



/*\fcnfh
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
					readinfo\_tli */ 
{
  //'fp' datafile file pointer.
  //'alloc' number of allocated line\_transition structures.
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
  PREC_LNDATA *ltgf;
  PREC_LNDATA *ltelow;
  PREC_LNDATA *ltwl;
  short *ltisoid;

  PREC_LNDATA iniw=li->wavs.i*li->wavs.fct/TLI_WAV_UNITS;
  PREC_LNDATA finw=li->wavs.f*li->wavs.fct/TLI_WAV_UNITS;
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
  ltgf=(PREC_LNDATA *)calloc(alloc,sizeof(PREC_LNDATA));
  ltwl=(PREC_LNDATA *)calloc(alloc,sizeof(PREC_LNDATA));
  ltelow=(PREC_LNDATA *)calloc(alloc,sizeof(PREC_LNDATA));
  ltisoid=(short *)calloc(alloc,sizeof(short));

  if(!ltgf||!ltwl||!ltelow||!ltisoid)
    transiterror(TERR_CRITICAL|TERR_ALLOC,
		 "Cannot allocate memory for linetran\n"
		 "structure array of length %i, in function readdatarng()\n"
		 ,alloc);

  //find starting point in datafile.\par
  //if it is TLI-ascii file then do a sequential search
  if(li->asciiline){
    //go to where readinfo\_tli left off and skip all comments, since
    //then, 'li->asciiline' won't increase
    fseek(fp,li->endinfo,SEEK_SET);
    while((rc=fgetupto_err(lp=line,maxline,fp,&linetoolong,tr->f_line,
			   li->asciiline++)) =='#'||rc=='\n');
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
      while((rc=fgetupto_err(lp=line,maxline,fp,&linetoolong,tr->f_line,
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
    rn=sizeof(short)+3*sizeof(PREC_LNDATA);
    nfields=((j-offs)/rn);
    if(nfields*rn+offs!=j){
      transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
		   "Data file does not have an integer number of records\n"
		   "Initial byte %i, final %i, record size %i\n"
		   ,offs,j,rn);
      return -3;
    }

    //Do the binary search in units of TLI
    datafileBS(fp, offs, nfields, iniw, &j, rn);
    transitDEBUG(21,verblevel,"Beginning found at position %li ",j);
    //check whether we need to start reading from records further back
    //because of repetition of wavelength
    if (j){
      do{
	fseek(fp,offs+--j*rn,SEEK_SET);
	fread(&wltmp,sizeof(PREC_LNDATA),1,fp);
      }while(wltmp>=iniw);
      j++;
    }
    //seek file to starting point.
    transitDEBUG(21,verblevel,"and then slide to %li\n",j);
    fseek(fp,offs+j*rn,SEEK_SET);
  }

  /* Main loop to read all the data. */
  //read the data in the main loop
  i=0;
  while(1){
    rn=rc=1;
    //enlarge structure if allocated spaced got filled up
    if(i==alloc){
      alloc<<=1;
      ltisoid=(short *)      realloc(ltisoid,alloc*sizeof(short)      );
      ltgf=   (PREC_LNDATA *)realloc(ltgf   ,alloc*sizeof(PREC_LNDATA));
      ltwl=   (PREC_LNDATA *)realloc(ltwl   ,alloc*sizeof(PREC_LNDATA));
      ltelow= (PREC_LNDATA *)realloc(ltelow ,alloc*sizeof(PREC_LNDATA));

      if(!ltgf||!ltwl||!ltelow||!ltisoid)
	transiterror(TERR_CRITICAL|TERR_ALLOC,
		     "Cannot enlarge memory allocation area\n"
		     "for array of linetran structure to length %i in\n"
		     "function readdatarng\n",alloc);
    }
    //if ascii, after skipping comments read the 4 fields: center, isoid,
    //lowE, loggf
    if(li->asciiline){
      while((rc=fgetupto_err(lp=line,maxline,fp,&linetoolong,
			     tr->f_line,li->asciiline+offs++))
	    =='#'||rc=='\n');
      //if it is not end of file, read the records.
      if(rc){
	ltwl[i]=strtod(lp,&lp2);
	if(lp==lp2) 
	  return invalidfield(line,tr->f_line,li->asciiline+offs
			      ,1,"central wavelength");
	ltisoid[i]=strtol(lp2,&lp,0);
	if(lp==lp2)
	  return invalidfield(line,tr->f_line,li->asciiline+offs
			      ,2,"isotope ID");
	ltelow[i]=strtod(lp,&lp2);
	if(lp==lp2)
	  return invalidfield(line,tr->f_line,li->asciiline+offs
			      ,3,"lower energy level");
	ltgf[i]=strtod(lp2,&lp);
	if(lp==lp2)
	  return invalidfield(line,tr->f_line,li->asciiline+offs
			      ,4,"log(gf)");
      }
    }
    //If binary read a data structure
    else{
      rn=fread(ltwl+i,sizeof(PREC_LNDATA),1,fp);
      fread(ltisoid+i,sizeof(short),1,fp);
      fread(ltelow+i,sizeof(PREC_LNDATA),1,fp);
      fread(ltgf+i,sizeof(PREC_LNDATA),1,fp);
      transitDEBUG(26,verblevel,"Wavelength:%.8f iso:%i\n",ltwl[i],
		   ltisoid[i]);
    }

    //Warn if EOF was found. This should be OK if you want to read
    //until the very last line.
    if(!rn||!rc){
      transiterror(TERR_WARNING,
		   "End-of-file in datafile '%s'.\n"
		   "Last wavelength read (%f) was in record %i.\n"
		   "If you are reading the whole range, you can safely ignore\n"
		   "this warning.\n\n"
		   ,tr->f_line,ltwl[i-1],i);
      break;
    }

    /* TD: Skip isotope at read time */

    //End if we passed maximum required wavelength
    if(ltwl[i]>finw)
      break;

    i++;
  }
  transitDEBUG(21,verblevel,
	       "Number of lines just read: %li\n"
	       ,i);

  //Realloc final strucuture array to its required size and set the
  //pointer in the transit structure.
  alloc=i;
  struct line_transition *lt=&li->lt;
  lt->isoid=(short *)      realloc(ltisoid,alloc*sizeof(short)      );
  lt->gf   =(PREC_LNDATA *)realloc(ltgf   ,alloc*sizeof(PREC_LNDATA));
  lt->wl   =(PREC_LNDATA *)realloc(ltwl   ,alloc*sizeof(PREC_LNDATA));
  lt->elow =(PREC_LNDATA *)realloc(ltelow ,alloc*sizeof(PREC_LNDATA));

  if(!ltgf||!ltwl||!ltelow||!ltisoid){
    transiterror(TERR_CRITICAL|TERR_ALLOC,
		 "Cannot diminish memory allocation area\n"
		 "for array of linetran structure to length %i in\n"
		 "function readdatarng\n",alloc);
    return -1;
  }  

  //store number of lines
  li->n_l=i;

  //close file, set progress indicator and exit the number of read line.
  fclose(fp);
  tr->pi|=TRPI_READDATA;
  return i;
}



/* TD: Accept isotope list */
/* TD: Error return codes */
/* TD: check that strtok doesn't overwrite first string */
/* \fcnfh
   readlineinfo: It only function is to call all the functions regarding
   reading TLI information.

   @returns 0 on success
*/
int readlineinfo(struct transit *tr) /* General parameters and
					hints  */
{
  long rn;
  struct transithint *th=tr->ds.th;
  static struct lineinfo st_li;
  static struct isotopes st_iso;
  memset(&st_li,0,sizeof(struct lineinfo));
  memset(&st_iso,0,sizeof(struct isotopes));
  tr->ds.li=&st_li;
  tr->ds.iso=&st_iso;

  //Try to read hinted info file
  transitprint(1,verblevel, "\nReading info file '%s'...\n"
	       ,th->f_line);
  if((rn= readinfo_tli(tr,&st_li))!=1)
    transiterror(TERR_SERIOUS,
		 "readinfo_tli() returned an error code %i!\n"
		 ,rn);
  transitprint(1,verblevel, " done%c\n",'.');

  //Check the remainder (margin and range) of the hinted values
  //related to line database reading.
  if((rn=checkrange(tr,&st_li))<0)
    transiterror(TERR_SERIOUS,
		 "checkrange() returned error code %i!\n"
		 ,rn);
  //output status so far if the verbose level is enough
  if(rn>0&&verblevel>1)
    transiterror(TERR_WARNING,
		 "checkrange() modified the suggested parameters,\n"
		 "it returned code 0x%x\n\n"
		 ,rn);
  double fct=st_li.wavs.fct;
  double fct_to_microns = fct/1e-4;

  transitprint(2,verblevel,
	       "   After checking limits, the wavelength range to be\n"
	       "  used is %g to %g, including a margin of %g.\n"
	       ,fct*tr->ds.li->wavs.i,fct*tr->ds.li->wavs.f
	       ,tr->margin);

  //read data file
  transitprint(1,verblevel, "\nReading data...\n");
  if((rn=readdatarng(tr,&st_li))<1)
    transiterror(TERR_SERIOUS,
		 "readdatarng() returned an error code %li\n"
		 ,rn);
  transitprint(1,verblevel, " done.\n");

  //Status so far
  transitprint(2,verblevel,
	       "OK, So far:\n"
	       " * I read %li records from the datafile.\n"
	       " * The wavelength range read was %.8g to %.8g microns.\n"
	       " * Current margin is %.4g microns.\n"
	       " * Usable range is thus %.8g to %.8g.\n"
	       , st_li.n_l
	       , st_li.wavs.i*fct_to_microns,st_li.wavs.f*fct_to_microns
	       , tr->margin*1e4
	       , st_li.wavs.i*fct_to_microns+tr->margin
	       , st_li.wavs.f*fct_to_microns-tr->margin);



#ifndef NODEBUG_TRANSIT
  rn=1; //Some random number to test
  struct line_transition *lt=&tr->ds.li->lt;
  transitDEBUG(21,verblevel,
	       " * And the record %li has the following info\n"
	       "Wavelength: %.10g\n"
	       "Lower Energy Level: %.10g\n"
	       "Log(gf): %.10g\n"
	       "Isotope: %i\n"
	       ,rn,lt->wl[rn],lt->elow[rn], lt->gf[rn]
	       ,lt->isoid[rn]);
#endif

  transitDEBUG(21,verblevel,
	       "database min and max: %.10g(%.10g) and %.10g(%.10g)\n"
	       ,st_li.wi,tr->ds.li->wi
	       ,st_li.wf,tr->ds.li->wf);

  return 0;
}



/* \fcnfh
   frees lineinfo structure 

   @returns 0 on success
*/
int
freemem_isotopes(struct isotopes *iso,
		 long *pi)
{
  int i;

  //free structures
  for(i=0;i<iso->n_i;i++)	/* Allocated in readlineinfo */
    free_isof(iso->isof+i);
  free_isof(iso->isof+iso->n_i); /* Allocated in at_file */
  for(i=0;i<iso->n_db;i++)
    free_db(iso->db+i);
  free_isov(iso->isov);

  //free arrays
  free(iso->isodo);
  free(iso->isof);
  free(iso->isov);
  free(iso->db);

  //unset appropiate flags
  *pi&=~(TRPI_READINFO|TRPI_READDATA|TRPI_CHKRNG|
	 TRPI_GETATM);
  return 0;
}



/* \fcnfh
   frees lineinfo structure 

   @returns 0 on success
*/
int
freemem_lineinfotrans(struct lineinfo *li,
		   long *pi)
{
  int i;

  //free the four arrays of lt
  struct line_transition *lt=&li->lt;
  free(lt->wl);
  free(lt->elow);
  free(lt->gf);
  free(lt->isoid);

  //free isov, dbnoext and samp in li
  free_isov(li->isov);
  free(li->isov);

  for(i=0;i<li->ndb;i++)
    free_dbnoext(li->db+i);
  free(li->db);

  free_samp(&li->wavs);

  //zero all the structure
  memset(li,0,sizeof(struct lineinfo));

  //unset appropiate flags.
  *pi&=~(TRPI_READDATA|TRPI_READINFO|TRPI_CHKRNG);
  return 0;
}


/* \fcnfh
   Saves line information
*/
void
saveline(FILE *fp,
	 struct lineinfo *li)
{
}


#ifdef DBGREADLINEINFO
/* \fcnfh
   main function for debugging only
*/
int main(int argc, char **argv)
{
  struct transit tr;
  struct transithint th;
  struct lineinfo *li;
  int i,ti1,nbins,ans;
  PREC_LNDATA *ltgf;
  PREC_LNDATA *ltelow;
  PREC_LNDATA *ltwl,*twl;
  short *ltisoid,*tisoid;

  tr.ds.th=&th;
  th.na=0;

  verblevel=20;

  th.m=0.001;
  th.na|=TRH_WM;
  char defile_line[]="./res/lineread.tli";
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
  ltgf=tr.ds->lt.gf;
  ltwl=tr.ds->lt.wl;
  ltisoid=tr.ds->lt.isoid;
  ltelow=tr.ds->lt.elow;

  ti1=(int)(log10(li->n_l)+1);

  printf("Done reading the file.\n\n"
	 "dbread_pands() test results:\n");
  printf("Chosen wavelength range was from %.10g to %.2f [nm]\n"
	 " %*li lines read\n"
	 " Choosing %i equal-sized bins, the result is\n"
	 ,li->wi,li->wf,ti1,li->n_l,nbins);

  long qb[tr.n_i];
  float szb=(li->wf-li->wi)/nbins;
  double endb;
 
  twl=ltwl;
  tisoid=ltisoid;
  if(!nbins)
    Pprintf(1,"  hmmm, you chose 0 bins!\n");
  for(i=0;i<nbins;i++){
    memset(qb,0,sizeof(*qb)*4);
    endb=li->wi+(i+1)*szb;
    //    PPprintf(1,2,"KK %g %f\n",lp->wl,endb);
    while(*twl<endb&&twl-ltwl<li->n_l){
      qb[*tisoid++]++;
      twl++;
    }

    Pprintf(1," %*i = %i + %i + %i + %i lines shorter than %.3f\n"
	   ,ti1,qb[0]+qb[1]+qb[2]+qb[3],qb[0],qb[1],qb[2],qb[3],endb);
  }


  Pprintf(1,"\nWanna know the value of a single record?\n"
	  "If so, write record number (range 0 - %i), else "
	    "press ^C: "
	  ,li->n_l-1);

  while(Pgeti(0,&ans,(int)(log10(li->n_l))+1)>=0){
    if(ans<li->n_l&&ans>=0){
      Pprintf(1,"Wavelength: %.10g\n"
	      ,ltwl[ans]);
      Pprintf(1,"Lower Energy Level: %.10g\nLog(gf): %.10g\n"
	      , ltelow[ans], ltgf[ans]);
      printf("Isotope Name: %s\n"
	     ,tr.isof[ltisoid[ans]].n);
    }
    else
      Pprintf(1,"\nInvalid record number, so ...");

    Pprintf(1,"\nWanna know the value of another single record?\n"
	    "If so, write the record number (range 0 - %i), else just "
	    "press ^C: "
	    ,li->n_l-1);
  }
  
}

#undef checkprepost

#endif

