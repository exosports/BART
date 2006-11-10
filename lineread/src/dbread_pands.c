/*
 * dbread_pands.c - Driver to read Partridge & Schwenke for lineread.
 *              Part of lineread program.
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
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include <lineread.h>
#include <math.h>


#if 0
/* TD: Is not safe to use structre to read data, because it may
   contain unexpected padding */
/* TD: BS can be fooled by lines with identical line centers */


#define PANDS_RECLENGTH 8
#define PANDS_NCODIDX 32786
#define DRIVERNAME "pands"

/* Wavelength is stored in nanometers */
static double pands_fct=1e-7;
static const char *pands_fct_ac="nm";
short gabby_dbread;

static int isoname(char ***isotope, int niso);
static int read_zpands(char *filename, PREC_ZREC ***Z, PREC_ZREC **T,
		       PREC_CS ***CS, int *nT, int nIso);



/*
  reversebytes: Invert the byte order. It can be used for any variable
  type, just specify the number of bytes of its class through the use of
  sizeof().

  @requires pointer to be a pointer to any type and byte to be an
  unsigned int.
*/




PREC_NREC
bread_pands(char *filename,
	     struct linedb **lines, //2 pointers in order to be
				//able to allocate memory
	     float wlbeg,           //wavelengths in tli_fct
	     float wlend,           //units
	     /* Partition function data file */
	     char *Zfilename,
	     /* For the following 3 parameter, the memory is
		allocated in the dbread_* functions, and the
		size is returned in the last parameters. */
	     PREC_ZREC ***Z,        //Partition function:[iso][T]
	     PREC_ZREC **isomass,   //Isotopes' mass in AMU
	     PREC_CS ***CS,         //Isotope's cross-section: [iso][T]
	     PREC_ZREC **T,         //temps for Z and CS
	     int *nT,               //number of temperature
				//points 
	     int *nIso,             //number of isotopes
	     char ***isonames)      //Isotope's name
{
  PREC_NREC nrec;
  char *deffname  = "./oth/pands/h2ofast.bin";
  char *defzfname = "./oth/pands/h2opartfn.dat";

  struct linedb *line;
  struct recordstruct{
    PREC_BREC iwl;                    //Wavelength index
    short int ielow,igflog;           //Lower energy and log(gf) indices.
  }record;
  struct stat fs;
  double tablog[PANDS_NCODIDX+1];     //Array used to convert stored int
				      //to float
  double lnwl;                        //log of beginning wavelength to
				      //be used in BS
  PREC_NREC lnwl1=0, lnwl2=0, irec, frec; //upper and lower values of BS
  double ratiolog;
  PREC_NREC i;
  FILE *fp;

  wlbeg*=pands_fct/tli_fct;
  wlend*=pands_fct/tli_fct;

  *nIso=NUM_ISOT;
  if(!filename)
    filename=deffname;
  if(!Zfilename)
    Zfilename=defzfname;

  /*
    Following is to transform stored integer of gflog and ielow into their
    real floating point values
  */
  for(i=1;i<=PANDS_NCODIDX;i++)
    tablog[i]=pow(10,(i-16384)*0.001);

  ratiolog=log((double)1.0+(double)1.0/(2e+6));

  if(stat(filename,&fs)==-1){
    mperror(MSGP_USER|MSGP_ALLOWCONT,
	    "Data file '%s' cannot be accesed by stat() in "
	    "function dbread_pands().\nThis is important to obtain "
	    "its size and hence the number of lines to be\n "
	    "examinated\n"
	    ,filename);
    return -2;
  }
  nrec=fs.st_size;
  if(nrec/PANDS_RECLENGTH<nrec/(float)PANDS_RECLENGTH){
    mperror(MSGP_USER|MSGP_ALLOWCONT,
	    "Data file '%s' does not contain an integer number of "
	    "%i-bytes records!.\nAre you sure it is the right %s "
	    "file?\n"
	    ,filename,PANDS_RECLENGTH,pands_name);
    return -3;
  }
  nrec/=PANDS_RECLENGTH;

  if((fp=fopen(filename,"r"))==NULL){
    mperror(MSGP_USER|MSGP_ALLOWCONT,
	    "Data file '%s' cannot be opened\nas stream in function "
	    "dbread_pands(), stopping.\n"
	    ,filename);
    return -1;
  }
  fread(&lnwl1,4,1,fp);
  reversebytes(&lnwl1,4);

  fseek(fp,-PANDS_RECLENGTH,SEEK_END);
  fread(&lnwl2,4,1,fp);
  reversebytes(&lnwl2,4);

  lnwl=log(wlbeg)/ratiolog;
  if(lnwl2>0 && lnwl>lnwl2){
    mperror(MSGP_USER|MSGP_ALLOWCONT,
	    "Last wavelength in the file (%g %s) is shorter "
	    "than\nrequested initial wavelength (%g %s)\n"
	    ,exp(lnwl2*ratiolog), pands_fct_ac, wlbeg
	    ,pands_fct_ac);
    return -4;
  }
  if(lnwl1>0 && lnwl<lnwl1){
    wlbeg=exp(ratiolog*lnwl1);
    irec=0;
  }
  else
    dbreadBSf(fp,0,nrec,lnwl,&irec,8);
    //    pandsBS(0,nrec,lnwl,&irec);
  if(gabby_dbread>1)
    fprintf(stderr,
	    "Located beginning wavelength %g at position %li\n", wlbeg, irec);

  lnwl=log(wlend)/ratiolog;
  if(lnwl1>0 && lnwl<lnwl1){
    mperror(MSGP_USER|MSGP_ALLOWCONT,
	    "First wavelength in the file (%g %s) is longer "
	    "than\nrequested ending wavelength (%g %s)\n"
	    ,exp(lnwl1*ratiolog), pands_fct_ac, wlend
	    ,pands_fct_ac);
    return -5;
  }
  if(lnwl2>0 && lnwl>lnwl2){
    wlend=exp(ratiolog*lnwl1);
    frec=0;
  }
  else{
    dbreadBSf(fp,irec,nrec,lnwl,&frec,8);
    //    pandsBS(irec,nrec,lnwl,&frec);
    frec++;
  }
  if(gabby_dbread>1){
    fprintf(stderr,
	    "Located ending wavelength %g at position %li\n",wlend,frec);
    fprintf(stderr,
	    "About to initialize memory space to hold %li records.\n"
	    " I'll require %.2fMb of available memory.\n"
	    ,frec-irec
	    ,(frec-irec)*(float)sizeof(struct linedb)/1024.0/1024.0);
  }

  if(((*lines)=calloc(frec-irec,sizeof(struct linedb)))==NULL){
    mperror(MSGP_SYSTEM|MSGP_ALLOWCONT,
	    "Cannot allocate memory to hold all the data from %s\n"
	    "Required memory was %li bytes\n"
	    ,pands_name,(frec-irec)*sizeof(struct linedb));
    return -6;
  }
  if(gabby_dbread>1)
    fprintf(stderr,"Success in memory allocation\n");

  /* Main reading loop */
  if(gabby_dbread>0)
    fprintf(stderr,"reading... ");

  fseek(fp,irec*PANDS_RECLENGTH,SEEK_SET);

  i=0;
  do{
    line=(*lines)+i;

    /*TD: Is not safe to use structre to read data, because it may
      contain padding */

    fread(&record,1,8,fp);
    reversebytes(&(record.iwl),4);
    reversebytes(&(record.ielow),2);
    reversebytes(&(record.igflog),2);

#if 0
    int b;
    fprintf(stderr,"\n%08li  ",irec+i);
    for(b=0;b<8;b++)
      fprintf(stderr,"%02hhx ",(char)*(((char *)&record)+b));
#endif

    line->wl=exp(record.iwl*ratiolog)*pands_fct/tli_fct;

    //Isotopes (1h1h16o, 1h1h17o, 1h1h18o, 1h2h16o)
    if(record.ielow>0)
      line->isoid=record.igflog>0?0:1; //isotopes indices are Kurucz's - 1
    else
      line->isoid=record.igflog>0?2:3;
    line->recpos=i+irec;
    line->elow=(PREC_LNDATA)abs(record.ielow);
    line->gf=tablog[abs(record.igflog)];

    i++;
  }while(line->wl<wlend && i+irec<frec);
  if(gabby_dbread>0)
    fprintf(stderr, "done\n");

  fclose(fp);

  if(isonames != NULL)
    isoname(isonames, *nIso);

  if(isomass != NULL){
  }

  if(Z!=NULL)
    if((nrec=read_zpands(Zfilename,Z,T,CS,nT,*nIso))!=1)
      mperror(MSGP_USER,
	      "Function read_zpands() return error code '%i'\n",
	      nrec);


  return i;
}


#define MAX_LINE 200
/*\fcnfh
  read_zpands: Read from file 'filename'(char *), the partition function
  information into 'Z'(PREC_ZREC ***), it is a two dimensional array
  depending on isotope (1st dimension) and temperature(2nd
  dimension). Values of temperature are stored in 'T'(PREC_ZREC **), of
  size 'nT'(int *). 'nIso'(int) indicates the number of isotopes in
  consideration. It is assumed that the isotopes columns have the same
  order as the indices being used.
*/
static int read_zpands(char *filename, /* Doh! */
		       PREC_ZREC ***Z, /* Partition function */
		       PREC_ZREC **T,  /* Temperature points */
		       PREC_CS  ***CS, /* Cross section points */
		       int *nT,        /* Number of temp. points */
		       int nIso)       /* Number of isotopes */
{
  FILE *fp;
  int i,j,cnt;
  char line[MAX_LINE], *sp,*sp2;
  int ignorelines;

  ignorelines=5;		/* Header */
  (*nT)=8; 			/* Initial value for allocation of
				   temperature info */

  if((fp=fopen(filename,"r"))==NULL){
    mperror(MSGP_USER,
	    "Data file '%s' cannot be opened as stream in function "
	    "read_zpands(), stopping.\n"
	    ,filename);
  }

  for(i=0;i<ignorelines;i++){
    fgets(line,MAX_LINE,fp);
  }

  *Z   = (PREC_ZREC **)calloc(nIso, sizeof(PREC_ZREC *));
  *CS  = (PREC_CS   **)calloc(nIso, sizeof(PREC_CS *));
  **CS = (PREC_CS    *)calloc(nIso*(*nT), sizeof(PREC_CS));
  for(i=0 ; i<nIso ; i++){
    (*Z)[i]  = (PREC_ZREC *)calloc((*nT), sizeof(PREC_ZREC));
    (*CS)[i] = (PREC_CS   *)calloc((*nT), sizeof(PREC_CS));
  }
  *T=(PREC_ZREC *)calloc((*nT),sizeof(PREC_ZREC));

  cnt=0;

  while(fgets(line,MAX_LINE,fp)!=NULL){
    sp=line;
    (*T)[cnt]=strtod(sp,&sp2);

    for(i=0;i<nIso;i++){
      if(sp2==sp){
	mperror(MSGP_USER,
		"In function read_zpands(): line %i of file\n '%s'"
		" has %i columns instead of %i.\n",
		cnt+ignorelines,filename, i+1,nIso+1);
      }
      sp = sp2;
      (*Z)[i][cnt] = strtod(sp,&sp2);
    }

    if(++cnt==(*nT)){
      (*nT)<<=1;
      (*T)=(PREC_ZREC *)realloc((*T),(*nT)*sizeof(PREC_ZREC));
      for(i=0;i<nIso;i++)
	(*Z)[i]=(PREC_ZREC *)realloc((*Z)[i],(*nT)*sizeof(PREC_ZREC));
    }
  }
  (*nT)=cnt;
  (*T)=(PREC_ZREC *)realloc((*T),(*nT)*sizeof(PREC_ZREC));
  for(i=0;i<nIso;i++){
    (*Z)[i]  = (PREC_ZREC *)realloc((*Z)[i], (*nT)*sizeof(PREC_ZREC));
    (*CS)[i] = (PREC_CS   *)realloc((*CS)[i],(*nT)*sizeof(PREC_CS));
    for(j=0 ; j<(*nT) ; j++){
      (*CS)[i][j] = SIGWATER;
    }
  }

  fclose(fp);

  return 1;

}
#undef MAX_LINE

/*
  isoname: returns the name of the isotopes in the newly allocated array
  'isoname'(char ***)

  @todo    error check for calloc calls
  @returns 1 on success
*/
static int isoname(char ***isonames, int niso)
{
  int i;

  *isonames=(char **)calloc(niso,sizeof(char *));
  for(i=0 ; i<niso ; i++){
    (*isonames)[i] = (char *)calloc(strlen(isotope[i])+1,
				    sizeof(char));
    strcpy((*isonames)[i], isotope[i]);
  }

  return 1;
  
}

#ifdef TEST_RUN

/* \fcnfh
   
 */
static int
test_fill_defaults(char *file,
		   char *zfile,
		   char *dbread_wl_unit,
		   int *numisot,
		   char ***isonames,
		   PREC_NREC (*dbread_main)(char *,  struct linedb **,
					    float, float, char *, 
					    PREC_ZREC ***, PREC_ZREC **, 
					    PREC_ZREC **, PREC_CS ***,
					    int *, int *, char ***))
{
  strcpy(file,"../oth/pands/h2ofast.bin");
  strcpy(zfile,"../oth/pands/h2opartfn.dat");
  strcpy(dbread_wl_unit, pands_fct_ac);
  *numisot=NUM_ISOT;
  *isonames = (char **)calloc(NUM_ISOT, sizeof(char *));
  **isonames = (char *)calloc(sizeof(isotope)+NUM_ISOT, sizeof(char));

  fprintf(stderr, "Allocating %i bytes for isotope. Correct and check.\n"
	  , sizeof(isotope)+1);

  int i;
  for (i=0 ; i<NUM_ISOT ; i++){
    int len=strlen(strcpy(*isonames[i], isotope[i]))+1;
    if (i!=NUM_ISOT-1)
      *isonames[i+1] = *isonames[i] + len;
  }

  dbread_main = &dbread_pands;
  return 0;
}

#endif
#endif









/**********************************
 **********************************
 **********************************/


#include <lineread.h>

#ifndef __BYTE_ORDER 
#error __BYTE_ORDER not defined, include <endian.h>
#endif
#if __BYTE_ORDER == __BIG_ENDIAN
static inline void
reversebytes(void *pointer, int nbytes)
{
  char one;
  int n,tot;
  tot = nbytes--/2;
  for(n=0 ; n<tot ; n++){
    one = *((char *)(pointer)+n);
    *((char *)(pointer)+n) = *((char *)(pointer)+nbytes-n);
    *((char *)(pointer)+nbytes-n) = one;
  }
}
#else
#define reversebytes(pointer, nbytes) do{}while(0)
#endif   /* __BYTE_ORDER */



#define NUM_ISOT 4
#define RWATER   (3.2e-8/2.0)	//water molecule radius (cm)
#define PI       (3.141592653589793) //PI!
#define SIGWATER (PI*RWATER*RWATER) //water cross section (cm^2)
#define PREC_BREC long int		/* Type for the Big record  */

#define PANDS_RECLENGTH 8
#define PANDS_NCODIDX   32786

static const char *defaux = "h2opartfn.dat";
static const char *deftarget = "h2ofast.bin";
static const char *pands_name="Partridge & Schwenke (1997)";
static const char     *pands_iso[] = {"1H1H16O",   "1H1H17O",   "1H1H18O",   "1H2H16O"};
static const PREC_MASS pands_mass[] = {18.01056468, 19.01478156, 20.01481046, 19.01684143};

/* Wavelength is stored in nanometers */
static double pands_fct=1e-7;
static const char *pands_fct_ac="nm";

static char *dbfilename = NULL, *dbauxname = NULL;
static FILE *fp = NULL, *fpaux = NULL;
static _Bool partitionread = 0;
static int verbose_db = 10;


#define sfread(pointer, size, nmemb, f) sfread_fcn(pointer, size, nmemb, f, __LINE__)
static void
sfread_fcn(void *pointer, size_t size, size_t nmemb, FILE *f, long int line)
{
  size_t read = fread(pointer, size,  nmemb, f);
  if (read!=nmemb)
    fprintf(stderr,
	    "%s (%li): Number of elements read was only %li of %li\n"
	    , __FILE__, line, (long)read, (long)nmemb);
}

/* \fcnfh
   free memory used by this driver
 */
static void
pands_free()
{
  if (dbauxname)
    free(dbauxname);
}


/* \fcnfh
  dbreadBSf: Perform a binary search in file pointed by 'fpbs'(FILE *)
  between 'initial'(PREC_NREC) and 'final'(PREC_NREC) looking for
  'lookfor'(PREC_NREC) at the first item of a record, result is stored
  in 'resultp'(PREC_NREC *). Records are of length 'reclength'(int) and
  the first item of each of them is of type PREC_BREC.
*/
static inline void
dbreadBSf(FILE *fpbs,		/* File pointer */
	  PREC_NREC initial,	/* initial index */
	  PREC_NREC final,	/* last index */
	  double lookfor,	/* target value */
	  PREC_NREC *resultp,	/* result index */
	  int reclength)	/* Total length of record */
{
  int irec1,irec2;
  PREC_BREC temp;
  const int trglength=sizeof(PREC_BREC);

  irec1=initial;
  irec2=final;
  do{
    *(resultp)=(irec2+irec1)/2;
    fseek(fpbs,reclength*(*resultp),SEEK_SET);
    sfread(&temp,trglength,1,fpbs);
    reversebytes(&temp,trglength);
    if(lookfor>temp)
      irec1=*(resultp);
    else
      irec2=*(resultp);
  }while (irec2-irec1>1);
  *resultp=irec1;
}


/* \fncfh
   If this driver can read name, then return 1 otherwise 0
*/
static _Bool
db_find(const char *name)
{
  int len = strlen(name), lent = strlen(deftarget);
  

  messagep(verbose_db, "Comparing given '%s' with stored '%s'\n"
	   , name, deftarget);
  if (len >= lent &&
      strncmp(deftarget, name+len-lent, lent) == 0)
    return 1;

  return 0;
}


/********************************************************/


/* \fcnfh
   Open database and auxiliary file, set all FILE pointers.
*/
static int
db_open(char *dbname, 
	char *dbaux)
{
  //Open database
  if((fp = fopen(dbname,"r")) == NULL)
    mperror(MSGP_USER,
	    "Could not open file '%s' for reading\n"
	    , dbname);
  dbfilename = dbname;

  if (dbaux)
    dbauxname = dbaux;
  else{
    int len = strlen(dbname), 
      lent = strlen(deftarget), 
      lena = strlen(defaux);
    dbauxname = (char *)calloc(len-lent+lena+2, sizeof(char));
    strncpy(dbauxname, dbname, len-lent);
    strcpy(dbauxname+len-lent, defaux);    
  }

  //Optionally open auxiliary
  if((fpaux = fopen(dbauxname,"r")) == NULL)
    mperror(MSGP_USER,
	    "Could not open auxiliary file '%s' for reading\n"
	    , dbauxname);

  return LR_OK;
}


/********************************************************/


/* \fcnfh
   Close all files
 */
static int
db_close()
{
  if(fp)
    fclose(fp);
  if(fpaux)
    fclose(fpaux);

  pands_free();

  return LR_OK;
}


/********************************************************/

/* \fcnfh
   Reads partition info 
*/
static _Bool
db_part(char **name,
	unsigned short *nT,
	PREC_TEMP **T,
	unsigned short *niso,
	char ***isonames,
	PREC_MASS **mass,
	PREC_Z ***Z,
	PREC_CS ***CS)
{

  long maxline = 300;
  char line[maxline], *lp, *lp2;
  int ignorelines=5;		/* Header */

  unsigned short nIso = *niso = NUM_ISOT;
  int nt=0, nta = 8;		/* Initial value for allocation of
				   temperature info */


  for(int i=0 ; i<ignorelines ; i++)
    fgets(line, maxline, fpaux);

  PREC_Z  **z     = *Z  = (PREC_Z  **)calloc(nIso, sizeof(PREC_Z *));
  PREC_CS **cs    = *CS = (PREC_CS **)calloc(nIso, sizeof(PREC_CS *));
  PREC_TEMP *temp =      (PREC_TEMP *)calloc(nta,  sizeof(PREC_TEMP));
  *z  = (PREC_Z    *)calloc(nIso*nta, sizeof(PREC_Z));
  *cs = (PREC_CS   *)calloc(nIso*nta, sizeof(PREC_CS));
  for(unsigned short i=1 ; i<nIso ; i++){
    z[i]  = *z  + i*nta;
    cs[i] = *cs + i*nta;
  }


  while((lp=fgets(line,maxline,fpaux)) != NULL){
    if(nt == nta){
      unsigned int old = nta;
      *z   = (PREC_Z    *)realloc(*z,   nIso*(nta<<=1)*sizeof(PREC_Z));
      *cs  = (PREC_CS   *)realloc(*cs,  nIso*nta      *sizeof(PREC_CS));
      temp = (PREC_TEMP *)realloc(temp, nta           *sizeof(PREC_TEMP));
      for(unsigned short i=nIso-1 ; i>0 ; i--){
	z[i]  = *z  + i*nta;
	cs[i] = *cs + i*nta;
	memmove(z[i],  *z  + i*old, old);
	memmove(cs[i], *cs + i*old, old);
      }
    }	


    temp[nt]=strtod(lp,&lp2);

    for(int i=0 ; i<nIso ; i++){
      if(lp2==lp){
	mperror(MSGP_USER|MSGP_ALLOWCONT,
		"In function read_zpands(): line %i of file\n '%s':\n%s"
		" has %i columns (isotopes+1) instead of %i.\n",
		nt+ignorelines, dbauxname, line, i+1, nIso+1);
	free(*z);
	free(*cs);
	free(z);
	free(cs);
	free(temp);
	lineread_free();
	exit(EXIT_FAILURE);
      }
      lp = lp2;
      z[i][nt]  = strtod(lp,&lp2);
      cs[i][nt] = SIGWATER;
    }
  }

  *nT=nt;
  *T = (PREC_TEMP *)realloc(temp, nt*sizeof(PREC_TEMP));
  for(unsigned short i=0 ; i<nIso ; i++){
    (*Z)[i]  = *z  + i*nt;
    (*CS)[i] = *cs + i*nt;
    memmove((*Z)[i],  *z  + i*nta, nt);
    memmove((*CS)[i], *cs + i*nta, nt);
  }
  **Z  = (PREC_Z  *)realloc(**Z,  nIso*nt*sizeof(PREC_Z));
  **CS = (PREC_CS *)realloc(**CS, nIso*nt*sizeof(PREC_CS));

  if(fpaux)
    fclose(fpaux);
  fpaux=NULL;

  *name     = strdup(pands_name);
  *isonames = (char     **)calloc(nIso, sizeof(char *));
  *mass     = (PREC_MASS *)calloc(nIso, sizeof(PREC_MASS));
  for (int i=0 ; i<nIso ; i++){
    (*mass)[i] = pands_mass[i];
    (*isonames)[i] = strdup(pands_iso[i]);
  }

  partitionread = 1;
  return 0;
}


/********************************************************/


/* \fcnfh
   Reads transition info from DB

   @returns -1 if db_part has not been called before
*/
static long int
db_info(struct linedb **lineinfo,
	double wav1,
	double wav2)
{
  if (!partitionread)
    return -1;

  /*
    Following is to transform stored integer of gflog and ielow into their
    real floating point values
  */
  double tablog[PANDS_NCODIDX+1];
  for(int i=1 ; i<=PANDS_NCODIDX ; i++)
    tablog[i] = pow(10,(i-16384)*0.001);
  double ratiolog=log((double)1.0+(double)1.0/(2e+6));

  struct stat st;
  if(stat(dbfilename,&st) == -1){
    mperror(MSGP_USER|MSGP_ALLOWCONT,
	    "Data file '%s' cannot be accesed by stat() in "
	    "function dbread_pands().\nThis is important to obtain "
	    "its size and hence the number of lines to be\n "
	    "examinated\n"
	    ,dbfilename);
    exit(EXIT_FAILURE);
  }
  off_t nrec = st.st_size,
    zrec = ftell(fp);
  if((zrec+nrec)/PANDS_RECLENGTH < (zrec+nrec)/(float)PANDS_RECLENGTH){
    mperror(MSGP_USER|MSGP_ALLOWCONT,
	    "Data file '%s' does not contain an integer number of "
	    "%i-bytes records!.\nAre you sure it is the right '%s' "
	    "file?\n"
	    ,dbfilename, PANDS_RECLENGTH, pands_name);
    exit(EXIT_FAILURE);
  }
  nrec/=PANDS_RECLENGTH;
  zrec/=PANDS_RECLENGTH;


  wav1 /= pands_fct/tli_fct;
  wav2 /= pands_fct/tli_fct;
  PREC_LNDATA lndb1, lndb2, 
    lnwav1 = log(wav1)/ratiolog,
    lnwav2 = log(wav2)/ratiolog;
  messagep(verbose_db, "PANDS Driver: Going to look for wavelength range %g - %g (%s)\n"
	   , wav1, wav2, pands_fct_ac);


  /* Searching for the initial point from wherever the file currently is
     until the last byte.  First, check if we are not standing at the
     required point */
  off_t irec, frec;
  sfread(&lndb1, 4, 1, fp);
  reversebytes(&lndb1, 4);

  fseek(fp, -PANDS_RECLENGTH, SEEK_END);
  sfread(&lndb2, 4, 1, fp);
  fprintf(stderr, "lnRangde1 %g\n", lndb2);
  reversebytes(&lndb2, 4);
  fprintf(stderr, "lnRangde2 %g\n", lndb2);
  messagep(verbose_db, "lnRange %g - %g (%g)\n", lndb1, lndb2, ratiolog);
  messagep(verbose_db, "Remainder database range: %g - %g (%s)\n"
	   , exp(ratiolog*lndb1), exp(ratiolog*lndb2), pands_fct_ac);

  //Make sure we have acceptable boundaries and then find the records to
  //start and finish (probing for repeated wavelengths)
  if(lnwav1 > lndb2 || lnwav2 < lndb1)
    return 0;
  if(lnwav1 <= lndb1)
    irec = zrec;
  else{
    dbreadBSf(fp, zrec, nrec, lnwav1, &irec, 8);
    long temp;
    do{
      if (irec<zrec)
	break;
      fseek(fp, PANDS_RECLENGTH*(irec--), SEEK_SET);
      sfread(&temp, 4, 1, fp);
    }while(temp >= lnwav1);
    irec++;
  }
  if(lnwav2 >= lndb2)
    frec = nrec-1;
  else{
    dbreadBSf(fp, irec, nrec, lnwav2, &frec, 8);
    long temp;
    do{
      fseek(fp, PANDS_RECLENGTH*(frec++), SEEK_SET);
      sfread(&temp, 4, 1, fp);
    }while(temp<=lnwav2);
    frec--;
  }
  messagep(verbose_db, "Target initial and final records found "
	   "in relative positions %li and %li (of range %li-%li)\n"
	   , irec-zrec, frec-zrec, zrec, nrec);

  struct linedb *line;
  messagep(4,
	   "About to initialize memory space to hold %li records.\n"
	   " I'll require %.2fMb of available memory.\n"
	   ,frec-irec+1
	   ,(frec-irec+1)*(float)sizeof(struct linedb)/1024.0/1024.0);
  if(((*lineinfo)=(struct linedb *)calloc(frec-irec+1, sizeof(struct linedb))) == NULL){
    mperror(MSGP_SYSTEM|MSGP_ALLOWCONT,
	    "Cannot allocate memory to hold all the data from %s\n"
	    "Required memory was %li bytes\n"
	    ,pands_name, (frec-irec+1)*sizeof(struct linedb));
    exit(EXIT_FAILURE);
  }
  messagep(4, "Success in memory allocation\n");



  //Read the whole thing
  fseek(fp,irec*PANDS_RECLENGTH,SEEK_SET);

  struct recordstruct{
    PREC_BREC iwl;                    //Wavelength index
    short int ielow,igflog;           //Lower energy and log(gf) indices.
  }record;
  do{

    sfread(&(record.iwl), 4, 1, fp);
    reversebytes(&(record.iwl),4);
    sfread(&(record.ielow), 2, 1, fp);
    reversebytes(&(record.ielow),2);
    sfread(&(record.igflog), 2, 1, fp);
    reversebytes(&(record.igflog),2);


    line->wl=exp(record.iwl*ratiolog)*pands_fct/tli_fct;

    //Isotopes (1h1h16o, 1h1h17o, 1h1h18o, 1h2h16o)
    if(record.ielow>0)
      line->isoid = record.igflog>0?0:1; //isotopes indices are Kurucz's - 1
    else
      line->isoid = record.igflog>0?2:3;
    line->recpos = irec + line - *lineinfo;
    line->elow   = (PREC_LNDATA)abs(record.ielow);
    line->gf     = tablog[abs(record.igflog)];

    line++;
  }while(line-*lineinfo<frec-irec+1);

  return line-*lineinfo+1;
}


/********************************************************/



static const driver_func pdriverf_pands = {
  "Partridge & Schwenke (by Kurucz) driver",
  &db_find,
  &db_open,
  &db_close,
  &db_info,
  &db_part
};

driver_func *
initdb_pands()
{
  return (driver_func *)&pdriverf_pands;
}
