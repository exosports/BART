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
#define PREC_BREC int		/* Type for the Big record  */

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
  size_t read = fread(pointer, size, nmemb, f);
  if (read != nmemb)
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
  long int irec1,irec2;
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
	memmove(z[i],  *z  + i*old, old*sizeof(PREC_Z));
	memmove(cs[i], *cs + i*old, old*sizeof(PREC_CS));
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
    nt++;
  }

  *nT=nt;
  MESSAGEP(verbose_db, "P&S driver: found %i temperatures\n", nt);
  *T = (PREC_TEMP *)realloc(temp, nt*sizeof(PREC_TEMP));

  for(unsigned short i=0 ; i<nIso ; i++){
    (*Z)[i]  = *z  + i*nt;
    (*CS)[i] = *cs + i*nt;
    memmove((*Z)[i],  *z  + i*nta, nt*sizeof(PREC_Z));
    memmove((*CS)[i], *cs + i*nta, nt*sizeof(PREC_CS));
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
  PREC_BREC lndb1, lndb2;
  double lnwav1 = log(wav1)/ratiolog,
    lnwav2 = log(wav2)/ratiolog;
  MESSAGEP(verbose_db, "P&S Driver: Going to look for wavelength range %g - %g (%s)\n"
	   , wav1, wav2, pands_fct_ac);


  /* Searching for the initial point from wherever the file currently is
     until the last byte.  First, check if we are not standing at the
     required point */
  off_t irec, frec;
  sfread(&lndb1, 4, 1, fp);
  reversebytes(&lndb1, 4);

  fseek(fp, -PANDS_RECLENGTH, SEEK_END);
  sfread(&lndb2, 4, 1, fp);
  reversebytes(&lndb2, 4);
  MESSAGEP(verbose_db, "P&S driver: lnRange %i - %i (%g)\n", lndb1, lndb2, ratiolog);
  MESSAGEP(verbose_db, "P&S driver: remainder database range: %g - %g (%s)\n"
	   , exp(ratiolog*lndb1), exp(ratiolog*lndb2), pands_fct_ac);

  //Make sure we have acceptable boundaries and then find the records to
  //start and finish (probing for repeated wavelengths)
  if(lnwav1 > lndb2 || lnwav2 < lndb1)
    return 0;
  if(lnwav1 <= lndb1)
    irec = zrec;
  else{
    dbreadBSf(fp, zrec, nrec, lnwav1, &irec, 8);
    MESSAGEP(verbose_db,
	     "P&S driver: Found initial wavelength (%f) at record %li, checking twins..."
	     , lnwav1, irec);
    int temp;
    do{
      if (irec<zrec)
	break;
      fseek(fp, PANDS_RECLENGTH*(irec--), SEEK_SET);
      sfread(&temp, 4, 1, fp);
    }while(temp >= lnwav1);
    irec++;
    MESSAGEP(verbose_db, "done (%li)\n", irec);
  }
  if(lnwav2 >= lndb2)
    frec = nrec-1;
  else{
    dbreadBSf(fp, irec, nrec, lnwav2, &frec, 8);
    MESSAGEP(verbose_db,
	     "P&S driver: Found final wavelength (%f) at record %li, checking twins..."
	     , lnwav2, frec);
    int temp;
    do{
      fseek(fp, PANDS_RECLENGTH*(frec++), SEEK_SET);
      sfread(&temp, 4, 1, fp);
    }while(temp <= lnwav2);
    frec--;
    MESSAGEP(verbose_db, "done (%li)\n", frec);
  }
  MESSAGEP(verbose_db, "P&S driver: Target initial and final records found "
	   "in relative positions %li and %li (of range %li-%li)\n"
	   , irec-zrec, frec-zrec, zrec, nrec);

  struct linedb *line;
  messagep(5,
	   "\nP&S driver: About to initialize memory space to hold %li records.\n"
	   "P&S driver: I'll require %.2fMb of available memory.\n"
	   ,frec-irec+1
	   ,(frec-irec+1)*(float)sizeof(struct linedb)/1024.0/1024.0);
  if(((*lineinfo)=(struct linedb *)calloc(frec-irec+1, sizeof(struct linedb))) == NULL){
    mperror(MSGP_SYSTEM|MSGP_ALLOWCONT,
	    "P&S driver: Cannot allocate memory to hold all the data from %s\n"
	    "P&S driver: Required memory was %li bytes\n"
	    ,pands_name, (frec-irec+1)*sizeof(struct linedb));
    exit(EXIT_FAILURE);
  }
  line = *lineinfo;
  messagep(5, "P&S driver: Success in memory allocation\n");



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


    line->wl = exp(record.iwl*ratiolog) * pands_fct/tli_fct;

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

  return line-*lineinfo;
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
