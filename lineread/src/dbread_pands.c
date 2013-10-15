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

/* List of functions defined here:

static void sfread_fcn(void *pointer, size_t size, size_t nmemb, FILE *f,
                       long int line)
    Safe(?) read.

static inline void dbreadBSf(FILE *fpbs, PREC_NREC initial, PREC_NREC final,
                             double lookfor, PREC_NREC *resultp, int reclength)
    Binary search in file between initial and final records looking
    for the record which first value is lookfor.  Store found record's
    index in resultp.

static _Bool db_find(const char *name)
    Determines if a file can be read with this data-base reader.  If filename
    equals the default filename return 1, otherwise return 0.

static int db_open(char *dbname, char *dbaux)
    Opens data base and auxiliary file, set FILE pointers.

static int db_close()
    Close data bases files pointers and free memory.

static _Bool db_part(char **name, unsigned short *nT, PREC_TEMP **T,
                     unsigned short *niso, char ***isonames, PREC_MASS **mass,
                     PREC_Z ***Z, PREC_CS ***CS)
    Reads partition info, allocates and stores values in input arrays, 
    closes the auxiliary file pointer.

static long int db_info(struct linedb **lineinfo, double wav1, double wav2)
    Check the data base has the right format and size.  Set the
    indices of the initial and final line records to read.  Allocate
    linedb.  Read the data base and store values. Return the number of
    records read.

driver_func *initdb_pands()
    Initialize the data-base reader struct.
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
#define reversebytes(pointer, nbytes) do{}while(0)  /* Do nothing  */
#endif   /* __BYTE_ORDER */

#define NUM_ISOT 4                  /* Number of isotopes          */
#define RWATER (3.2e-8/2.0)         /* Water molecule radius (cm)  */
#define PI (3.141592653589793)      /* Pi                          */
#define SIGWATER (PI*RWATER*RWATER) /* Water cross section (cm^2)  */
#define PREC_BREC int               /* Type for the Big record (?) */

#define PANDS_RECLENGTH 8
#define PANDS_NCODIDX   32786

static const char *defaux = "h2opartfn.dat";
static const char *deftarget = "h2ofast.bin";
static const char *pands_name = "Partridge & Schwenke (1997)";
static const char *pands_iso[]      = {"1H1H16O",   "1H1H17O",   "1H1H18O",   "1H2H16O"};
static const PREC_MASS pands_mass[] = {18.01056468, 19.01478156, 20.01481046, 19.01684143};

/* Wavelength comes stored in nanometers */
static double      pands_fct    = 1e-7;
static const char *pands_fct_ac = "nm";

static char *dbfilename = NULL, *dbauxname = NULL;
static FILE *fp = NULL, *fpaux = NULL; /* File pointers to read */
static _Bool partitionread = 0;
static int verbose_db = 10;            /* FINDME: what am I? */


#define sfread(pointer, size, nmemb, f) sfread_fcn(pointer, size, nmemb, f, __LINE__)
/* FIDME: Safe(?) read  */
static void
sfread_fcn(void *pointer,   /* Pointer to memory where to store the data */
           size_t size,     /* Size of each element to be read           */
           size_t nmemb,    /* Number of element to be read.             */
           FILE *f,         /* Pointer to file to read from.             */
           long int line){  /* FINDME: what is this?                     */
  size_t read = fread(pointer, size, nmemb, f);
  /* Return an error if the read did not succeed:                        */
  if (read != nmemb)
    fprintf(stderr,
            "%s (%li): Number of elements read was only %li of %li\n",
            __FILE__, line, (long)read, (long)nmemb);
}

/* \fcnfh
   Free memory used by this driver  */
static void
pands_free(){
  if (dbauxname)
    free(dbauxname);
}


/* \fcnfh
  Binary search in file pointed by 'fpbs'(FILE *) between 'initial'(PREC_NREC)
   and 'final'(PREC_NREC) looking for the record which first value is 
  'lookfor'(PREC_NREC).  Store its the record index in 'resultp'(PREC_NREC *).
   Records are of length 'reclength'(int).                                    */
static inline void
dbreadBSf(FILE *fpbs,           /* File pointer           */
          PREC_NREC initial,    /* Initial index          */
          PREC_NREC final,      /* Last index             */
          double lookfor,       /* Target value           */
          PREC_NREC *resultp,   /* Index of target value  */
          int reclength){       /* Length of a record     */
  long int irec1, irec2;
  PREC_BREC temp;
  const int trglength = sizeof(PREC_BREC);

  irec1 = initial;
  irec2 = final;
  do{
    *(resultp) = (irec2+irec1)/2;                /* index of middle record */
    fseek(fpbs, reclength*(*resultp), SEEK_SET); /* Go there               */
    sfread(&temp, trglength, 1, fpbs);           /* Read the first item    */
    reversebytes(&temp, trglength);              /* Reverse if needed      */
    if(lookfor>temp)
      irec1 = *(resultp);
    else
      irec2 = *(resultp);
  }while (irec2-irec1>1);
  *resultp = irec1;
}


/* \fncfh
   Determines if a file can be read with this data-base reader.  If filename
    equals the default filename return 1, otherwise return 0.                 */
static _Bool
db_find(const char *name){
  /* Compares name (the filename, actually) to the default Partridge and
     Schwenke filename.  If they match, assumes it's a pands type database */
  int len = strlen(name), lent = strlen(deftarget);
  
  if (len >= lent &&
      strncmp(deftarget, name+len-lent, lent) == 0)
    return 1;

  return 0;
}


/* \fcnfh
   Opens data base and auxiliary file, set FILE pointers */
static int
db_open(char *dbname,  /* Data base filename           */
        char *dbaux){  /* Auxiliary data base filename */
  /* Open database, use fp to handle the file         */
  if((fp = fopen(dbname,"r")) == NULL)
    mperror(MSGP_USER, "Could not open file '%s' for reading\n", dbname);
  dbfilename = dbname;

  /* Set dbauxname, the auxiliar database filename: */
  if (dbaux)
    dbauxname = dbaux;
  else{
    int len  = strlen(dbname),    /* Input database filename length          */
        lent = strlen(deftarget), /* Default filename (without path) length  */
        lena = strlen(defaux);    /* Default partition-function filename len */
    dbauxname = (char *)calloc(len-lent+lena+2, sizeof(char));
    strncpy(dbauxname, dbname, len-lent); /* Copy path                       */
    strcpy(dbauxname+len-lent, defaux);   /* Copy part. function filename    */
  }

  /* Open auxiliary data base if there is one, use fpaux to handle the file  */
  if((fpaux = fopen(dbauxname,"r")) == NULL)
    mperror(MSGP_USER,
            "Could not open auxiliary file '%s' for reading\n",
            dbauxname);
  return LR_OK;
}


/* \fcnfh
   Closes data-base files pointers and frees memory */
static int
db_close(){
  if(fp)
    fclose(fp);
  if(fpaux)
    fclose(fpaux);
  pands_free();
  return LR_OK;
}


/* \fcnfh
   Reads partition info, allocates and stores values in input arrays, 
   closes the auxiliary file pointer                                   */
static _Bool
db_part(char **name,           /* FINDME: data-base name?       */
        unsigned short *nT,    /* Number of temperature samples */
        PREC_TEMP **T,         /* Temperatures array            */
        unsigned short *niso,  /* Number of isotopes            */
        char ***isonames,      /* Isotopes names                */
        PREC_MASS **mass,      /* Isotpes masses                */
        PREC_Z ***Z,           /* FIMDE: what am I?             */
        PREC_CS ***CS){        /* Cross sections                */

  long maxline = 300;            /* Max characters to read per line          */
  char line[maxline], *lp, *lp2; /* The lines being readed                   */
  int ignorelines=5;             /* Number of ignored lines at file begining */
  unsigned short nIso = *niso = NUM_ISOT; /* Default number of isotopes      */
  int nt=0,                      /* Read number of temperatures              */
      nta=8;                     /* Allocated number of temperatures         */

  /* Read through the header */
  for(int i=0; i<ignorelines; i++)
    fgets(line, maxline, fpaux);

  /* Allocate output arrays: */
  PREC_Z  **z     = *Z  = (PREC_Z  **)calloc(nIso, sizeof(PREC_Z *));
  PREC_CS **cs    = *CS = (PREC_CS **)calloc(nIso, sizeof(PREC_CS *));
  PREC_TEMP *temp =      (PREC_TEMP *)calloc(nta,  sizeof(PREC_TEMP));
  *z  = (PREC_Z    *)calloc(nIso*nta, sizeof(PREC_Z));
  *cs = (PREC_CS   *)calloc(nIso*nta, sizeof(PREC_CS));
  for(unsigned short i=1; i<nIso; i++){
    z[i]  = *z  + i*nta; /* Initialize the pointers elements (also pointers) */
    cs[i] = *cs + i*nta;
  }

  /* Read file lines while there is a line to read:                      */
  while((lp=fgets(line, maxline, fpaux)) != NULL){
    /* If read lines reach nta, reallocate to double the size of arrays */
    if(nt == nta){
      unsigned int old = nta;
      *z   = (PREC_Z    *)realloc(*z,   nIso*(nta<<=1)*sizeof(PREC_Z));
      *cs  = (PREC_CS   *)realloc(*cs,  nIso*nta      *sizeof(PREC_CS));
      temp = (PREC_TEMP *)realloc(temp, nta           *sizeof(PREC_TEMP));
      for(unsigned short i=nIso-1; i>0; i--){
        z[i]  = *z  + i*nta;
        cs[i] = *cs + i*nta;
        /* FINDME: From where to where am I moving the data? */
        memmove(z[i],  *z  + i*old, old*sizeof(PREC_Z));
        memmove(cs[i], *cs + i*old, old*sizeof(PREC_CS));
      }
    }

    /* Get temperature */
    temp[nt] = strtod(lp, &lp2);

    for(int i=0; i<nIso; i++){
      /* If strtod read an empty or blank string lp2 and lp points to the
         same character (end of line reached before reading all isotopes)  */
      if(lp2==lp){
        mperror(MSGP_USER|MSGP_ALLOWCONT,  /* FINDME: typo? read_Zpands ?? */
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
      /* Fill in Z and cross section: */
      lp = lp2;
      z[i][nt]  = strtod(lp, &lp2); /* Reads lp, advances lp2 */
      cs[i][nt] = SIGWATER;
    }
    nt++;
  }

  /* Reallocate to discard remaining unused space:              */
  *nT=nt; /* Set output number of temperatures read             */
  MESSAGEP(verbose_db, "P&S driver: found %i temperatures\n", nt);
  *T = (PREC_TEMP *)realloc(temp, nt*sizeof(PREC_TEMP));
  for(unsigned short i=0; i<nIso; i++){
    (*Z)[i]  = *z  + i*nt;
    (*CS)[i] = *cs + i*nt;
    memmove((*Z)[i],  *z  + i*nta, nt*sizeof(PREC_Z));
    memmove((*CS)[i], *cs + i*nta, nt*sizeof(PREC_CS));
  }
  **Z  = (PREC_Z  *)realloc(**Z,  nIso*nt*sizeof(PREC_Z));
  **CS = (PREC_CS *)realloc(**CS, nIso*nt*sizeof(PREC_CS));

  /* Close partition function file:                             */
  if(fpaux)
    fclose(fpaux);
  fpaux=NULL;

  /* FINDME: Fill in WHAT name?, isotopes names, and isotopes masses */
  *name     = strdup(pands_name);
  *isonames = (char     **)calloc(nIso, sizeof(char *));
  *mass     = (PREC_MASS *)calloc(nIso, sizeof(PREC_MASS));
  for (int i=0; i<nIso; i++){
    (*mass)[i] = pands_mass[i];
    (*isonames)[i] = strdup(pands_iso[i]);
  }

  partitionread = 1;
  return 0;
}


/* \fcnfh
   Reads transition info from the data bases */
static long int
db_info(struct linedb **lineinfo, /* Output linedb with lines info */
        double wav1,              /* Low-wavelength limit to read  */
        double wav2){             /* High-wavelength limit to read */
  /* Returns -1 if db_part has not been called before   */
  if (!partitionread)
    return -1;

  /* Factors to transform stored-integer values of gflog and elow
     into their real floating point values                         */
  double tablog[PANDS_NCODIDX+1];
  for(int i=1; i<=PANDS_NCODIDX; i++)
    tablog[i] = pow(10,(i-16384)*0.001);
  double ratiolog = log((double)1.0 + (double)1.0/(2e+6));

  struct stat st;  /* FINDME: Where is the struct stat defined? */
  /* Check if data can be accesed by stat: */
  if(stat(dbfilename, &st) == -1){
    mperror(MSGP_USER|MSGP_ALLOWCONT,
            "Data file '%s' cannot be accesed by stat() in "
            "function dbread_pands().\nThis is important to obtain "
            "its size and hence the number of lines to be\n examinated.\n",
            dbfilename);
    exit(EXIT_FAILURE);
  }
  off_t nrec = st.st_size, /* Byte position of the end of the data */
        zrec = ftell(fp);  /* Current byte position of the data    */
  /* Check that the record lenght is a factor of PANDS_RECLENGTH:  */
  if((zrec+nrec)/PANDS_RECLENGTH < (zrec+nrec)/(float)PANDS_RECLENGTH){
    mperror(MSGP_USER|MSGP_ALLOWCONT,
            "Data file '%s' does not contain an integer number of "
            "%i-bytes records!.\nAre you sure it is the right '%s' file?\n",
            dbfilename, PANDS_RECLENGTH, pands_name);
    exit(EXIT_FAILURE);
  }
  /* Divide by record-length to get indices for records: */
  nrec/=PANDS_RECLENGTH;
  zrec/=PANDS_RECLENGTH;

  /* Convert from TLI to pands wavelength units: */
  wav1 /= pands_fct/tli_fct;
  wav2 /= pands_fct/tli_fct;
  PREC_BREC lndb1, lndb2;             /* Big records (?)      */
  /* FINDME: this is a terrible name, what does lndb stands for? */
  double lnwav1 = log(wav1)/ratiolog, /* Wavelengths as given */
         lnwav2 = log(wav2)/ratiolog; /* in the data base     */
  MESSAGEP(verbose_db, "P&S Driver: Going to look for wavelength range "
                       "%g - %g (%s)\n", wav1, wav2, pands_fct_ac);

  /* Search between the current and the last byte position in the file.
     First, check if we are not standing at the required point.          */
  off_t irec, frec;         /* Initial and final record indices          */
  sfread(&lndb1, 4, 1, fp); /* Read current record's wavelength          */
  reversebytes(&lndb1, 4);  /* Reverse byte order if 
                               __BYTE_ORDER == __BIG_ENDIAN              */
  fseek(fp, -PANDS_RECLENGTH, SEEK_END); /* Go to last record in file    */
  sfread(&lndb2, 4, 1, fp);              /* Read last recod's wavelength */
  reversebytes(&lndb2, 4);
  MESSAGEP(verbose_db, "P&S driver: ln range %i - %i (%g)\n", 
           lndb1,               lndb2,               ratiolog);
  MESSAGEP(verbose_db, "P&S driver: remainder database range: %g - %g (%s)\n",
           exp(ratiolog*lndb1), exp(ratiolog*lndb2), pands_fct_ac);

  /* Set the record indices boundaries (probing for repeated wavelengths): */
  /* The range I want and the range I have don't overlap:                  */
  if(lnwav1 > lndb2 || lnwav2 < lndb1)
    return 0;
  /* The lowest data wavelength is greater than the starting wavelength
     I want:                                                               */
  if(lnwav1 <= lndb1)
    irec = zrec;
  /* Binary search in [zrec, nrec] to find record index for lnwav1 (irec)  */
  else{
    /* FINDME: Edited code, changed hard-coded 8 to PANDS_RECLENGTH        */
    dbreadBSf(fp, zrec, nrec, lnwav1, &irec, PANDS_RECLENGTH);
    MESSAGEP(verbose_db, "P&S driver: Found initial wavelength (%f) at "
                         "record %li, checking twins...", lnwav1, irec);
    int temp;
    /* FINDME: Why doing this DO_WHILE? irec can't be smaller than zrec */
    do{
      if (irec<zrec)
        break;
      fseek(fp, PANDS_RECLENGTH*(irec--), SEEK_SET); /* Move the file pointer to irec */
      sfread(&temp, 4, 1, fp);
    }while(temp >= lnwav1);
    irec++;
    MESSAGEP(verbose_db, "done (%li)\n", irec);
  }
  /* The largest data wavelength is shorter than the last wavelength
     I want:                                                              */
  if(lnwav2 >= lndb2)
    frec = nrec-1;
  /* Binary search to find record index for lnwav2 (frec):                */
  else{
    dbreadBSf(fp, irec, nrec, lnwav2, &frec, 8);
    MESSAGEP(verbose_db, "P&S driver: Found final wavelength (%f) at "
                         "record %li, checking twins...", lnwav2, frec);
    int temp;
    do{
      fseek(fp, PANDS_RECLENGTH*(frec++), SEEK_SET);
      sfread(&temp, 4, 1, fp);
    }while(temp <= lnwav2);
    frec--;
    MESSAGEP(verbose_db, "done (%li)\n", frec);
  }
  MESSAGEP(verbose_db, "P&S driver: Target initial and final records found in "
                       "relative positions %li and %li (of range %li-%li)\n",
                       irec-zrec, frec-zrec, zrec, nrec);

  /* Attempt to allocate memory for lines: */
  struct linedb *line;
  messagep(5, "\nP&S driver: About to initialize memory space to hold %li "
              "records.\nP&S driver: I'll require %.2fMb of memory.\n",
           frec-irec+1,
           (frec-irec+1)*(float)sizeof(struct linedb)/1024.0/1024.0);
  if(((*lineinfo) = (struct linedb *)calloc(frec-irec+1, 
                                            sizeof(struct linedb))) == NULL){
    mperror(MSGP_SYSTEM|MSGP_ALLOWCONT,
            "P&S driver: Cannot allocate memory to hold all the data from %s\n"
            "P&S driver: Required memory was %li bytes\n",
            pands_name, (frec-irec+1)*sizeof(struct linedb));
    exit(EXIT_FAILURE);
  }
  /* Declare the pointer to *lineinfo: */
  line = *lineinfo;
  messagep(5, "P&S driver: Success in memory allocation\n");

  /* Read the whole thing */
  fseek(fp, irec*PANDS_RECLENGTH, SEEK_SET);

  /* FINDME: Why is this defined here and not outside? */
  struct recordstruct{
    /* FINDME: These are not indices, these are values! */
    PREC_BREC iwl;            /* Wavelength index                  */
    short int ielow, igflog;  /* Lower energy and log(gf) indices  */
  }record;

  do{
    /* Read a record:  */
    sfread(&(record.iwl),    4, 1, fp);
    reversebytes(&(record.iwl),    4);
    sfread(&(record.ielow),  2, 1, fp);
    reversebytes(&(record.ielow),  2);
    sfread(&(record.igflog), 2, 1, fp);
    reversebytes(&(record.igflog), 2);

    /* Store info (convert to correct values first):                        */
    line->wl = exp(record.iwl*ratiolog) * pands_fct/tli_fct;
    /* Isotopes (1h1h16o, 1h1h17o, 1h1h18o, 1h2h16o)                        */
    if(record.ielow>0)
      line->isoid = record.igflog>0 ? 0:1; /* Isotopes ID's are Kurucz's-1  */
    else
      line->isoid = record.igflog>0 ? 2:3;
    line->recpos = irec + line - *lineinfo;
    line->elow   = (PREC_LNDATA)abs(record.ielow);
    line->gf     = tablog[abs(record.igflog)]; /* convert to correct values */

    line++; /* Move record pointer to next position                         */
  }while(line-*lineinfo<frec-irec+1);

  /* Return number of records read: */
  return line-*lineinfo; 
}


/* Array of functions to be parsed as arguments in initdb_pands: */
/* FINDME: EXPLAIN MY NAME, WHAT DOES 'PDRIVERF' STANDS FOR? */
static const driver_func pdriverf_pands = {
  "Partridge & Schwenke (by Kurucz) driver",
  &db_find,
  &db_open,
  &db_close,
  &db_info,
  &db_part
};

/* Initialize the data base reader structure */
driver_func *
initdb_pands(){
  return (driver_func *)&pdriverf_pands;
}
