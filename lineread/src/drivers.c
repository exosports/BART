/*********************************************************************
 * 
 * drivers.c - File that finds and executes the appropriate driver for
 * each database. Part of lineread program.
 *
 * Copyright (C) 2006 Patricio Rojo
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 *
 **********************************************************************/

/* List of functions defined here: 

static inline void pfwrite(void *pointer, size_t size, size_t nmemb, FILE *fp)
    Print to screen  and write to file the specified data

static inline int find_dbd(char *file, unsigned short nfcn)
    Find the appropriate driver to read 'file', returns its index
    or -1 if none is found.

int db_drivers(struct hints *hint)
    Main driver routine to read data bases and store required info into a
    TLI file.  

void drivers_free(struct hints *hint)
    Close drivers, free pointers 

int find_alldrivers(struct hints *hint, unsigned short nfcn)
    Initialize drivers and find which one reads each data base in hint

int setdriversnoutput(struct hints *hint)
    Opens the database and auxiliar partition function files 
    for reading.  Opens the output TLI file for writing.  
    Prints basic info to the output file.

int readwritepartition(unsigned short *acum)
    Read partition function file.  Store info into output file.
    Count number of isotopes per data base. 

int readwritetransition(unsigned short *acum, double ini, double fin,
                        double del)
    Read transitions from data bases in chunks, sort them, and write
    them to the output file.
*/


#include <lineread.h>
#define _VERSION_ONLY_DECLARATIONS_
#include <version_lr.h>

/* Add here a content struct for each driver, which should also have an
   "extern" entry in lineread.h as well as its proto_ file included. */
/* Declare and initializes an array of the different driver_func     */
/* The internal structure is defined in include/structures_lr.h      */
/* Each driver type is defined in src/dbread_type.c                  */
static driver_func * (*init_driver[])() = {
  &initdb_debug,
  //  &initdb_hitran4,
  &initdb_pands,
  &initdb_text
};

static driver_func **drivers=NULL; /* Pointer to the drivers in init_driver */
static char *outfilename = NULL;   /* Output TLI filename                   */
static unsigned short ndb;         /* Number of data bases in hint          */
static _Bool dry;                  /* Is it a dry run?                      */
static unsigned short *DBdriver = NULL; /* driver index to reads each db    */
static FILE *fpout = NULL;         /* stream to write to the outputfile     */
static int verbose_TLIout = 10;    /* Comment me */
static int verbose_TLIout2 = 5;    /* Comment me */


/* Print to screen and write to file the specified data.                 */
static inline void
pfwrite(void *pointer, /* Pointer to the array of elements to be written */
        size_t size,   /* Size of each element to be writen              */
        size_t nmemb,  /* Number of elements to be writen                */
        FILE *fp){     /* Pointer to file to write in                    */
  if(size*nmemb)
    MESSAGEP(verbose_TLIout, "%04hx%s\n",
             *((unsigned short *)pointer), 
             size==sizeof(unsigned short)?"":" ...");
  if(!dry)
    fwrite(pointer, size, nmemb, fp);
}


/* \fcnfh
    Find the appropriate driver to read 'file', returns its
    index or -1 if none is found.                               */
static inline int
find_dbd(char *file,            /* Data base filename           */
         unsigned short nfcn){  /* Total number of driver types */
  for (int i=0 ; i<nfcn ; i++)
    if((*(drivers[i]->find))(file))
      return i;
  return -1;
}


/* \fcnfh
   Main driver routine to read data bases and store required info
   into a TLI file.                                                    */
int
db_drivers(struct hints *hint){
  /* Number of driver_func types                                       */
  unsigned short nfcn = sizeof(init_driver)/sizeof(driver_func *);
  ndb = hint->ndb;  /* Number of data bases in hint                    */
  dry = hint->dry;  /* is it a dry run?                                */
  /* Driver index to read each data base:                              */
  DBdriver = (unsigned short *)calloc(ndb, sizeof(unsigned short));

  /* Initialize the drivers and find the right one for each data base: */
  find_alldrivers(hint, nfcn);

  /* Set up drivers & output file:                                     */
  setdriversnoutput(hint);
  
  /* Cummulative number of isotopes per database:                      */
  unsigned short acum[ndb+1]; 
  /* Read & write the partition info:                                  */
  readwritepartition(acum);

  /* Read & write the line transitions:                                */
  readwritetransition(acum, hint->iniw, hint->finw, hint->delw);

  return LR_OK;
}


/* \fcnfh
   Close drivers, free pointers              */
void  
drivers_free(struct hints *hint) {
  /* Free the drivers:                       */
  for (int i=0; i<hint->ndb; i++)
    (*(drivers[DBdriver[i]]->close))();
  /* Free the drivers' indices:              */
  if(DBdriver)
    free(DBdriver);
  /* Free the TLI file pointer               */
  if (fpout != NULL && fpout != stdout && !dry)
    fclose(fpout);
  /* Free the pointer to drivers             */
  if(drivers)
    free(drivers);
  /* Free the pointer to the output filename */
  if(outfilename)
    free(outfilename);
}


/* \fcnfh
   Initialize drivers and find which one reads each data base in hint  */
int
find_alldrivers(struct hints *hint,
                unsigned short nfcn){
  /* Initialize drivers, get their names:                              */
  int allocnm = 8, /* Size of array containing the drivers names       */
     lennm    = 0; /* Length of all drivers names + 2 chars x each one */
  /* Array concatenating all the names of the dirvers:                 */
  char *alldb = (char *)calloc(allocnm, sizeof(char));
  /* Initialize pointer to drivers:                                    */
  drivers = (driver_func **)calloc(nfcn, sizeof(driver_func *));
  for (int i=0; i<nfcn; i++){
    /* Initialize the drivers:                                         */
    drivers[i] = (*(init_driver[i]))();
    /* Extract drivers' name.  Add it's length's to lennm:             */
    char *name = (char *)drivers[i]->name;
    lennm += strlen(name) + 2;
    /* Double the value of allocnm until it's larger than lennm:       */
    while (lennm >= allocnm)
      alldb = (char *)realloc(alldb, (allocnm<<=1)*sizeof(char));
    strcat(alldb, name);
    if(i<nfcn-1) strcat(alldb, ", ");
  }
  /* Realloc alldb to have size lennm.  FINDME: Why to make it so
     complicated then??                                                */
  alldb = (char *)realloc(alldb, lennm*sizeof(char));

  /* Find the appropriate driver for each data base:           */
  messagep(2, "Finding drivers for %i database(s)", ndb);
  messagep(3, ":\n");
  int result; /* Index of driver that reads a data base        */
  for (int i=0; i<ndb; i++){
    if((result = find_dbd(hint->db[i], nfcn)) < 0){
      /* Throw an error if no driver works for a data base     */
      mperror(MSGP_USER|MSGP_ALLOWCONT,
              "The file '%s' could not be associated to any "
              "supported database.  Currently, lineread can read: %s.\n",
               hint->db[i], alldb);
      lineread_free();
      free(alldb);
      exit(EXIT_FAILURE);
    }
    messagep(3, " %s: found %s", hint->db[i], drivers[result]->name);
    messagep(2, ".");
    messagep(3, "\n");
    /* Store the index in DBdriver and hint.dbd */
    DBdriver[i] = hint->dbd[i] = result;
  }
  free(alldb);
  if(verblevel <3) messagep(2, " ");
  messagep(2, "done\n");

  return LR_OK;
}


/* \fcnfh
   Opens the database and auxiliar partition function files 
   for reading.  Opens the output TLI file for writing.  
   Prints basic info to the output file.                         */
int
setdriversnoutput(struct hints *hint){
  /* Get name of output file:                                    */
  /* FINDME: WHY COMPARE UP TO 2 CHARACTERS? */
  outfilename = strdup(strncmp("-", hint->datafile, 2)==0 ?
                       "Standard Output":hint->datafile);

  /* Open the database and partition function files for reading: */
  for (int i=0; i<ndb; i++)
    (*(drivers[DBdriver[i]]->open))(hint->db[i], hint->dbaux[i]);

  /* Open hint.datafile for writing (if not a dry run):          */
  if (!dry){
    if(strcmp("-", hint->datafile) == 0)
      fpout = stdout;
    else if((fpout=fopen(hint->datafile, "w"))==NULL){
      mperror(MSGP_USER|MSGP_ALLOWCONT,
              "Data file '%s' cannot be opened for writing.\n",
              hint->datafile);
      lineread_free();
      exit(EXIT_FAILURE);
    }
  }
  else
    messagep(2, "Dry-");
  messagep(2, "Opened input DB%s (%i), output stream (%s), and writing "
              "TLI header:\n", ndb>0?"s":"", ndb, outfilename); 

  /* Magic number, which in big endian would be:
     abb3b6ff = {(char)0xff-'T',(char)0xff-'L',(char)0xff-'I',(char)0xff},
     or in little endian: ffb6b3ab                                         */
  int32_t magic=((0xff-'T')<<24)|((0xff-'L')<<16)|((0xff-'I')<<8)|(0xff);
  char *undefined_string = "";
  unsigned short rn = strlen(undefined_string);

  /* Print basic info:  */
  MESSAGEP(verbose_TLIout2, "Magic number:  %i\n", magic); 
  pfwrite(&magic,           sizeof(int32_t),         1, fpout);
  MESSAGEP(verbose_TLIout2, "TLIVersion:    %u\n", TLIversion); 
  pfwrite(&TLIversion,      sizeof(unsigned short),  1, fpout); 
  MESSAGEP(verbose_TLIout2, "Version:       %u\n", version); 
  pfwrite(&version,         sizeof(unsigned short),  1, fpout); 
  MESSAGEP(verbose_TLIout2, "Revision:      %u\n", revision); 
  pfwrite(&revision,        sizeof(unsigned short),  1, fpout);
  MESSAGEP(verbose_TLIout2, "Initial wl: %7g\n", hint->iniw);
  pfwrite(&(hint->iniw),    sizeof(double),          1, fpout);
  MESSAGEP(verbose_TLIout2, "Final   wl: %7g\n", hint->finw);
  pfwrite(&(hint->finw),    sizeof(double),          1, fpout);
  MESSAGEP(verbose_TLIout2, "Length undefined: %u\n", rn); 
  pfwrite(&rn,              sizeof(unsigned short),  1, fpout); 
  MESSAGEP(verbose_TLIout2, "Undefined:        %s", rn==0?"\n":""); 
  pfwrite(undefined_string, sizeof(char),           rn, fpout);

  messagep(2, "done\n");

  return 1;
}


/* \fcnfh
   Read partition function file.  Store info into output file.
   Count number of isotopes per data base.                      */
int
readwritepartition(unsigned short *acum){
  /* Set to 0 all values of acum */
  memset(acum, 0, (ndb+1)*sizeof(unsigned short));
  /* rdb is the "real DB" counter, while i (below) is the counter 
     of the number of files that are read (FINDME: WHAT'S THE DIFFERENCE?): */
  int rdb=0;

  /* Store the number of data bases:                               */
  messagep(2, "Reading and writing partition function information:\n");
  long ndbpos = ftell(fpout); /* Get current file pointer position */
  MESSAGEP(verbose_TLIout, "Number of DBs:    ");
  pfwrite(&ndb,  sizeof(unsigned short), 1, fpout);

  /* Read info:                                                       */
  for (unsigned short i=0; i<ndb; i++, rdb++){
    char *name, **isonames;  /* FINDME: WHAT name? and isotopes names */
    unsigned short nT, niso; /* Number of temperatures and isotopes   */
    /* Declare variables for partition function info:                 */
    PREC_MASS *mass; /* Isotopes masses                               */
    PREC_TEMP *T;    /* Temperatures sampled                          */
    PREC_Z    **Z;   /* ??                                            */
    PREC_CS   **CS;  /* Cross section                                 */
    /* Read the partition function file:                              */
    _Bool moredb = (*(drivers[DBdriver[i]]->part))
                   (&name, &nT, &T, &niso, &isonames, &mass, &Z, &CS);

    /* Store partition function info into output file: */
    MESSAGEP(verbose_TLIout2," For DB #%i (file #%i):\n", rdb, i);
    unsigned short rn = strlen(name);
    MESSAGEP(verbose_TLIout2," Length of name:  %u\n", rn);
    pfwrite(&rn,    sizeof(unsigned short),  1, fpout); 
    MESSAGEP(verbose_TLIout2," Name:  %s\n", name);
    pfwrite(name,   sizeof(char),           rn, fpout);  /* FINDME: Should  be &name ?? */
    MESSAGEP(verbose_TLIout2," Number of temperatures: %u\n", nT);
    pfwrite(&nT,    sizeof(unsigned short),  1, fpout);
    MESSAGEP(verbose_TLIout2," Number of isotopes:     %u\n", niso);
    pfwrite(&niso, sizeof(unsigned short),  1, fpout);
    MESSAGEP(verbose_TLIout," Temperatures:           ");
    pfwrite(T,      sizeof(PREC_TEMP),      nT, fpout); /* FINDME: Should be &T ?? */
    MESSAGEP(verbose_TLIout+1, "Each T: >>");
    for(int k=0; k<nT ; k++)
      MESSAGEP(verbose_TLIout+1, "%g__",T[k]);
    MESSAGEP(verbose_TLIout+1, "<<\n");
    
    /* Store each isotopes info: */
    for(int j=0; j<niso; j++){
      MESSAGEP(verbose_TLIout,"  For isotope %i:\n", j);
      rn = strlen(isonames[j]);
      MESSAGEP(verbose_TLIout,"   Length of name: ");
      pfwrite(&rn,         sizeof(unsigned short),  1, fpout);
      MESSAGEP(verbose_TLIout,"   Name:           ");
      pfwrite(isonames[j], sizeof(char),           rn, fpout); /* FINDME: do these need an & ?? */
      MESSAGEP(verbose_TLIout,"   Masses:         ");
      pfwrite(mass+j,      sizeof(PREC_MASS),       1, fpout);
      MESSAGEP(verbose_TLIout,"   Partition:      ");
      pfwrite(Z[j],        sizeof(PREC_Z),         nT, fpout);
      MESSAGEP(verbose_TLIout+1, "Each Z: >>");
      for(int k=0; k<nT; k++)
        MESSAGEP(verbose_TLIout+1, "%g__",Z[j][k]);
      MESSAGEP(verbose_TLIout+1, "<<\n");
      MESSAGEP(verbose_TLIout,"   Cross sections: ");
      pfwrite(CS[j],       sizeof(PREC_CS),        nT, fpout);
      MESSAGEP(verbose_TLIout+1, "Each CS: >>");
      for(int k=0; k<nT; k++)
        MESSAGEP(verbose_TLIout+1, "%g__",CS[j][k]);
      MESSAGEP(verbose_TLIout+1, "<<\n");
    }

    MESSAGEP(verbose_TLIout, " DB correlative number: ");
    pfwrite(&rdb, sizeof(unsigned short),  1, fpout);

    /* Free memory: */
    free(name);
    free(T);
    for(int j=0 ; j<niso ; j++)
      free(isonames[j]);
    free(isonames);
    free(mass);
    free(*Z);
    free(Z);
    free(*CS);
    free(CS);

    /* Cumulative number of isotopes per database */
    if(acum[i+1])
      acum[i+1] += niso;
    else
      acum[i+1]  = acum[i] + niso;

    // messagep(2,".");
    if (moredb) i--;
  }

  /* Print out info: */
  long ndbpos2 = ftell(fpout);
  fseek(fpout, ndbpos, SEEK_SET);
  MESSAGEP(verbose_TLIout, "Corrected Number of DBs:    "); 
  pfwrite(&rdb, sizeof(unsigned short), 1, fpout); 
  fseek(fpout, ndbpos2, SEEK_SET);

  MESSAGEP(verbose_TLIout, "Total number of isotopes: ");
  pfwrite(acum+ndb, sizeof(unsigned short), 1, fpout);
  MESSAGEP(verbose_TLIout, "--------------------------\n");

  messagep(2, "done.\n");

  return LR_OK;
}


/* \fcnfh
   Read line transitions, sort it, and write it on output file */
int
readwritetransition(unsigned short *acum, /* Cumul. number of isotopes/DB    */
                    double ini,           /* Lowest wavelength limit to read */
                    double fin,           /* Upper wavelength limit to read  */
                    double del){          /* Wavelength chunk sizes to read  */
  double wav1 = ini;                        /* Low wavelength limit of chunk */
  struct linedb *lineinfo[ndb], *curr[ndb]; /* Cumulative and current linedb */
  unsigned long nrec=0;                     /* Number of lines stored in TLI */

  messagep(2, "Reading and writing line transitions:\n");
  /* Read data in chunks between wav1 and wav2=wav1+del:  */
  while (wav1 < fin){
    double wav2 = wav1 + del;
    if (wav2>fin) wav2 = fin;
    long int ncurr[ndb]; /* Number of transitions read    */
    messagep(2, " Wavelength range %g to %g microns: reading...", wav1, wav2);

    unsigned short left=0;  /* Number of data bases with lines read */
    for (int i=0; i<ndb; i++){
      ncurr[left] = (*(drivers[DBdriver[i]]->info))(lineinfo+left, wav1, wav2);
      if (ncurr[left] == 0) /* No lines were read         */
        continue;
      curr[left] = lineinfo[left];
      left++; /* Count only of lines were read            */
    }

    messagep(2, " sorting & writing...");
    while (left){
      /* Search for the lowest-wavelength record:               */
      unsigned short mindb = 0;  /* Data base of min wavelength */
      double wmin = curr[0]->wl; /* Minimum wavelength          */
      for (int i=1; i<left; i++)
        if (curr[i]->wl < wmin)
          wmin = curr[mindb=i]->wl;

      /* Write lowest-wavelength transition: */
      unsigned short id = acum[mindb] + curr[mindb]->isoid;
      verbose_TLIout += 2;
      MESSAGEP(verbose_TLIout, "Wavelength:  %g - ", curr[mindb]->wl);
      pfwrite(&(curr[mindb]->wl), sizeof(PREC_LNDATA),    1, fpout);
      MESSAGEP(verbose_TLIout, " Isotope ID: ");
      pfwrite(&id,                sizeof(unsigned short), 1, fpout);
      MESSAGEP(verbose_TLIout, " Elow:       ");
      pfwrite(&curr[mindb]->elow, sizeof(PREC_LNDATA),    1, fpout);
      MESSAGEP(verbose_TLIout, " gf:         ");
      pfwrite(&curr[mindb]->gf,   sizeof(PREC_LNDATA),    1, fpout);
      verbose_TLIout -= 2;

      nrec++; /* Count number of lines records */

      curr[mindb]++; /* Advance pointer of data base that contained wmin  */
      /* If this was the last transition in this database, discard this
         database from sorting                                            */
      if(!--ncurr[mindb]){
        while(++mindb < left){
          MESSAGEP(verbose_TLIout,
                   "changing down %i (n: %li, c: %lx, a: %i)\n", mindb,
                   ncurr[mindb], (unsigned long)curr[mindb], acum[mindb]);
          /* Shift the data-base indices of those > mindb */
          ncurr[mindb-1] = ncurr[mindb];
          curr[mindb-1]  = curr[mindb];
          acum[mindb-1]  = acum[mindb];
        }
      left--;
      }
    }
    messagep(2," done\n");

    /* Free lineinfo for this wavelength range and go for the next one */
    for (int i=0; i<ndb; i++)
      free(lineinfo[i]);
    wav1 = wav2;
  }

  messagep(3, "\nSuccessfully written %li records in file '%s'.\n",
           nrec, outfilename);

  return LR_OK;
}
