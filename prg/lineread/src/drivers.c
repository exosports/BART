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

#include <lineread.h>
#define _VERSION_ONLY_DECLARATIONS_
#include <version_lr.h>

//Add here a content struct for each driver, which should also have an
//"extern" entry in lineread.h as well as its proto_ file included.
static 
driver_func * (*init_driver[])() = 
{
  &initdb_debug,
  //  &initdb_hitran4,
  &initdb_pands,
  &initdb_text
};

static driver_func **drivers=NULL;

static char *outfilename = NULL;
static unsigned short ndb;
static _Bool dry;
static unsigned short *DBdriver = NULL;
static FILE *fpout = NULL;
static int verbose_TLIout = 10;
static int verbose_TLIout2 = 5;


static inline void
pfwrite(void *pointer,
	size_t size,
	size_t nmemb,
	FILE *fp)
{
  if(size*nmemb)
    MESSAGEP(verbose_TLIout, "%04hx%s\n",
	     *((unsigned short *)pointer), size==sizeof(unsigned short)?"":" ..."); 
  if(!dry)
    fwrite(pointer, size, nmemb, fp);
}


/* \fcnfh
   Finds the DB driver appropriate to read file 'file', or -1 if none is
   found.
*/
static inline int
find_dbd(char *file,
	 unsigned short nfcn)
{
  for (int i=0 ; i<nfcn ; i++)
    if((*(drivers[i]->find))(file))
      return i;
  return -1;
}




/* \fcnfh
 Returns the name of the ith kind of database
*/
int
db_drivers(struct hints *hint)
{
  unsigned short nfcn = sizeof(init_driver)/sizeof(driver_func *);
  ndb = hint->ndb;
  dry = hint->dry;
  DBdriver = (unsigned short *)calloc(ndb, sizeof(unsigned short));


  /* finding the right driver for each database */
  find_alldrivers(hint, nfcn);


  /* Setting up drivers & output file */
  setdriversnoutput(hint);

  
  /* Reading & writing partition info */
  unsigned short acum[ndb+1];
  readwritepartition(acum);


  /* Reading & writing line transitions */
  readwritetransition(acum, hint->iniw, hint->finw, hint->delw);

  return LR_OK;
}


/*********************************************/


/* \fcnfh
   Close drivers & output file
*/
void  
drivers_free(struct hints *hint)
{
  for (int i=0 ; i<hint->ndb ; i++)
    (*(drivers[DBdriver[i]]->close))();
  if(DBdriver)
    free(DBdriver);
  if (fpout != NULL && fpout != stdout && !dry)
    fclose(fpout);
  if(drivers)
    free(drivers);
  if(outfilename)
    free(outfilename);
}


/*********************************************/


/* \fcnfh
   Set-up drivers and open output file
*/
int
find_alldrivers(struct hints *hint,
		unsigned short nfcn)
{
  int allocnm = 8, lennm=0;
  char *alldb = (char *)calloc(allocnm, sizeof(char));

  drivers = (driver_func **)calloc(nfcn, sizeof(driver_func *));

  for (int i=0 ; i<nfcn ; i++){
    drivers[i] = (*(init_driver[i]))();
    char *name = (char *)drivers[i]->name;
    lennm += strlen(name) + 2;
    while (lennm >= allocnm)
      alldb = (char *)realloc(alldb, (allocnm<<=1)*sizeof(char));
    strcat(alldb, name);
    if(i<nfcn-1) strcat(alldb, ", ");
  }
  alldb = (char *)realloc(alldb, lennm*sizeof(char));

  messagep(2, "Finding drivers for %i databases", ndb);
  messagep(3, ":\n");
  int result;
  for (int i=0 ; i<ndb ; i++){
    if((result = find_dbd(hint->db[i], nfcn)) < 0){
      mperror(MSGP_USER|MSGP_ALLOWCONT,
	      "The file '%s' could not be associated to any "
	      "supported database.  Currently, lineread can read: %s.\n"
	      , hint->db[i], alldb);
      lineread_free();
      free(alldb);
      exit(EXIT_FAILURE);
    }
    messagep(3, " %s: found (%s)", hint->db[i], drivers[result]->name);
    messagep(2, ".");
    messagep(3, "\n");
    DBdriver[i] = hint->dbd[i] = result;
  }
  free(alldb);
  if(verblevel <3) messagep(2, " ");
  messagep(2, "done\n");

  return LR_OK;
}


/*********************************************/


/* \fcnfh
   Set-up drivers and open output file
*/
int
setdriversnoutput(struct hints *hint)
{
  outfilename = strdup(strncmp("-",hint->datafile,2)==0?
		       "Standard Output":hint->datafile);


  for (int i=0 ; i<ndb ; i++)
    (*(drivers[DBdriver[i]]->open))(hint->db[i], hint->dbaux[i]);

  if(!dry){
    if(strcmp("-", hint->datafile) == 0)
      fpout = stdout;
    else if((fpout=fopen(hint->datafile, "w"))==NULL){
      mperror(MSGP_USER|MSGP_ALLOWCONT,
	      "Data file '%s' cannot be opened for writing.\n"
	      , hint->datafile);
      lineread_free();
      exit(EXIT_FAILURE);
    }
  }
  else
    messagep(2, "Dry-");
  messagep(2, "Opened input DB%s (%i), output stream (%s),\n "
	   "and writing TLI header... "
	     , ndb>0?"s":"", ndb, outfilename); 

  /* Magic number, which in big endian would be abb3b6ff =
    {(char)0xff-'T',(char)0xff-'L',(char)0xff-'I',(char)0xff},
    or in little endian: ffb6b3ab */
  int32_t magic=((0xff-'T')<<24)|((0xff-'L')<<16)|((0xff-'I')<<8)|(0xff);
  char *undefined_string = "";
  unsigned short rn = strlen(undefined_string);

  MESSAGEP(verbose_TLIout2, "Magic number:     "); 
  pfwrite(&magic,           sizeof(int32_t),        1, fpout);
  MESSAGEP(verbose_TLIout2, "TLIVersion:       "); 
  pfwrite(&TLIversion,      sizeof(unsigned short),  1, fpout); 
  MESSAGEP(verbose_TLIout2, "Version:          "); 
  pfwrite(&version,         sizeof(unsigned short),  1, fpout); 
  MESSAGEP(verbose_TLIout2, "Revision:         "); 
  pfwrite(&revision,        sizeof(unsigned short),  1, fpout);
  MESSAGEP(verbose_TLIout2, "IniW (%7g):   ", hint->iniw); 
  pfwrite(&(hint->iniw),    sizeof(double),          1, fpout);
  MESSAGEP(verbose_TLIout2, "FinW (%7g):   ", hint->finw); 
  pfwrite(&(hint->finw),    sizeof(double),          1, fpout);
  MESSAGEP(verbose_TLIout2, "Length undefined: "); 
  pfwrite(&rn,              sizeof(unsigned short),  1, fpout); 
  MESSAGEP(verbose_TLIout2, "Undefined:        %s", rn==0?"\n":""); 
  pfwrite(undefined_string, sizeof(char),           rn, fpout);

  messagep(2, "done\n");

  return 1;
}


/**********************************************/



/* \fcnfh
Read partition info and write it on output file 
*/

int
readwritepartition(unsigned short *acum)
{
  memset(acum, 0, (ndb+1)*sizeof(unsigned short));
  int rdb=0;

  messagep(2, "Reading and writing partition information");
  long ndbpos = ftell(fpout);
  MESSAGEP(verbose_TLIout, "Number of DBs:    "); 
  pfwrite(&ndb,  sizeof(unsigned short), 1, fpout); 

  //rdb is the "real DB" counter, while i is the counter of the number
  //of files that are read.
  for (unsigned short i=0 ; i<ndb ; i++, rdb++){
    char *name, **isonames;
    unsigned short nT, niso;
    PREC_MASS *mass;
    PREC_TEMP *T;
    PREC_Z    **Z;
    PREC_CS   **CS;
    _Bool moredb = (*(drivers[DBdriver[i]]->part))
      (&name, &nT, &T, &niso, &isonames, &mass, &Z, &CS);

    MESSAGEP(verbose_TLIout," For DB #%i (file #%i):\n", rdb, i);
    unsigned short rn = strlen(name);
    MESSAGEP(verbose_TLIout," Length of name:         ");
    pfwrite(&rn,    sizeof(unsigned short),  1, fpout); 
    MESSAGEP(verbose_TLIout," Name:                   ");
    pfwrite(name,   sizeof(char),           rn, fpout);
    MESSAGEP(verbose_TLIout," Number of temperatures: ");
    pfwrite(&nT,    sizeof(unsigned short),  1, fpout);
    MESSAGEP(verbose_TLIout," Number of isotopes:     ");
    pfwrite(&niso, sizeof(unsigned short),  1, fpout);
    MESSAGEP(verbose_TLIout," Temperatures:           ");
    pfwrite(T,      sizeof(PREC_TEMP),      nT, fpout);
    MESSAGEP(verbose_TLIout+1, "Each T: >>");
    for(int k=0; k<nT ; k++)
      MESSAGEP(verbose_TLIout+1, "%g__",T[k]);
    MESSAGEP(verbose_TLIout+1, "<<\n");
    
    for(int j=0 ; j<niso ; j++){
      MESSAGEP(verbose_TLIout,"  For isotope %i:\n", j);
      rn = strlen(isonames[j]);
      MESSAGEP(verbose_TLIout,"   Length of name: ");
      pfwrite(&rn,         sizeof(unsigned short),  1, fpout);
      MESSAGEP(verbose_TLIout,"   Name:           ");
      pfwrite(isonames[j], sizeof(char),           rn, fpout);
      MESSAGEP(verbose_TLIout,"   Masses:         ");
      pfwrite(mass+j,      sizeof(PREC_MASS),       1, fpout);
      MESSAGEP(verbose_TLIout,"   Partition:      ");
      pfwrite(Z[j],        sizeof(PREC_Z),         nT, fpout);
      MESSAGEP(verbose_TLIout+1, "Each Z: >>");
      for(int k=0; k<nT ; k++)
	MESSAGEP(verbose_TLIout+1, "%g__",Z[j][k]);
      MESSAGEP(verbose_TLIout+1, "<<\n");
      MESSAGEP(verbose_TLIout,"   Cross sections: ");
      pfwrite(CS[j],       sizeof(PREC_CS),        nT, fpout);
      MESSAGEP(verbose_TLIout+1, "Each CS: >>");
      for(int k=0; k<nT ; k++)
	MESSAGEP(verbose_TLIout+1, "%g__",CS[j][k]);
      MESSAGEP(verbose_TLIout+1, "<<\n");
    }

    MESSAGEP(verbose_TLIout, " DB correlative number: ");
    pfwrite(&rdb, sizeof(unsigned short),  1, fpout);

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

    //Set accumulated number of  for a multi- or single-DB read
    if(acum[i+1])
      acum[i+1] += niso;
    else
      acum[i+1]  = acum[i] + niso;

    messagep(2,".");
    if (moredb) i--;
  }

  long ndbpos2 = ftell(fpout);
  fseek(fpout, ndbpos, SEEK_SET);
  MESSAGEP(verbose_TLIout, "Corrected Number of DBs:    "); 
  pfwrite(&rdb,  sizeof(unsigned short), 1, fpout); 
  fseek(fpout, ndbpos2, SEEK_SET);


  MESSAGEP(verbose_TLIout, "Total number of isotopes: ");
  pfwrite(acum+ndb, sizeof(unsigned short), 1, fpout);
  MESSAGEP(verbose_TLIout, "--------------------------\n");

  messagep(2, " done\n");

  return LR_OK;
}


/**********************************************/



/* \fcnfh
   Read line transition, sort it, and write it on output file 
*/
int
readwritetransition(unsigned short *acum,
		    double ini,
		    double fin,
		    double del)
{
  double wav1 = ini;
  struct linedb *lineinfo[ndb], *curr[ndb];
  unsigned long nrec=0;

  messagep(2, "Reading and writing line transitions:\n");

  while (wav1 < fin){
    double wav2 = wav1 + del;
    if (wav2>fin) wav2 = fin;
    long int ncurr[ndb];

    messagep(2, " Wavelength range %g to %g microns: reading..."
	     , wav1, wav2);
    //Get line info from the DBs
    unsigned short left=0;
    for (int i=0 ; i<ndb ; i++){
      ncurr[left] = (*(drivers[DBdriver[i]]->info))(lineinfo+left, wav1, wav2);
      if (ncurr[left] == 0)
	continue;
      curr[left] = lineinfo[left];
      left++;
    }

    messagep(2, " sorting & writing...");
    while (left){

      //Search for the lowest wavelength among the remaining DBs
      unsigned short mindb = 0;
      double wmin = curr[0]->wl;
      for (int i=1 ; i<left ; i++)
	if (curr[i]->wl < wmin)
	  wmin = curr[mindb=i]->wl;


      //Write lowest-wavelength transition
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

      nrec++;

      //If this was the last transition in this database, in this range,
      //then discard this database from sorting.
      curr[mindb]++;
      if(!--ncurr[mindb]){
	while(++mindb < left){
	  MESSAGEP(verbose_TLIout,
		   "changing down %i (n: %li, c: %lx, a: %i)\n"
		   , mindb, ncurr[mindb], (unsigned long)curr[mindb], acum[mindb]);
	  ncurr[mindb-1] = ncurr[mindb];
	  curr[mindb-1]  = curr[mindb];
	  acum[mindb-1]  = acum[mindb];
	}
	left--;
      }
    }
    messagep(2," done\n");

    //free line information from this range and go for the next
    //wavelength range
    for (int i=0 ; i<ndb ; i++)
      free(lineinfo[i]);

    wav1 = wav2;
  }

  messagep(3, "\nSuccessfully written %li records in file '%s'.\n"
	   , nrec, outfilename);

  return LR_OK;
}





