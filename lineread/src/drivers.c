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
#include <version_lr.h>


#define pfwrite(...) do{if (!dry) fwrite(__VA_ARGS__);}while(0)

//Add here a content struct for each driver, which should also have an
//"extern" entry in lineread.h as well as its proto_ file included.
static const
driver_func (*drivers)[] = 
{
  driverf_debug
};

static unsigned short ndb = 0;
static _Bool dry = 0;
static unsigned short *DBdriver = NULL;


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
  unsigned short nfcn = sizeof(initdriver)/sizeof(driver_func);
  ndb = hint->ndb;
  dry = hint->dry;
  DBdriver = (unsigned short *)calloc(ndb, sizeof(unsigned short));


  /* finding the right driver for each database */
  find_alldrivers(hint, nfcn);


  /* Setting up drivers & output file */
  FILE *fpout = setdriversnoutput(hint);


  /* Reading & writing partition info */
  unsigned short niso[ndb];
  readwritepartition(niso, fpout);


  /* Reading & writing line transitions */
  readwritetransition(niso, fpout, hint->iniw, hint->finw, hint->delw);


  /* Closing up drivers & output file */
  for (int i=0 ; i<ndb ; i++)
    (*(drivers[DBdriver[i]]->close))();
  if (fpout != stdout && !dry)
    fclose(fpout);
  free(DBdriver);
  
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


  for (int i=0 ; i<nfcn ; i++){
    (*(drivers[i])->init)();
    char *name = drivers[i]->name;
    lennm += strlen(name) + 2;
    while (lennm >= allocnm)
      alldb = (char *)realloc(alldb, (allocnm<<=1)*sizeof(char));
    strcat(alldb, name);
    if(i<nfcn-1) strcat(alldb, ", ");
  }
  alldb = (char *)realloc(alldb, lennm*sizeof(char));


  for (int i=0 ; i<ndb ; i++)
    if((DBdriver[i] = hint->dbd[i] = find_dbd(hint->db[i], nfcn)) < 0)
      mperror(MSGP_USER,
	      "The file '%s' could not be associated to any "
	      "supported database.  Currently, lineread can read: %s."
	      , hint->db[i], alldb);
}


/*********************************************/


/* \fcnfh
   Set-up drivers and open output file
*/
FILE *
setdriversnoutput(struct hints *hint)
{
  FILE *fpout = NULL;
  for (int i=0 ; i<ndb ; i++)
    (*(drivers[DBdriver[i]]->open))(hint->db[i], hint->dbaux[i]);

  if(!dry){
    if(strcmp("-", hint->datafile) == 0)
      fpout = stdout;
    else if((fpout=fopen(hint->datafile, "w"))==NULL){
      mperror(MSGP_USER,
	      "Data file '%s' cannot be opened for writing.\n"
	      , hint->datafile);
    }
  }

  /* Magic number, which in big endian would be abb3b6ff =
    {(char)0xff-'T',(char)0xff-'L',(char)0xff-'I',(char)0xff},
    or in little endian: ffb6b3ab */
  long int magic=((0xff-'T')<<24)|((0xff-'L')<<16)|((0xff-'I')<<8)|(0xff);
  char *undefined_string = "";
  unsigned short rn = strlen(undefined_string);

  pfwrite(&magic,           sizeof(long int),        1, fpout);
  pfwrite(&version,         sizeof(unsigned short),  1, fpout); 
  pfwrite(&revision,        sizeof(unsigned short),  1, fpout);
  pfwrite(&(hint->iniw),    sizeof(double),          1, fpout);
  pfwrite(&(hint->finw),    sizeof(double),          1, fpout);
  pfwrite(&rn,              sizeof(unsigned short),  1, fpout); 
  pfwrite(undefined_string, sizeof(char),           rn, fpout);

  return fpout;
}


/**********************************************/



/* \fcnfh
   Read line transition, sort it, and write it on output file 
*/
int
readwritetransition(unsigned short *niso,
		    FILE *fpout,
		    double ini,
		    double fin,
		    double del)
{
  double wav1 = ini;
  struct linedb *lineinfo[ndb], *curr[ndb];

  while (wav1 < fin){
    double wav2 = wav1 + del;
    long int ncurr[ndb], icurr[ndb];
    unsigned short acum[ndb];

    //Get line info from the DBs
    for (int i=0 ; i<ndb ; i++){
      ncurr[i] = (*(drivers[DBdriver[i]]->info))(lineinfo+i, wav1, wav2);
      curr[i] = lineinfo[i];
      acum[i] = icurr[i] = 0;
    }
    //Set accumulated number of isotopes
    for (int i=1 ; i<ndb ; i++)
      for (int j=0 ; j<i ; j++)
	acum[i] += niso[j];

    unsigned short left=ndb;
    while (left){

      //Search for the lowest wavelength among the remaining DBs
      unsigned short mindb = 0;
      double wmin = curr[0]->wl;
      for (int i=1 ; i<left ; i++)
	if (curr[i]->wl < wmin)
	  wmin = curr[mindb=i]->wl;


      //Write lowest-wavelength transition
      unsigned short id = acum[mindb] + curr[mindb]->isoid;
      pfwrite(&(curr[mindb]->wl), sizeof(PREC_LNDATA),    1, fpout);
      pfwrite(&id,                sizeof(unsigned short), 1, fpout);
      pfwrite(&curr[mindb]->elow, sizeof(PREC_LNDATA),    1, fpout);
      pfwrite(&curr[mindb]->gf,   sizeof(PREC_LNDATA),    1, fpout);


      //If this was the last transition in this database, in this range,
      //then discard this database from sorting.
      if(++icurr[mindb] == ncurr[mindb]){
	while(++mindb < left){
	  icurr[mindb-1] = icurr[mindb];
	  ncurr[mindb-1] = ncurr[mindb];
	  curr[mindb-1]  = curr[mindb];
	  acum[mindb-1]  = acum[mindb];
	}
	left--;
      }
    }

    //free line information from this range and go for the next
    //wavelength range
    for (int i=0 ; i<ndb ; i++)
      free(lineinfo[i]);

    wav1 = wav2;
  }
}


/**********************************************/



/* \fcnfh
Read partition info and write it on output file 
*/
int
readwritepartition(unsigned short *niso,
		   FILE *fpout)
{
  for (unsigned short i=0 ; i<ndb ; i++){
    char *name, **isonames;
    unsigned short nT;
    PREC_MASS *mass;
    PREC_TEMP *T;
    PREC_Z    **Z;
    PREC_CS   **CS;
    (*(drivers[DBdriver[i]]->part))(&name, &nT, &T, niso+i, &isonames, &mass, &Z, &CS);

    unsigned short rn = strlen(name);
    pfwrite(&rn,    sizeof(unsigned short),  1, fpout); 
    pfwrite(name,   sizeof(char),           rn, fpout);
    pfwrite(&nT,    sizeof(unsigned short),  1, fpout);
    pfwrite(niso+i, sizeof(unsigned short),  1, fpout);
    pfwrite(T,      sizeof(PREC_TEMP),      nT, fpout);

    for(int j=0 ; j<niso[i] ; j++){
      rn = strlen(isonames[j]);
      pfwrite(&rn,         sizeof(unsigned short),  1, fpout);
      pfwrite(isonames[j], sizeof(char),           rn, fpout);
      pfwrite(mass+j,      sizeof(PREC_MASS),       1, fpout);
      pfwrite(Z[j],        sizeof(PREC_Z),         nT, fpout);
      pfwrite(CS[j],       sizeof(PREC_CS),        nT, fpout);
    }

    pfwrite(&i,     sizeof(unsigned short),  1, fpout);
  }
}
