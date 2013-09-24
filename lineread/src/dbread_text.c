/*
 * dbread_text.c - Driver to read line information from a text file.
 *              Part of lineread program.
 *
 * Copyright (C) 2005-2006 Patricio Rojo (pato@astro.cornell.edu)
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

#define DRIVERNAME "text"

static long int currline = 0;
static long maxline=300;
static unsigned short ndb;
static char *dbfilename;



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
  mperror(MSGP_USER|MSGP_ALLOWCONT,
	       "Line %i of file '%s': Field %i (%s) has\n"
	       " not a valid value:\n%s\n"
	       ,nmb,file,fld,fldn,line);
  return -5;
}


/*****************************************************/


/*************************
 *************************
 *************************/

#define checkprepost(pointer,pre,omit,post) do{                       \
   if(pre)                                                            \
     mperror(MSGP_USER,                                               \
             "Pre-condition failed on line %i(%s)\n while reading:\n" \
	     "%s\n\nTLI_Ascii format most likely invalid\n"           \
             ,__LINE__,__FILE__,line);                                \
   while(omit)                                                        \
     pointer++;                                                       \
   if(post)                                                           \
     mperror(MSGP_USER,                                               \
             "Post-condition failed on line %i(%s)\n while reading:\n"\
	     "%s\n\nTLI_Ascii format most likely invalid\n"           \
             ,__LINE__,__FILE__,line);                                \
                                             }while(0)

/* \fcnfh
   It outputs error. Used when EOF is found before expected
*/
static void
earlyend(char *file, long lin)
{
  mperror(MSGP_USER|MSGP_ALLOWCONT,
	       "readlineinfo:: EOF unexpectedly found at line %i in\n"
	       "ascii-TLI linedb info file '%s'\n"
	       ,lin,file);
  lineread_free();
  exit(EXIT_FAILURE);
}


void 
linetoolong_text(int max,		/* Maxiumum length of an accepted line
					 */ 
		 char *file,		/* File from which we were reading */
		 int line)		/* Line who was being read */

{
  char fname[strlen(file)+1];
  strcpy(fname, file);
  lineread_free();
  linetoolong(max, fname, line);
}


static FILE *fp = NULL;
static _Bool partitionread = 0;
static int verbose_db = 15;


/* \fncfh
   If this driver can read name, then return 1 otherwise 0
*/
static _Bool
db_find(const char *name)
{
  FILE *fp;
  if((fp=fopen (name, "r")) != NULL){
    int maxlen=50;
    char line[maxlen];
    const char *id = "#TLI-ASCII";

    if(fgets(line, maxlen-1, fp) && 
       strncmp(line, id, strlen(id)) == 0){
      fclose(fp);
      return 1;
    }
    fclose(fp);
  }

  return 0;
}


/********************************************************/


/* \fcnfh
   Open database and auxiliary file, set all FILE pointers.
*/
static int
db_open(char *dbfilenameo, 
	char *dbaux)
{
  dbfilename = strdup(dbfilenameo);

  //Open database
  if((fp = fopen(dbfilename,"r")) == NULL)
    mperror(MSGP_USER,
	    "Could not open file '%s' for reading\n"
	    , dbfilename);

  int maxlen=50;
  char line[maxlen], *lp;
  const char *id = "#TLI-ASCII";

  if(fgets(line, maxlen-1, fp) && 
     strncmp(line, id, strlen(id)) != 0)
    mperror(MSGP_SYSTEM,
	    "File '%s' does not have the proper "
	    "TLI-ASCII heading, but it was approved "
	    "by db find(?)\n"
	    , dbfilename);

  currline = 1;
  size_t pos = ftell(fp);
  char rc;
  settoolongerr(&linetoolong_text,dbfilename,&currline);


  //skip comments and blank lines
  while((rc=fgetupto(line,maxline,fp)) == '#' || 
	rc == '\n')
    currline++;
  if(!rc) earlyend(dbfilename, currline);

  //get number of database which needs to be one at this point. If
  //omitted number of databases, it is assumed to be 1
  if (line[0] == 'd'){
    ndb = strtol(line+1, &lp, 0);
    if (lp == line+1) 
      mperror(MSGP_USER,
	      "Invalid number of databases in TLI-ASCII '%s'\n"
	      , dbfilename);
  }
  else{
    currline = ndb = 1;
    fseek(fp, pos, SEEK_SET);
  }

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

  free(dbfilename);
  freetoolongerr();

  return LR_OK;
}


/********************************************************/

/* \fcnfh
   Reads partition info 
*/
static _Bool
db_part(char **name,
	unsigned short *nT,
	PREC_TEMP **Tp,
	unsigned short *niso,
	char ***isonamesp,
	PREC_MASS **massp,
	PREC_Z ***Zp,
	PREC_CS ***CSp)
{
  char line[maxline], *lp, *lp2, rc;


  //skip comments and blank lines
  while((rc=fgetupto(line,maxline,fp)) == '#' || 
	rc == '\n')
    currline++;
  if(!rc) earlyend(dbfilename, currline);


  //read name, number of temps, and number of isotopes
  if((*name = readstr_sp_alloc(line, &lp, '_')) == NULL)
    mpallocerror(0);
  checkprepost(lp, 0, *lp==' '||*lp=='\t',*lp=='\0');
  int rn;
  rn = getnl(2, ' ', lp, niso, nT);
  checkprepost(lp, rn!=2, 0, 0);

  //allocate temperature, mass, cross section and partition info
  PREC_CS **CS;
  PREC_Z  **Z;
  char **inames;
  PREC_MASS *mass;
  PREC_TEMP *T;
  inames = *isonamesp = (char     **)calloc(*niso,    sizeof(char *));
  mass   = *massp     = (PREC_MASS *)calloc(*niso,    sizeof(PREC_ZREC));
  T      = *Tp        = (PREC_TEMP *)calloc(*nT,      sizeof(PREC_ZREC));
  CS     = *CSp       = (PREC_CS  **)calloc(*niso,    sizeof(PREC_CS *));
  Z      = *Zp        = (PREC_Z   **)calloc(*niso,    sizeof(PREC_ZREC *));
  *Z     = (PREC_Z    *)calloc(*niso**nT,sizeof(PREC_ZREC));
  *CS    = (PREC_CS   *)calloc(*niso**nT,sizeof(PREC_CS));

  //get isotope name and mass
  while((rc=fgetupto(lp2=line,maxline,fp)) == '#' || rc == '\n')
    currline++;
  if(!rc) earlyend(dbfilename, currline);
  for(int i=0 ; i<*niso ; i++){
    Z[i]  = Z[0]  + *nT*i;
    CS[i] = CS[0] + *nT*i;
    if((inames[i]=readstr_sp_alloc(lp2,&lp,'_')) == NULL)
      mpallocerror(0);
    //get mass
    mass[i] = strtod(lp,&lp2);

    if(i == *niso-1)
      break;

    checkprepost(lp2, lp==lp2, *lp2==' '||*lp2=='\t',*lp2=='\0');
  }
  //Last isotope has to be followed by an end of string
  checkprepost(lp2, 0, *lp2==' '||*lp2=='\t', *lp2!='\0');

  //get for each temperature a line with temperature, partfcn and
  //cross-sect info
  for(long t=0 ; t<*nT ; t++){
    while((rc=fgetupto(lp=line, maxline, fp)) == '#' || 
	  rc == '\n')
      currline++;
    if (!rc) earlyend(dbfilename, currline);
    while (*lp==' ') lp++;
    T[t] = strtod(lp, &lp);
    checkprepost(lp, *lp=='\0', *lp==' '||*lp=='\t', *lp=='\0');

    for(unsigned short i=0 ; i<*niso ; i++){
      Z[i][t] = strtod(lp, &lp);
      checkprepost(lp, *lp=='\0', *lp==' '||*lp=='\t', *lp=='\0');
    }
    unsigned short i;
    for(i=0 ; i<*niso-1 ; i++){
      CS[i][t] = strtod(lp, &lp);
      checkprepost(lp, *lp=='\0', *lp==' '||*lp=='\t', *lp=='\0');
    }
    CS[i][t] = strtod(lp, &lp);
    checkprepost(lp, 0, *lp==' '||*lp=='\t', *lp!='\0');
  }

  partitionread = 1;

  if (--ndb) return 1;
  else return 0;
}


/********************************************************/


/* \fcnfh
   Reads transition info from DB

   @returns -1 if db_part has not been called before
*/
static long int
db_info(struct linedb **lineinfop,
	double wav1,
	double wav2)
{
  if (!partitionread)
    return -1;

  MESSAGEP(verbose_db, "Driver: Going to look for wavelength range %g - %g\n"
	   , wav1, wav2);
  char line[maxline], *lp, *lp2;

  long int i=0, alloc=8;
  struct linedb *lineinfo =  (struct linedb *)calloc(alloc, 
						     sizeof(struct linedb));

  size_t pos = ftell(fp);
  long posline = currline;
  do{
    char rc;
    while((rc=fgetupto(lp=line,maxline,fp))=='#' || 
	  rc=='\n')
      currline++;
    //if it is not end of file, read the records.
    if(!rc)
      break;

    currline++;
    if (i==alloc)
      lineinfo = (struct linedb *)realloc(lineinfo, (alloc<<=1) * 
					  sizeof(struct linedb));
    double wavl = strtod(lp, &lp2);
    if(lp==lp2)
      return invalidfield(line, dbfilename, currline
			  , 1, "central wavelength");
    if (wavl < wav1) continue;
    if (wavl > wav2){
      fseek(fp, pos, SEEK_SET);
      currline = posline;
      break;
    }

    lineinfo[i].recpos = i;
    lineinfo[i].wl     = wavl;
    lineinfo[i].isoid  = strtol(lp2, &lp, 0);
    if(lp==lp2)
      return invalidfield(line, dbfilename, currline
			  , 2, "isotope ID");
    lineinfo[i].elow=strtod(lp, &lp2);
    if(lp==lp2)
      return invalidfield(line, dbfilename, currline
			  , 3, "lower energy level");
    lineinfo[i++].gf=strtod(lp2, &lp);
    if(lp==lp2)
      return invalidfield(line, dbfilename, currline
			  , 4, "log(gf)");

    pos = ftell(fp);
    posline = currline;
  }while(1);

  *lineinfop  = (struct linedb *)realloc(lineinfo, 
					 i*sizeof(struct linedb));
 
  return i;
}


/********************************************************/


static const driver_func pdriverf_text = {
  "TLI-ASCII driver",
  &db_find,
  &db_open,
  &db_close,
  &db_info,
  &db_part
};


driver_func *
initdb_text()
{
  return (driver_func *)&pdriverf_text;
}

