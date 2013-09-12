/*********************************************************************
 * 
 * dbread_template.c - Template driver to read a database with lineread
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


static FILE *fp = NULL, *fpaux = NULL;
static _Bool partitionread = 0;
static int verbose_db = 15;

/* \fncfh
   If this driver can read name, then return 1 otherwise 0
*/
static _Bool
db_find(const char *name)
{
  "Verify that name is readable";

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

  //Optionally open auxiliary
  if (dbaux)
    "Open Auxiliary file";

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

  return LR_OK;
}


/********************************************************/

/* \fcnfh
   Reads partition info.  Returns true if there are more databases to read.
*/
static _Bool
db_part(char **name,
	unsigned short *nT,
	PREC_TEMP **T,
	unsigned short *niso,
	char ***isonames,  	/*  Each name on its individually
				    allocated string */
	prec_MASS **mass,
	PREC_Z ***Z,
	PREC_CS ***CS)
{

  "Read all partition info";

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

  messagep(verbose_db, "Driver: Going to look for wavelength range %g - %g\n"
	   , wav1, wav2);
  int maxline = 200;
  char line[maxline], *lp;

  long int i=0, alloc=8;
  *lineinfo = (struct linedb *)calloc(alloc, 
				      sizeof(struct linedb));

  "Read line information";
  
  return i;
}


/********************************************************/



static const driver_func pdriverf_TEMPLATE = {
  "TEMPLATE driver",
  &db_find,
  &db_open,
  &db_close,
  &db_info,
  &db_part
};

driver_func *
initdb_TEMPLATE()
{
  return (driver_func *)&pdriverf_TEMPLATE;
}
