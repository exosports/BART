/*
 * lineread.h - header for lineread program related sorces.
 *              Part of Transit program.
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

#ifndef _LINEREAD_H_
#define _LINEREAD_H_

/* DO NOT change the following structures, not even the order of its
   components. Which is used in lineread:: main():: writing loop */

extern struct linedb{
  PREC_NREC recpos;         //Record position in the file
  unsigned :0;              //Just some padding
  PREC_LNDATA wl;           //Wavelength in nm.
  PREC_LNDATA elow;         //Lower energy level in cm-1
  PREC_LNDATA lgf;          //Log(gf) value
  short isoid;              //Isotope ID (Assumed to be in range)
} linedb;

short gabby_dbread=0;

#include <lineread_proto.h>

#ifdef LINEFCN_RUNTIME
/* To be implemented: Run-time check of available drivers.*/

int check_nlinefcn(PREC_NREC *(*(*linefcn))(),char ***name);

#else

/*
  dbread_*: Reading of different line databases.

  @returns number of lines read on success, or
           -1 if data file 'filename' cannot be opened. 
           -2 if data file cannot be stat()ed
           -3 if the file doesn't contain an integer number of records
	   -4 if initial wavelength is too long.
	   -5 if final wavelength is too short.
           -6 if memory cannot be allocated
*/

#if 0
#define LINEREAD_DEF(database)                                           \
PREC_NREC dbread_ ## database(                                           \
                    /* If NULL then default value is used */             \
                       char *filename,                                   \
		    /* 2 pointers level in order to be able to allocate  
		       memory */                                         \
		       struct linedb **lines,                            \
		    /* wavelengths in nanometer */                       \
		       float wlbeg,                                      \
		       float wlend,                                      \
	            /* Filename for Partition function data, if NULL
		       then default value is used */                     \
		       char *Zfilename,                                  \
		    /* For the following 3 parameter, the memory is
		       allocated in the dbread_* functions, and the
		       size is returned in the last parameters. That
		       also implies an extra level of pointers */        \
		    /* Partition function(isotope, temperature) */       \
		       PREC_ZREC ***Z,                                   \
		    /* Temperatures where Z was sampled. */              \
		       PREC_ZREC **T,                                    \
		    /* Isotopes mass in AMUs */                          \
		       PREC_ZREC **isomass,                              \
		    /* number of temperature points */                   \
		       int *nT,                                          \
		    /* number of isotopes */                             \
		       int *nIso,                                        \
		    /* Names of isotopes */                              \
                       char ***isonames)

/*
  For every extra function you want lineread to read: Update
  DBREAD_NFCN, add a LINEREAD_DEF() line and an entry in (*linefcn[]) 
*/

LINEREAD_DEF(pands);
#endif
#define dbread_nfcn 1
PREC_NREC (*linefcn[dbread_nfcn])(char *,struct linedb **,float,float
				  ,char *, PREC_ZREC ***,PREC_ZREC **
				  ,PREC_ZREC **, int *, int *,char ***)={
				    dbread_pands
				  };


#endif

#endif


