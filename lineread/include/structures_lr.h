/*********************************************************************
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

#ifndef _STRUCTURES_LR_H_
#define _STRUCTURES_LR_H_

/* DO NOT change the following structures, not even the order of its
   components. Which is used in lineread:: main():: writing loop */

struct linedb{
  PREC_NREC recpos;         /* Record position in the file         */
  unsigned :0;              /* Just some padding                   */
  PREC_LNDATA wl;           /* Wavelength in nm.                   */ 
                            /* FINDME: in nm?? */
  PREC_LNDATA elow;         /* Lower energy level in cm-1          */
  PREC_LNDATA gf;           /* GF value                            */
  short isoid;              /* Isotope ID (Assumed to be in range) */
};

struct hints{
  int ndb;         /* Number of data bases                                */
  char **db;       /* Data-base filenames                                 */
  char **dbaux;    /* Auxiliar data-base filenames (partition function)   */
  int *dbd;        /* WHAT IS THIS?                                       */
  double iniw,     /* Lower-limit wavelength to read (microns)            */
         finw,     /* Higher-limit wavelength to read (microns)           */
         delw;     /* Range of wavelengths to read at each time (microns) */
  char *datafile;  /* Output TLI filename                                 */
  _Bool dry;       /* Is this a run without outputs?                      */
};

typedef struct driver_func{
  const char *name;
  _Bool (*find)(const char *name); /* Returns true if the driver can read
                                      the given file                        */
  int (*open)(char *dbname,
              char *dbaux);        /* Open the data base dbname and
                                      auxiliary daba base dbaux for reading */
  int (*close)();                  /* Close the driver and ree memory       */
  long int (*info)(struct linedb **lineinfo,
                   double wav1,
                   double wav2);   /* Read from wav1 to wav2, allocate and
                                      store in lineinfo.  Returns number of
                                      fields read.                          */
  _Bool (*part)(char **name,
                unsigned short *nT,
                PREC_TEMP **T,
                unsigned short *niso,
                char ***isonames,
                PREC_MASS **mass,
                PREC_Z ***Z,
                PREC_CS ***CS);    /* Read partition function information
                                      from database and return the required
                                      components.  Returns true while more
                                      databases need to be read.            */
} driver_func;

#endif /* _STRUCTURES_LR_H_ */
