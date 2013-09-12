/*
 * lineread.h - header for lineread program related sorces.
 *              Part of Transit program.
 *
 * Copyright (C) 2003 Patricio Rojo (pato@astro.cornell.edu)
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
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

#ifndef _LINEREAD_H_
#define _LINEREAD_H_

#include <math.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <endian.h>
#include <sys/stat.h>
#include <string.h>
#include <ctype.h>
#include <pu/iomisc.h>

#include <messagep.h>

#include <flags_lr.h>
#include <types_lr.h>
#include <structures_lr.h>


extern short gabby_dbread;
extern double tli_fct;
extern struct linedb linedb;




#include <proto_lineread.h>
#include <proto_argum.h>
#include <proto_drivers.h>

#include <proto_dbread_debug.h>
#include <proto_dbread_pands.h>
#include <proto_dbread_text.h>
//#include <proto_dbread_hitran4.h>

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

#endif


