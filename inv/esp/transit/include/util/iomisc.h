/*
 * iomisc.h - Miscellaneous input/output utilities header file
 *
 * Copyright (C) 2004 Patricio Rojo (pato@astro.cornell.edu)
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

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <termios.h>
#include <stdarg.h>
#include <string.h>

#define __PADTO_COLUMNS 78

extern char *linepad_break;
extern int printpad_columns;
extern char printpad_endpadchar;

#define readd(fp,cp) readds(fp,cp,NULL,0)

#include <proto_iomisc.h>
