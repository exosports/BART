/*
 * at_common.c - Read atmospheric info common files. Component of the
 *                Transit program.
 *
 * Copyright (C) 2004-2006 Patricio Rojo (pato@astro.cornell.edu)
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

#include <readatm.h>


/* \fcnfh
   Skip to the next field pointed by 'lp', omits leading and trailing
   spaces.

   @returns pointer to the beginning of next non-space
*/
inline static char *
nextfield(char *lp)
{
  while(*lp==' '||*lp=='\t') lp++;
  while(*lp!=' '&&*lp!='\0') lp++;
  while(*lp==' '||*lp=='\t') lp++;

  return lp;
}


