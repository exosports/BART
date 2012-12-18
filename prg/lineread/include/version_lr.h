/*********************************************************************
 * 
 * version_lr.h - Header that contains the version of lineread.
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

#ifndef _VERSION_LR_H_
#define _VERSION_LR_H_


#ifdef _VERSION_ONLY_DECLARATIONS_

extern unsigned short version;
extern unsigned short revision;
extern int version_rc;

extern int TLIversion;

#else  /*  _VERSION_ONLY_DECLARATIONS_ */

unsigned short version = 4;
unsigned short revision = 0;
int version_rc = 0;

int TLIversion = 4;

#endif  /*  _VERSION_ONLY_DECLARATIONS_ */


#endif /* _VERSION_LR_H_ */
