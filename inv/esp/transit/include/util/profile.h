/*
 * profile.h - Spectral line shape functions
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


#ifndef _PROFILE_H
#define _PROFILE_H

#define PREC_VOIGT float
#define PREC_VOIGTP double

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define VOIGT_QUICK 0x00001   //Quick integration.
extern int _voigt_maxelements;

#include <proto_voigt.h>

#endif /* _PROFILE_H */

