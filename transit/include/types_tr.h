/*
 * types_tr.h - types for the transit program.
 *
 * Copyright (C) 2003-2006 Patricio Rojo (pato@astro.cornell.edu)
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

#ifndef _TYPES_TR_H
#define _TYPES_TR_H

/*  Types  */
#define PREC_NSAMP int      /* Type for radius and wavelength indices */
#define PREC_NREC long      /* Type for record indices                */
#define PREC_ZREC double    /* Type for the partition info            */
#define PREC_LNDATA double  /* Type for the line data output          */
#define PREC_RES double     /* Type for every partial result          */
#define PREC_ATM double	    /* Type for atmospheric data              */
#define PREC_CS  double	    /* Type for cross-section                 */
#define PREC_CIA float      /* Type for collision induced absorption  */

#endif /* _TYPES_TR_H */
