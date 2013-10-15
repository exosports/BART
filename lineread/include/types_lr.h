/*
 * types_lr.h - types for the lineread program.
 *
 * Copyright (C) 2006 Patricio Rojo (pato@astro.cornell.edu)
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

#ifndef _TYPES_LR_H
#define _TYPES_LR_H

/*  Types  */
#define PREC_NREC   long    /* Type for record indices       */
#define PREC_ZREC   double  /* Type for the partition info   */
#define PREC_LNDATA double  /* Type for the line data output */
#define PREC_RES    double  /* Type for every partial result */
#define PREC_Z      double  /* Type for WHAT?                */
#define PREC_CS     double  /* Type for cross-section        */
#define PREC_TEMP   double  /* Type for temperature          */
#define PREC_MASS   double  /* Type for isotope mass         */

#endif /* _TYPES_LR_H */
