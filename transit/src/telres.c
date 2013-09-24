/*
 * telres.c    - Convolves to telescope resolution.  Component of the
 *               Transit program. 
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

static inline PREC_ZREC interpolateTR(PREC_ZREC *x,
				      PREC_ZREC *y,
				      int n,
				      PREC_ATM val);

static inline int makegauTR(int *nsmr,
			    double **smr,
			    double width,
			    double dl,
			    float mgau);

static inline double integrateTR(PREC_RES *spectra,
				 int osw);



