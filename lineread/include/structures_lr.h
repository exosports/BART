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
  PREC_NREC recpos;         //Record position in the file
  unsigned :0;              //Just some padding
  PREC_LNDATA wl;           //Wavelength in nm.
  PREC_LNDATA elow;         //Lower energy level in cm-1
  PREC_LNDATA gf;           //GF value
  short isoid;              //Isotope ID (Assumed to be in range)
};


#endif /* _STRUCTURES_LR_H_ */
