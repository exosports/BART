/*
 * version_tr.h - Header containing version information with
 *                variable initialization.
 *
 * Copyright (C) 2003 Patricio Rojo (pato@astro.cornell.edu)
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

#ifndef _VERSION_TR_H
#define _VERSION_TR_H

/* Version history:
   0.3: First light, it only calculates Kappa values of it and of
        lorentz widening seem rather high. 102703. Patricio Rojo
   0.5: Lorentz and doppler width work all right now. 102803. PR
   1.0: One dimensional extinction coefficient working with input from
        command line or atmopshere file. 032504. PMR
   1.1: Multi radius working for extinction calculation. 032804. PMR
   1.2: (1.1 patch) Exportable and fixes on one P,T point. 033104. PMR
   2.0: Modulation is fully working, not much debugging chase done,
        though. 051804. PMR
   2.1: Tau seems to be completely debugged, spline integration
        works. 070804. PMR
   2.2: Modulation seems to be completely debugged, spline integration
        works. 071504. PMR
   2.3: NEVER COMPLETED. Only rc1 was out when important changes in
        reading atmospheric file were done, hence I decided to jump
        straight into version 3.0. 072204. PMR
   3.0: Extinction due only to molecules is working fine!!
   3.1: Extinction due only to molecules is working fine, but now is
        true!
   3.2: Extinction due to fixed tau clouds and CIA is now
        supported. 111604. PMR
 */
int version=3;
int revision=2;
int version_rc=-1;


#endif /* _VERSION_TR_H */

