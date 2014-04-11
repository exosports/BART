/*
 * constants_tr.h - constants for the transit program.
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

#ifndef _CONSTANTS_TR_H
#define _CONSTANTS_TR_H


/*  Constants    */
/* Units in cgs: */
#if defined(AMU)       || defined(EC)          || defined(LS)          ||  \
    defined(ME)        || defined(KB)          || defined(H)           ||  \
    defined(PI)        || defined (SIGCTE)     || defined(EXPCTE)      ||  \
    defined(WNU_O_WLU) || defined(AU)          || defined(SUNMASS)     ||  \
    defined(SUNRADIUS) || defined(HOUR)        || defined(ONEOSQRT2PI) ||  \
    defined(SQRTLN2)
#error Some of the preprocessor constants were already defined!
#endif

#define RHOSTP (1.29e-3)         /* Density at standard temp and press       */
#define PI (3.141592653589793)   /* PI                                       */
#define DEGREES (PI/180.0)       /* to make degrees to radian                */
#define GGRAV (6.673e-8)         /* Gravitational constant (erg * cm / g^2)  */
#define HOUR (3600.0)            /* 1 hour (s)                               */
#define AU (14959786896040.492)  /* 1 Astronomical unit (cm)                 */
#define ANGSTROM (1e-8)          /* 1 Angstrom (cm)                          */
#define SUNMASS (1.9891e33)      /* Solar mass (g)                           */
#define SUNRADIUS (6.95508e10)   /* solar volumetric mean radius (cm)        */
#define AMU (1.66053886e-24)     /* Atomic Mass unit (g)                     */
#define LO (2.686763e19)         /* Loschmidt constant (cm-3)                */
#define EC (4.8032068e-10)       /* electronic charge (statcoulomb)          */
#define LS (2.99792458e10)       /* Light Speed (cm / s)                     */
#define ME (9.1093897e-28)       /* Electron mass (g)                        */
#define KB (1.380658e-16)        /* Boltzmann constant (erg / K)             */
#define H (6.6260755e-27)        /* Planck's constant (erg * s)              */
#define HC (H*LS)                /* for lower energy conversion (erg * cm)   */
#define SIGCTE (PI*EC*EC/LS/LS/ME/AMU) /* Cross-sec constant (cm / g)        */
#define EXPCTE (H*LS/KB)         /* Exponent constant (cm * K)               */

#define ONEOSQRT2PI (0.3989422804)         /* 1.0/sqrt(2pi)                  */
#define SQRTLN2  (0.83255461115769775635)  /* sqrt(ln(2))                    */

#define MAXNAMELEN 20

#ifdef __LITTLE_ENDIAN
/* {0xff-'t',0xff-'r',0xff-'s',0xff-'f'} */
#define __TR_SAVEFILE_MN__      "\xb5\xb7\xb6\xbd"
#else
#define __TR_SAVEFILE_MN__      "\xbd\xb6\xb7\xb5"
#endif /* Little or big endian */

#endif /* _CONSTANTS_TR_H */
