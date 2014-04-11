/*
 * flags_tr.h - Various flags for the transit program
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

#ifndef _FLAGS_TR_H
#define _FLAGS_TR_H

/*****   Flags   *****/

/* transithint parameters not accepted or changed: */
/* File flags:                  */
#define TRH_FL          0x00000001 /* Line info file        */
#define TRH_FA          0x00000002 /* Atmospheric info file */
#define TRH_FO          0x00000004 /* Output file           */
/* Margin and resolution flags: */
#define TRH_WNM         0x00000008 /* Wavenumber margin    */
#define TRH_TR          0x00000010 /* Telescope resolution */
#define TRH_WM          0x00000020 /* Wavelength margin    */
/* Line broadening flags:      */
#define TRH_VF          0x00000040 /* Voigt fine-binning */
#define TRH_TA          0x00000080 /* Times of alpha     */
#define TRH_DR          0x00000100 /* Max doppler ratio  */

#define TRH_MASS        0x00000200 /* mass abundance? */
#define TRH_TOOMUCH     0x00000400 /* Limit optical depth, above this is
                                      just set to this value. */
#define TRH_TAUISO      0x00000800 /* Optical depth is going to be for
                                      only this isotope. To show all of
                                      them, is it either -1, or 0 if
                                      TRU_EXTINPERISO disabled */
#define TRH_ST          0x00001000 /* Solution type */

#define TRH_WAVO        0x01000000 /* Wavelength oversampling       */
#define TRH_WNO         0x02000000 /* Wavenumber oversampling       */

/* Sampling flags: */
#define TRH_RAD         0x10000000 /* Radius sample specified       */
#define TRH_WAV         0x20000000 /* Wavelength sample specified   */
#define TRH_WN          0x40000000 /* Wavenumber sample specified   */
#define TRH_IPRM        0x80000000 /* Impact param sample specified */
#define TRH_NAME(n) (n==TRH_RAD?"radius":                          \
                    (n==TRH_WAV?"wavelength":                      \
                    (n==TRH_WN?"wavenumber":                       \
                    (n==TRH_IPRM?"impact parameter":               \
                                 "unknown(a.k.a bad 'fl' value)"))))

/* Internal flags: */
#define TRF_NOOVERSAMP  0x00000001 /* No oversampling in printsample files */ 
#define TRF_NOVALUE     0x00000002 /* Do not print each of the values */

/* Flags for interpolation mode: */
#define TRU_SAMPLIN     0x00000001 /* Linear interpolation */
#define TRU_SAMPSPL     0x00000002 /* Spline interpolation */
#define TRU_SAMPBITS    0x00000003 /* FINDME */

#define TRU_CNVFIX      0x00000000 /* Fix width: accepts one width(float)    */
#define TRU_CNVLINEAR   0x00000010 /* Linear change of width, accepts 
                                      initial(float) and ending width(float) */
#define TRU_CNVGIVEN    0x00000020 /* Resolution width is given, accepts 
                                      array(float *) and number(int) of 
                                      elements on it                         */
#define TRU_CNVFFT      0x00000030 /* FFT of a fixed width convolver         */

#define TRU_CNVGAUSSIAN 0x00000100 /* Gaussian convolution   */
#define TRU_CNVBOX      0x00000200 /* Square box convolution */

#define TRU_CNVMODBITS  0x000000f0 /* Bits that define the mode   used */ 
#define TRU_CNVMTHBITS  0x00000f00 /* Bits that define the method used */

/* User defined behavior flags: */
#define TRU_ATMNODEF    0x00000000 /* Don't use atmospheric defaults, fail
                                      if that condition is reached           */
#define TRU_ATMHARDC1P  0x00001000 /* Use hard coded values                  */
#define TRU_ATMASK1P    0x00002000 /* Obtain one-point atm. pars from stdin  */
#define TRU_ATM1PBITS   0x0000f000 /* One-point atmospheric behavior         */

/* Bits that define atmospheric parameters:   */
#define TRU_ATMBITS (TRU_ATM1PBITS|TRU_SAMPBITS)

#define TRU_EXTINPERISO 0x00100000 /* There won't be a calculation of
                                      extinction in a per isotope array,
                                      all of them should be combined      */
#define TRU_EXTBITS     0x00f00000

#define TRU_OUTTAU      0x01000000 /* Print out optical depth */
#define TRU_TAUBITS     0x0f000000 /* FINDME */

/* Progress indicator flags:  */
#define TRPI_READINFO     0x000001 /* readinfofile()   completed */
#define TRPI_READDATA     0x000002 /* readdatarng()    completed */
#define TRPI_CHKRNG       0x000004 /* checkrange()     completed */
#define TRPI_GETATM       0x000008 /* getatm()         completed */
#define TRPI_MAKERAD      0x000010 /* makeradsample()  completed */
#define TRPI_MAKEWAV      0x000020 /* makewavsample()  completed */
#define TRPI_MAKEWN       0x000040 /* makewnsample()   completed */
#define TRPI_MAKEIP       0x000080 /* makeipsample()   completed */
#define TRPI_IDXREFRAC    0x000100 /* idxrefrac()      completed */
#define TRPI_EXTWN        0x000200 /* extwn()          completed */
#define TRPI_TAU          0x000400 /* tau()            completed */
#define TRPI_GEOMETRY     0x000800 /* setgeom()        completed */
#define TRPI_GEOMETRYHINT 0x001000 /* setgeomhint()    completed */
#define TRPI_MODULATION   0x002000 /* modulation()     completed */
#define TRPI_CIA          0x004000 /* interpolatecia() completed */

/* flags for transiterror:  */
#define TERR_MESSAGE      0x000000
#define TERR_CRITICAL     0x000001
#define TERR_SERIOUS      0x000002
#define TERR_WARNING      0x000003
#define TERR_NOFLAGBITS   0x00000f
#define TERR_ALLOWCONT    0x000010
#define TERR_NOPREAMBLE   0x000020
#define TERR_ALLOC        0x000040
#define TERR_DBG          0x000080

#endif /* _FLAGS_TR_H */
