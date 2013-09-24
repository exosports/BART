/*
 * messagep.h   - Header for common messaging routines.
 *
 * Copyright (C) 2003-2007 Patricio Rojo (pato@das.uchile.cl)
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
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *
 */

#ifndef _MESSAGEP_H_
#define _MESSAGEP_H_

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/unistd.h>
#include <errno.h>


/****************************
 *    flags for mperror   *
 ****************************/

#define MSGP_MESSAGE      0x000000
#define MSGP_SYSTEM       0x000001
#define MSGP_USER         0x000002
#define MSGP_WARNING      0x000003

#define MSGP_NOFLAGBITS   0x00000f

#define MSGP_ALLOWCONT    0x000010
#define MSGP_NOPREAMBLE   0x000020
#define MSGP_ALLOC        0x000040
#define MSGP_DBG          0x000080
#define MSGP_NODBG        0x000100


/****************************
 *         Variables        *
 ****************************/

/* keeps tracks of number of errors that where allowed to continue. */
extern int msgp_nowarn;
extern int verblevel;


//#define DEBUG_ERROR
#ifdef  DEBUG_ERROR
#define DBGERR | MSGP_DBG
#else
#define DBGERR
#endif /* DEBUG_ERROR */


/****************************
 *         Macros           *
 ****************************/

#define mperror(flag, ...) \
    mperror_fcn(flag DBGERR, __FILE__, __LINE__, __VA_ARGS__)
#define vmperror(flag, ...) \
    vmperror_fcn(flag DBGERR, __FILE__, __LINE__, __VA_ARGS__)
#define msgpassert(a,...) if(a) mperror(MSGP_SYSTEM,__VA_ARGS__)



#define messagep(thislevel, ...) do{       \
  if(thislevel <= verblevel)               \
    fprintf(stderr,__VA_ARGS__);           \
                                 }while(0)

#define mpallocerror(nmb)                                     \
        mperror(MSGP_SYSTEM,                                  \
	        " %s: Allocation failed for %i allocation\n"  \
	        "units in line %i. Impossible to continue.\n" \
                ,__FILE__,nmb,__LINE__)


#ifdef NODEBUG_MSGP
#define MSGPASSERT(...) ((void)0)
#define MESSAGEP(...) ((void)0)
#else
#define free(x)  do{free(x);x=NULL;}while(0)
#define MSGPASSERT(a,...) if(a) mperror(MSGP_SYSTEM,__VA_ARGS__)
#define MESSAGEP(...) messagep(__VA_ARGS__)
#endif


#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/messagep.c */
extern void messagep_name P_((char *name));
extern void messagep_free P_((void));
extern inline void mpdot P_((int thislevel));
extern int mperror_fcn P_((int flags, const char *file, const long line, const char *str, ...));
extern int vmperror_fcn P_((int flags, const char *file, const long line, const char *str, va_list ap));
extern int fileexistopen P_((char *in, FILE **fp));
extern FILE *verbfileopen P_((char *in, char *desc));
extern void linetoolong P_((int max, char *file, int line));

#undef P_

#endif /* _MESSAGEP_H_ */
