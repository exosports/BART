/*
 * transit.h - Common headers for the Transit program.
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

#ifndef _TRANSIT_H
#define _TRANSIT_H

#include <stdarg.h>
#include <math.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <pu/sampling.h>
#include <pu/profile.h>
#include <pu/iomisc.h>
#include <pu/numerical.h>
#include <pu/xmalloc.h>
#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <alloca.h>
#ifdef _USE_GSL
#include <gsl/gsl_spline.h>
#endif

#define compattliversion 4

#include <flags_tr.h>
#include <constants_tr.h>
#include <types_tr.h>

/*****   Macros   *****/
static __inline__ void
printextprogress(long wi, long wnn){
}

/* 

*/
static __inline__ double
stateeqnford(_Bool mass, /* Abundance by mass (1) or by number (0)  */
             double q,   /* Abundance                               */
             double ma,  /* Mean molecular mass (in AMU)            */
             double mi,  /* Molecular mass of the specie (in AMU)   */
             double p,   /* Pressure                                */
             double t){  /* Temperature                             */
  const double rho = AMU * q * p / KB / t;
  if(mass)
    return rho * ma;
  return rho * mi;
}

#define DEBUG_ERROR
#ifdef  DEBUG_ERROR
#define DBGERR | TERR_DBG
#else
#define DBGERR
#endif /* DEBUG_ERROR */


#define transiterror(flag, ...) \
    transiterror_fcn(flag DBGERR, __FILE__, __LINE__, __VA_ARGS__)
#define vtransiterror(flag, ...) \
    vtransiterror_fcn(flag DBGERR, __FILE__, __LINE__, __VA_ARGS__)


long fw_status;
#define fw(fcn, failurecondition, ...) do{               \
    if((fw_status=fcn(__VA_ARGS__)) failurecondition)    \
      transiterror(TERR_SERIOUS,                         \
                   #fcn "() returned error code %li\n",  \
                   fw_status);                           \
                                       }while(0)


#define transitassert(a, ...) if(a) transiterror(TERR_CRITICAL, __VA_ARGS__)

#define transitprint(thislevel, verblevel, ...) do{                         \
  if(thislevel <= verblevel)  fprintf(stderr, __VA_ARGS__); }while(0)

/* Add flag to transit if present in hint and change hint flag value: */
#define transitacceptflag(transit, hint, flag) do{                          \
        transit |= hint&flag; hint &= ~(flag); }while(0)

#define transitallocerror(nmb)                                              \
        transiterror(TERR_CRITICAL,                                         \
                     "transit:: %s: Allocation failed for %i allocation\n"  \
                     "units in line %i. Impossible to continue.\n",         \
                     __FILE__, nmb, __LINE__)
#define free_null(x) do{free(x);x=NULL;}while(0)

#ifdef NODEBUG_TRANSIT
#define transitDEBUG(...) ((void)0)
#define transitASSERT(...) ((void)0)
#else
#define free(x)  do{free(x);x=NULL;}while(0)
#define transitASSERT(a,...) if(a) transiterror(TERR_CRITICAL,__VA_ARGS__)
#define transitDEBUG(...) transitprint(__VA_ARGS__)
#endif

/* Maximum name length: */
#define maxeisoname 20

#ifndef HAVE__BOOL
#define _Bool short
#endif

extern int transit_nowarn;
extern int verblevel;      /* Verbose level > 10 is only for debuging */
extern int maxline;        /* Initialized in transitstd */
extern int version;
extern int revision;

enum isodo {unclear=0, atmfile, ignore, fixed, factor};

#include <structures_tr.h>

extern const transit_ray_solution slantpath;

extern const eclipse_ray_solution eclipsepath;

/***** Prototypes *****/
#include <proto_transit.h>
#include <proto_readlineinfo.h>
#include <atmosphere/proto_readatm.h>
#include <proto_transitstd.h>
#include <proto_makesample.h>
#include <proto_extinction.h>
#include <proto_idxrefraction.h>
#include <proto_tau.h>
#include <proto_argum.h>
#include <proto_geometry.h>
#include <proto_observable.h>
#include <proto_cia.h>
#include <proto_eclipse.h>

#endif /* _TRANSIT_H */
