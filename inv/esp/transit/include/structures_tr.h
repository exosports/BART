/*
 * transit_structures.h - Structures definition for the transit program
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

#ifndef _TRANSIT_STRUCTURES_H
#define _TRANSIT_STRUCTURES_H

/***** Structures *****/

/* Forward declarations */
struct transit;
struct geometry;

typedef struct {
  const char *name;
  PREC_RES (*tauperb)		/* Optical depth per impact
				   parameter */
       (PREC_RES b,		/* impact parameter */
	PREC_RES *rad,		/* radius array */
	PREC_RES *refr,		/* refractivity index */
	PREC_RES ***ex,		/* extinction[rad][iso][nwn] */
	long nrad,		/* number of radii elements */
	short iso,		/* isotope chosen */
	long wn,		/* wavenumber looked */
	PREC_RES *dt,		/* differential optical depth [rad].
				   Auxiliary array */
	gsl_interp_accel *acc);	/* accelerating pointer. Auxiliary array
				 */
  PREC_RES (*obsperwn)		/* Quantity obtained from
				   integration of optical depth
				*/ 
       (PREC_RES *tau,
	PREC_RES *b,
	long nb,
	struct geometry *star,
	gsl_interp_accel *acc);
} transit_ray_solution;
extern const transit_ray_solution slantpath;


typedef struct {          	/* One item per sampling element */
  PREC_NREC n;			/* number of elements */
  PREC_RES d;			/* Spacing */
  PREC_RES i;			/* initial value */
  PREC_RES f;			/* final value */
  int o;			/* oversampling */
  PREC_RES *v;			/* values of the sampling */
  double fct;			/* Factor to multiply v with to obtain a
				   cgs value */
} prop_samp;


typedef struct {          	/* One item per isotope and
				   miscellaneous conditions, usually
				   radius or temperature */ 
  PREC_ZREC *z;            	/* Partition function [radius or temp] */
  PREC_CS *c;              	/* Cross section [radius or temp] */
  PREC_ATM *d;			/* Environment: Density [radius], not
				   used in lineinfo structure */ 
  PREC_ATM *q;			/* Abundance [radius], not
				   used in lineinfo structure */ 
} prop_isov;


typedef struct {          	/* One item per isotope */
  int d;			/* Database to which they belong */
  char *n;			/* Name */
  PREC_ZREC m;			/* Mass */
} prop_isof;


typedef struct {		/* One item per atmospheric conditions
				   in the atmosphere */
  PREC_ATM *p;			/* Pressure (cgs=dyne/cm2) */
  PREC_ATM *t;			/* Temperature [radius] */
} prop_atm;


typedef struct {          	/* One item per database */
  char *n;			/* Name */
  int i;			/* Number of isotopes */
  int s;			/* Isotopes starting from this index
				   belong to this database */
} prop_db;


typedef struct {          	/* One item per database */
  int t;			/* Number of temperatures */
  PREC_ZREC *T;			/* Temperatures */ 
} prop_dbnoext;


struct lineinfo {		/* Used to keep parameters in
				   readlineinfo() */
  int twii_ver;			/* TWII version */
  int twii_rev;			/* TWII revision */
  prop_samp wavs;		/* wavelength sampling extracted */
  double wi,wf;			/* initial and final wavelength in the
				   database */
  long endinfo;			/* position at the end of the info part
				   of the info file */
  int asciiline;		/* line number in an TWII-ascii file
				   being read, it is zero if a binary
				   file. And the maximum value it gets
				   is the first line of the transition
				   info. */
  prop_isov *isov;		/* Variable isotope information (w/temp)
				   [iso] */
  prop_dbnoext *db;		/* Temperature info from databases [DB]
				   */
};


struct line_transition {	/* One item per transition */
  PREC_LNDATA *wl;		//Wavelength
  PREC_LNDATA *elow;		//Lower energy level
  PREC_LNDATA *gf;		//gf value
  short *isoid;			//Isotope ID (Assumed to be in range)
  double wfct;			//'.wl' multiplied by this factor yields
				//cgs.
  double efct;			//'.elow' multiplied by this factor
				//yields cgs.
};


struct atm_data{		/* Keeps parameters in readatminfo() */
  prop_samp rads;		/* radius sampling */
  prop_isov *isov;		/* variable isotope info [isoext] */
  prop_atm atm;			/* Atmospheric properties */
  int n_niso;			/* Number of new isotopes */
  double *mm;			/* Mean molecular mass [rad] */
  _Bool mass;			/* whether the abundances in 'isov' are
				   mass abundances or not */
  char **n;			/* Name for isotopes in atmfile order */
  PREC_ZREC *m;			/* Mass for isotopes in file order [iso]
				   */
  int *isoeq;			/* Isotope to which each atmosphere
				   datafile column corresponds [iso] */
  enum isodo *isodo;		/* What is required from each isotope,
				   it can be given, ignore, or fixed */
  int n_aiso;			/* Number of isotopes in the atmosphere
				   file */
  char *info;			/* Optional atmosphere file information
				   or label */
  int begline;			/* line of beginning of radius dependent
				   info */
  long begpos;			/* position of beginning of radius
				   dependent info*/
};


struct extinction{
  PREC_RES ***e;		/* Extinction value [rad][iso:][wav]*/
  float maxratio;		/* Maximum Doppler width ratio between
				   current and last calculated profile.
				   If the value is greater than this,
				   then recalculate */
  int vf;			/* Number of fine-bins of the Voigt
				   function */
  float ta;			/* number of alphas that have to be
				   contained in the profile */
};


struct idxref {
  PREC_RES *n;			/* Index of refraction [rad] */
};


struct onept {
  double p,t;			/* pressure, temperature */
  double *q;			/* abundances for isotopes */
  int nq;			/* number of given abundances */
  _Bool one;			/* One point is required? */
  char **n;			/* Names of extra isotopes */
  PREC_ZREC *m;			/* Mass of extra isotopes */
  int nm;			/* number of given name and masses */
  int ne;			/* Number of extra isotopes */
};


struct optdepth {
  PREC_RES **t;			/* Optical depth [wn][ip] */
  short iso;			/* Isotope from which to calculate the
				   optical depth */
  long *first;			/* Index of the lowest impact parameter
				   value, lower than this the optical
				   depth is greater than '.toomuch'. It
				   is naturally assumed that optical
				   depth increases inward the
				   planet. [wn] */
  double toomuch;		/* Optical depth values greater than
				   this won't be calculated: the
				   extinction is assumed to be zero. */
};


struct geometry {
  float smaxis;			/* Semimajor axis */
  double smaxisfct;		/* 'smaxis' times this gives cgs
				   units. */
  double time;			/* this value is 0 when in the middle of
				   the eclipse */
  double timefct;		/* 'time' times this gives cgs units */
  float incl;			/* inclination of the planetary orbit
				   with respect to the observer, 90
				   degrees is edge on */
  float inclfct;		/* Units to convert inclination to
				   radians */
  double ecc;			/* eccentricty */
  double eccfct;		/* eccentricity's units */
  double lnode;			/* longitud of the ascending node */
  double lnodefct;		/* longitud of the ascending node units */
  double aper;			/* argument of the pericenter */
  double aperfct;		/* argument of the pericenter units */


  double starmass;		/* Mass of the star */
  double starmassfct;		/* 'starmass' times this gives cgs
				   units. */

  double starrad;		/* Star's radius */
  double starradfct;		/* 'starrad' times this gives cgs
				   units. */

  double x0,y0;			/* coordinates of the center of the
				   planet with respect to the
				   star. 'fct' to convert to cgs is
				   found in rads.fct. These fields are
				   not hinted. */
};


struct transithint {		/* Structure with user hinted data that
				   should go to the 'struct transit'
				   upon approval */
  char *f_atm,*f_line,*f_out;	/* Filenames */
  prop_samp rads,wavs,wns;	/* Sampling properties of
				   radius, wavelength and
				   wavenumber */
  prop_samp ips;		/* Impact parameter sampling, at what
				   radius sampling does the user wants
				   ray optical depth to be calculated */
  float allowrq;		/* How much less than one is accepted,
				   and no warning is issued if
				   abundances don't ad up to that */
  PREC_RES t;			/* Telescope resolution */
  PREC_RES m;			/* Amount of nanometers not trusted at
				   the boundaries, also how much out of
				   requested range it have to look for
				   transitions */
  PREC_RES wnm;			/* Same as above, but for wavenumbers */
  float maxratio_doppler;	/* Maximum doppler width deformation
				   ratio before recalculate profile */
  float timesalpha;		/* Number of alphas that have to be
				   contained in a calculated profile,
				   one side only */
  int voigtfine;		/* Fine-binning for Voigt function in
				   kapwl(), if accepted it goes to
				   tr.ds.op.vf */
  int verbnoise;		/* noisiest verbose level in a non
				   debugging run */ 
  _Bool mass;			/* whether the abundances read by getatm
				   are by mass or number */
  long fl;			/* flags */
  long na;			/* flags of non-accepted or just changed
				   hints */
  struct onept onept;		/* Parameters for onept atmosphere */
  double toomuch;		/* Optical depth values greater than
				   this won't be calculated: the
				   extinction is assumed to be zero. */
  short tauiso;			/* Whether user want to calculate
				   optical depth for all or some
				   isotopes, TRU_EXTPERISO has to be
				   on. */
  char *solname;		/* Name of the type of solution */
  struct geometry sg;		/* System geometry */
};


struct transit {		/* Main data structure */
  char *f_atm,*f_line,*f_out;	/* Filenames */
  FILE *fp_atm,*fp_out,*fp_line;/* Filepointers */
  float allowrq;		/* How much less than one is accepted,
				   so that no warning is issued if
				   abundances don't ad up to that */
  PREC_RES telres;		/* Telescope resolution */
  PREC_RES m;			/* Amount of nanometers not trusted at
				   the boundaries, also how much out of
				   requested range it have to look for
				   transitions */
  PREC_RES wnmi;		/* Amount of cm-1 not trusted at the
				   beginning */
  PREC_RES wnmf;		/* Amount of cm-1 not trusted at the end
				   */
  prop_samp rads,wavs,wns;	/* Sampling properties of radius,
				   wavelength and wavenumber */
  prop_samp ips;		/* Impact parameter sampling, at what
				   radius sampling does the user wants
				   ray optical depth to be calculated */
  prop_atm atm;			/* Sampled atmospheric data. Height in
				   kilometers. */
  prop_isof *isof;		/* Fixed isotope information
				   [isoextended] */
  prop_isov *isov;		/* Variable isotope information
				   [isoextended] */
  prop_db *db;			/* Database's info [DB] */
  int n_db,n_i,n_e;		/* Number of databases, of regular
				   isotopes, of extended isotopes, and
				   of lines in database */
  PREC_NREC n_l;		/* Number of lines in database */
  enum isodo *isodo;		/* What to do with every isotope */

  long fl;			/* flags */
  long pi;			/* progress indicator */

  struct line_transition lt;	/* line transition */

  transit_ray_solution *sol;	/* Solution type */
  PREC_RES *outpret;		/* Output dependent on wavelength only
				   as it travels to Earth before
				   telescope */
  PREC_RES *out;		/* Output as seen after interaction with
				   telescope */

  struct {			/* data structures pointers, this is
				   data that is not required for the
				   final computation */
    struct iso_noext *in;
    struct atm_data *at;
    struct transithint *th;
    struct lineinfo *li;
    struct extinction *ex;
    struct optdepth *tau;
    struct idxref *ir;
    struct geometry *sg;
  }ds;
};


#endif /* _TRANSIT_STRUCTURES_H */
