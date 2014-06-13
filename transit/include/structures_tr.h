/*
 * structures_tr.h - Structures definition for the transit program
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

#ifndef _TRANSIT_STRUCTURES_H
#define _TRANSIT_STRUCTURES_H

/*  Structures  */

/* Forward declarations */
struct transit;
struct geometry;

/* Type of ray solution, eclipse or transit */
typedef enum {transit, eclipse}  RaySol;

/* Structure definitions */
typedef struct {    /* Sampling struct                 */
  PREC_NREC n;      /* Number of elements              */
  PREC_RES d;       /* Spacing                         */
  PREC_RES i;       /* Initial value                   */
  PREC_RES f;       /* Final value                     */
  int o;            /* Oversampling                    */
  PREC_RES *v;      /* Values of the sampling          */
  double fct;       /* v units factor to cgs           */
} prop_samp;


typedef struct {    /* Isotope's variable (per layer) information: */
  unsigned int n;   /* Arrays' length                              */
  double *z;        /* Partition function [radius or temp]         */
  double *c;        /* Cross section      [radius or temp]         */
} prop_isov;


typedef struct {    /* Isotope's fixed information:  */
  int d;            /* Database to which they belong */
  char *n;          /* Isotope name                  */
  PREC_ZREC m;      /* Isotope mass                  */
} prop_isof;


typedef struct{    /* Molecule's information: */
  int n;           /* Number of elements      */
  PREC_ATM *d;     /* Density   [n]           */
  PREC_ATM *q;     /* Abundance [n]           */
} prop_mol;


typedef struct {    /* Atmospheric conditions:          */
  double *mm;       /* Mean molecular mass [rad]        */
  PREC_ATM *p;      /* Pressure    [rad]                */
  PREC_ATM *t;      /* Temperature [rad]                */
  PREC_ATM pfct;    /* p units factor to cgs (dyne/cm2) */
  PREC_ATM tfct;    /* t units factor to cgs (Kelvin)   */
} prop_atm;


typedef struct {    /* Database properties:             */
  char *n;          /* Database name                    */
  unsigned int i;   /* Number of isotopes               */
  int s;            /* Cumulative first isotope's index */
} prop_db;


typedef struct {    /* One item per database  */
  unsigned int t;   /* Number of temperatures */
  double *T;        /* Temperatures           */ 
} prop_dbnoext;


typedef struct {         /* Ray solution properties:              */
  const char *name;      /* Ray solution name                     */
  const char *file;      /* Ray solution filename (FINDME)        */
  const short monoip;    /* Request equispaced inpact parameter?  */
  PREC_RES (*tauperb)    /* Optical depth per impact parameter:   */
       (PREC_RES b,      /*  Impact parameter                     */
        PREC_RES *rad,   /*  Radius array                         */
        PREC_RES *refr,  /*  Refractivity index                   */
        PREC_RES *ex,    /*  Extinction[rad]                      */
        long nrad,       /*  Number of radii elements             */
        int exprlevel);  /*  FINDME                               */
  PREC_RES (*obsperwn)         /* Integrated optical depth:          */
        (PREC_RES *tau,        /*  Optical depth                     */
        long last,             /*  index where tau exceeded toomuch  */
        PREC_RES toomuch,      /*  Cutoff optical depth to calculate */ 
        prop_samp *ip,         /*  Impact parameter                  */
        struct geometry *star, /*  Geometry structure                */
        int exprlevel);        /*  Modulation level                  */
  const int nobs;        /* Number of levels of details as it can
                            be handled by obsperwn                */
} transit_ray_solution;



typedef struct {              /* Ray solution properties:                 */
  const char *name;           /* Ray solution name                        */
  const char *file;           /* Ray solution filename (FINDME)           */
  PREC_RES (*tauEclipse)      /* Optical depth for eclipse per wavenumber */
       (PREC_RES *rad,        /* Radius array                             */
        PREC_RES *ex,         /* Extinction[rad]                          */
        PREC_RES angle,       /* Ray grid                                 */
        long nrad);           /* Number of radii elements                 */
  PREC_RES (*eclIntenWn)      /* Integrated optical depth:                */
        (struct transit *tr,  /* Main structure                           */
        PREC_RES *tau,        /* Optical depth                            */
        PREC_RES w,           /* Current wavenumber value                 */
        long last,            /* Index where tau exceeded toomuch         */
        PREC_RES toomuch,     /* Cutoff optical depth to calculate        */
        prop_samp *rad);      /* Impact parameter                         */
} eclipse_ray_solution;



struct atm_isoprop{    /* Proportional-abundance isotopic parameters: */
  double f;            /* Fractional abundance                        */
  double m;            /* Isotope mass                                */
  int eq;              /* Isotope index from transit.ds.isotopes      */
  char n[maxeisoname]; /* Isotope name                                */
  char t[maxeisoname]; /* Molecule name                               */
};


struct line_transition{  /* Line transition parameters:         */
  PREC_LNDATA *wl;       /* Wavelength                          */
  PREC_LNDATA *elow;     /* Lower energy level                  */
  PREC_LNDATA *gf;       /* gf value                            */
  short *isoid;          /* Isotope ID (Assumed to be in range) */
  double wfct;           /* wl units factor to cgs              */
  double efct;           /* elow units factor to cgs            */
};


struct lineinfo{             /* Line information parameters:    */
  struct line_transition lt; /* Line transitions                            */
  unsigned short tli_ver;    /* TLI version                                 */
  unsigned short lr_ver;     /* lineread version                            */
  unsigned short lr_rev;     /* lineread revision                           */
  prop_samp wavs;            /* Wavelength sampling                         */
  double wi, wf;             /* Initial and final wavelength in database    */
  long endinfo;              /* Position at the end of the info part
                                of the info file                            */
  int asciiline;             /* TLI line of first transition.
                                Zero if a binary file.                      */
  int ni;                    /* Number of isotopes                          */
  int ndb;                   /* Number of databases                         */
  prop_isov *isov;           /* Variable isotope information (w/temp) [iso] */
  prop_dbnoext *db;          /* Temperature info from databases [DB]        */
  PREC_NREC n_l;             /* Number of lines in database                 */
};


struct atm_data{     /* Atmospheric file parameters:                    */
  int n_aiso;        /* Number of molecules in atmosphere file          */
  prop_samp rads;    /* Radius sampling                                 */
  prop_atm atm;      /* Atmospheric properties                          */
  prop_mol *molec;   /* Molecular information [n_aiso]                  */
  double *mm;        /* Mean molecular mass [rad]                       */
  char *info;        /* Optional atmosphere file information or label   */
  _Bool mass;        /* Abundances in isov by mass (1) of by number (0) */
  int begline;       /* Line of first radius dependent info             */
  long begpos;       /* Position of first radius dependent info         */
};


struct extinction{
  PREC_RES ***e;     /* Extinction value [iso:][rad][wav]                 */
  float maxratio;    /* Maximum Doppler width ratio between current and
                        last calculated profile.  If the value is greater
                        than this, then recalculate                       */
  int vf;            /* Number of fine-bins of the Voigt function         */
  float ta;          /* Number of alphas that have to be contained in
                        the profile                                       */
  _Bool periso;      /* Extinction per isotope                            */
  _Bool *computed;   /* Whether the extinction at the given radius was
                        computed [rad]                                    */
  double minelow;    /* Only use transitions with this minimum low
                        energy (in cm-1)                                  */
};


struct idxref{
  PREC_RES *n;       /* Index of refraction [rad] */
};


struct onept{
  double p, t;   /* Pressure and temperature        */
  int ne;        /* Number of extra isotopes        */
  _Bool one;     /* One PT required?                */
  double *q;     /* Isotopes abundances             */
  int nq;        /* Number of input abundances      */
  char **n;      /* Names of extra isotopes         */
  PREC_ZREC *m;  /* Mass of extra isotopes          */
  int nm;        /* Number of input mass-name pairs */
};


#if 0
struct savefiles {
  char *ext;   /* saves extinction */
  char *tau;   /* after tau() savefile */
  char *modulation;  /* after modulation() savefile */
};
#endif


struct optdepth{
  PREC_RES **t;     /* Optical depth [wn][ip]                             */
  long *last;       /* Index of the lowest impact parameter value, lower
                       than this the optical depth is greater than
                       '.toomuch'.  It is naturally assumed that optical 
                       depth increases inward the planet. [wn]            */
  double toomuch;   /* Optical depth values greater than this won't be
                       calculated: the extinction is assumed to be zero.  */
};

struct grid{
  PREC_RES **a;      /* Intensity grid, 2D, [an][wnn]                      */
};



struct geometry{
  float smaxis;       /* Semimajor axis                                    */
  double smaxisfct;   /* 'smaxis' times this gives cgs units.              */
  double time;        /* this value is 0 when in the middle of the eclipse */
  double timefct;     /* 'time' times this gives cgs units                 */
  float incl;         /* inclination of the planetary orbit with respect
                         to the observer, 90 degrees is edge on            */
  float inclfct;      /* Units to convert inclination to radians           */
  double ecc;         /* eccentricty                                       */
  double eccfct;      /* eccentricity's units                              */
  double lnode;       /* longitud of the ascending node                    */
  double lnodefct;    /* longitud of the ascending node units              */
  double aper;        /* argument of the pericenter                        */
  double aperfct;     /* argument of the pericenter units                  */

  double starmass;    /* Mass of the star                                  */
  double starmassfct; /* 'starmass' times this gives cgs units.            */

  double starrad;     /* Star's radius                                     */
  double starradfct;  /* 'starrad' times this gives cgs units.             */

  double x, y;        /* Coordinates of the center of the planet with
                         respect to the star. 'fct' to convert to cgs is
                         found in rads.fct. These fields are not hinted.   */

  _Bool transpplanet; /* If true, set maximum optical depth to toomuch     */
};


struct isotopes{
  prop_isof *isof;    /* Fixed isotope information      [n_i] */
  prop_isov *isov;    /* Variable isotope information   [n_i] */
  double *isoratio;   /* Isotopic abundance ratio       [n_i] */
  int *imol;          /* Molecule index for this isotope[n_i] */
  prop_db *db;        /* Database's info [n_db]               */
  int n_db,           /* Number of databases                  */
      n_i;            /* Number of isotopes                   */
};

struct molecules{
  int nmol;        /* Number of molecules  */
  prop_mol *molec; /* Molecular properties */
  char **name;     /* Molecules' names     */
  PREC_ZREC *mass; /* Molecules' masses    */
  double *radius;  /* Molecules' radii     */
};


struct outputray{
  PREC_RES *o;     /* Output as seen before interaction with telescope */
};

struct extcloud{
  double maxe;     /* Maximum opacity in [cm-1]                          */
  double rini;     /* Radius at which clouds start                       */
  double rfin;     /* Radius at which clouds has it maximum thickness
                      'maxe'. rfin < rini                                */
  double rfct;     /* Factor that will make the two radius values above
                      into cgs                                           */
};

struct extscat{
  double prm;
};

struct saves{
  char *ext;
};

/* Struct to store requested ext, tau, or cia detailed information:  */
struct detailfld{
  int n;         /* Number of requested wavenumber samples */
  PREC_RES *ref; /* Array of wavenumbers requested         */
  char file[80]; /* Output filename                        */
  char name[30]; /* Name of field                          */
};


struct detailout{
  struct detailfld ext, tau, cia;
};


struct cia{
  PREC_CIA **e;   /* Extinction from all CIA sources [wn][tmp] */
  char **file;
  int n;
};

/* Structure with user hinted data that should go to the 'struct
   transit' upon approval                                         */
struct transithint{  
  char *f_atm,          /* Atmosphere filename      */
       *f_line,         /* TLI filename             */
       *f_out,          /* Output (main) filename   */
       *f_toomuch,      /* Output toomuch filename  */
       *f_outsample;    /* Output sample filename   */
  PREC_NREC ot;         /* Radius index at which to print output from tau    */
  prop_samp rads, wavs, wns; /* Sampling properties of radius, wavelength
                                and wavenumber                               */
  RaySol  path;         /* Eclipse or transit ray solution.                  */
  long int ann;         /* Number of angles                                  */
  PREC_RES angles[10];  /* Angles                                            */
  prop_samp ips;        /* Impact parameter sampling, at what radius
                           sampling does the user wants ray optical depth to
                           be calculated                                     */
  float allowrq;        /* How much less than one is accepted, and no warning
                           is issued if abundances don't ad up to that       */
  PREC_RES margin;      /* Amount not trusted at the boundaries (uses
                           .wavs.fct to convert to cgs), also how much out of
                           requested range it have to look for transitions   */
  PREC_RES wnm;         /* Same as above, but for wavenumbers                */
  float maxratio_doppler; /* Maximum doppler width deformation ratio before
                             recalculate profile                             */
  float timesalpha;     /* Number of alphas that have to be contained in a
                           calculated profile, one side only                 */
  int voigtfine;        /* Fine-binning for Voigt function in kapwl(), if
                           accepted it goes to tr.ds.op.vf                   */
  int verbnoise;        /* Noisiest verbose level in a non debugging run     */ 
  _Bool mass;           /* Whether the abundances read by getatm are by
                           mass or number                                    */
  long fl;              /* flags                                             */
  _Bool userefraction;  /* Whether to use variable refraction                */

  double toomuch;       /* Optical depth values greater than this won't be
                           calculated: the extinction is assumed to be zero  */
  short tauiso;         /* Whether user want to calculate optical depth for
                           all or some isotopes, TRU_EXTPERISO has to be on. */
  double blowex;        /* Blow extinction by this amount before computing
                           tau, this option has no physical meaning, but
                           mostly debugging                                  */
  int taulevel;         /* Tau integration level of precision                */
  int modlevel;         /* Modulation integration level of precision         */
  char *solname;        /* Name of the type of solution                      */
  struct geometry sg;   /* System geometry                                   */
  struct onept onept;   /* Parameters for onept atmosphere                   */
  struct saves save;    /* Saves indicator of program stats                  */

  struct extcloud cl;
  struct detailout det;

  double minelow;       /* Only use transitions with this minimum low
                           energy (in cm-1)                                  */
  char **ciafile;
  int ncia;

};

/* Main data structure */
struct transit{  
  char *f_atm,       /* Atmosphere filename      */
       *f_line,      /* TLI filename             */
       *f_out,       /* Output (main) filename   */
       *f_toomuch,   /* Output toomuch filename  */
       *f_outsample; /* Output sample filename   */
  PREC_NREC ot;      /* Radius index at which to print output from tau       */

  FILE *fp_atm, *fp_out, *fp_line; /* Pointers to files                      */
  float allowrq;    /* How much less than one is accepted, so that no warning
                       is issued if abundances don't ad up to that           */
  PREC_RES telres;  /* Telescope resolution                                  */
  PREC_RES margin;  /* Wavelength amount not trusted at the boundaries in
                       cgs units, also how much out of requested range it
                       have to look for transitions                          */
  PREC_RES wnmi;    /* Amount of cm-1 not trusted at the beginning           */
  PREC_RES wnmf;    /* Amount of cm-1 not trusted at the end                 */
  long int angleIndex; /* Index of the current angle                         */
  PREC_RES *Flux;   /* Flux for eclipse                                      */
  prop_samp rads, wavs, wns; /* Sampling properties of radius, wavelength
                                and wavenumber                               */
  prop_samp ips;    /* Impact parameter sampling, at what radius sampling does
                       the user wants ray optical depth to be calculated     */
  prop_atm atm;     /* Sampled atmospheric data                              */
  short tauiso;     /* Isotope from which to calculate the optical depth     */
  double blowex;    /* Blow extinction by this amount before computing tau,
                       this option has no physical meaning, but mostly
                       debugging                                             */
  int taulevel;     /* Tau integration level of precision                    */
  int modlevel;     /* Modulation integration level of precision             */

  long fl;          /* flags                                                 */
  long pi;          /* progress indicator                                    */

  transit_ray_solution *sol; /* Transit solution type                        */
  eclipse_ray_solution *ecl; /* Eclipse solution type                        */
  PREC_RES *outpret; /* Output dependent on wavelength only as it travels
                        to Earth before telescope                            */

  struct saves save; /* Saves indicator of program stats                     */

  struct {          /* Data structures pointers, this is data that is not
                       required for the final computation                    */
    struct transithint *th;
    struct lineinfo    *li;
    struct atm_data    *at;
    struct extinction  *ex;
    struct grid        *intens;
    struct optdepth    *tau;
    struct idxref      *ir;
    struct geometry    *sg;
#if 0
    struct savefiles   *sf;
#endif
    struct isotopes    *iso;
    struct molecules   *mol;
    struct outputray   *out;
    struct extcloud    *cl;
    struct extscat     *sc;
    struct detailout   *det;
    struct cia         *cia;
  }ds;
};

#endif /* _TRANSIT_STRUCTURES_H */
