#ifndef _TRANSIT_H
#define _TRANSIT_H

#define TRANSIT

#include <P/Putil.h>
#include <stdarg.h>
#include <math.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <util/sampling.h>
#include <util/profile.h>
#include <util/iomisc.h>
#include <strings.h>


/*****   Flags   *****/
/* Flags for hintable parameters. Its precense indicates that the
   parameter was changed or not accepted */
#define TRH_FL          0x00000001 /* Line info file */
#define TRH_FA          0x00000002 /* Atmospheric info file */
#define TRH_FO          0x00000004 /* Output file */

#define TRH_WNM         0x00000008 /* Wavenumber margin */
#define TRH_TR          0x00000010 /* Telescope resolution */
#define TRH_WM          0x00000020 /* Wavelength margin */

#define TRH_VF          0x00000040 /* Voigt fine-binning */
#define TRH_TA          0x00000080 /* times of alpha */
#define TRH_DR          0x00000100 /* max doppler ratio */

#define TRH_MASS        0x00000200 /* mass abundance? */

/* Generic flags for sampling to be shifted as below to indicate radius,
   wavenumber and wavelength spacing, respectively */
/* According to the current value of TRH\_RADSFT, there are 15 bits used
   for sampling flags, that leaves 17 bits for other purposes (i.e. the
   lower half plus one */
#define TRH_RADSFT         (32-15) /* Number of bits to shift for radius
				      sample flags */
#define TRH_WAVSFT         (32-10) /* Number of bits to shift for
				      wavelength sample flags*/
#define TRH_WNSFT           (32-5) /* Number of bits to shift for
				      wavelength sample flags*/
#define TRH_SI          0x00000001 /* Initial value */
#define TRH_SF          0x00000002 /* Final value */
#define TRH_SD          0x00000004 /* value spacing */
#define TRH_SN          0x00000008 /* Value data number. Implies a
				      given array */
#define TRH_SO          0x00000010 /* Value oversampling */
#define TRH_SBITS       0x0000001f /* Mask for all sample flags */
#define TRH_SFTNAME(n) (n==TRH_RADSFT?"radius":         \
			 (n==TRH_WAVSFT?"wavelength":    \
			  (n==TRH_WNSFT?"wavenumber":    \
			   "unknown(a.k.a bad 'bitsshift' value)")))
			 

#define TRH_RI  (TRH_SI<<TRH_RADSFT) /* Initial radius */
#define TRH_RF  (TRH_SF<<TRH_RADSFT) /* Final radius */
#define TRH_RD  (TRH_SD<<TRH_RADSFT) /* Radius spacing */
#define TRH_RN  (TRH_SN<<TRH_RADSFT) /* Radius data number. Implies a
					given array */
#define TRH_RO  (TRH_SO<<TRH_RADSFT) /* Radius oversampling */

#define TRH_WI  (TRH_SI<<TRH_WAVSFT) /* Initial wavelength */
#define TRH_WF  (TRH_SF<<TRH_WAVSFT) /* Final wavelength */
#define TRH_WD  (TRH_SD<<TRH_WAVSFT) /* Wavelength spacing */
#define TRH_WN  (TRH_SN<<TRH_WAVSFT) /* Wavelength data number. Implies a
					given array */
#define TRH_WO  (TRH_SO<<TRH_WAVSFT) /* Wavelength oversampling */

#define TRH_WNI  (TRH_SI<<TRH_WNSFT) /* Initial wavenumber */
#define TRH_WNF  (TRH_SF<<TRH_WNSFT) /* Final wavenumber */
#define TRH_WND  (TRH_SD<<TRH_WNSFT) /* Wavenumber spacing */
#define TRH_WNN  (TRH_SN<<TRH_WNSFT) /* Wavenumber data number. Implies a
					given array */
#define TRH_WNO  (TRH_SO<<TRH_WNSFT) /* Wavenumber oversampling */



/* Flags for mode of telresconv */
#define TRU_CNVFIX      0x00000000 /* Fix width: accepts one
				      width(float) */ 
#define TRU_CNVLINEAR   0x00000001 /* Linear change of width,
				      accepts initial(float) and
				      ending width(float) */
#define TRU_CNVGIVEN    0x00000002 /* Resolution width is given,
				      accepts array(float *) and
				      number(int) of elements on it */
#define TRU_CNVFFT      0x00000003 /* FFT of a fixed width convolver */

#define TRU_CNVGAUSSIAN 0x00000010 /* Gaussian convolution */
#define TRU_CNVBOX      0x00000020 /* Square box convolution */

#define TRU_CNVMODBITS  0x0000000f /* Bits that define the mode used */ 
#define TRU_CNVMTHBITS  0x000000f0 /* Bits that define the method used
				      */

/* flags to indicate various user defined behaviors */
#define TRU_ATMNODEF    0x00000000 /* Don't use atmospheric defaults,
				      fail if that condition is
				      reached */
#define TRU_ATMHARDC1P  0x00000100 /* Use hard coded values */
#define TRU_ATMGIVEN1P  0x00000200 /* Atmosphere's one point is given
					in hint */
#define TRU_ATMASK1P    0x00000300 /* Atmosphere's one point is
					obtained from stdin */
#define TRU_ATM1PBITS   0x00000f00 /* Bits that define one point
				      atmospheric behavior */
#define TRU_ATMBITS     0x0000ff00 /* Bits that define atmospheric
				      parameters */


#define TRU_SAMPBITS    0x00007000
#define TRU_SAMPLIN	0x00001000
#define TRU_SAMPSPL     0x00002000

#define TRU_EXTBITS     0x00038000
#define TRU_EXTINPERISO 0x00008000 /* There won't be a calculation of
				      extinction in a per isotope array,
				      all of them should be combined */


/* Progress indicator flags */
#define TRPI_READINFO   0x000001 /* readinfofile() completed */
#define TRPI_READDATA   0x000002 /* readdatarng() completed */
#define TRPI_CHKRNG     0x000004 /* checkrange() completed */
#define TRPI_GETATM     0x000008 /* getatm() completed */
#define TRPI_MAKERAD    0x000010 /* makeradsample() completed */
#define TRPI_MAKEWAV    0x000020 /* makewavsample() completed */
#define TRPI_MAKEWN     0x000040 /* makewnsample() completed */


/* flags for transiterror */
#define TERR_MESSAGE    0x000000
#define TERR_CRITICAL   0x000001
#define TERR_SERIOUS    0x000002
#define TERR_WARNING    0x000003

#define TERR_NOFLAGBITS 0x00000f

#define TERR_ALLOWCONT  0x000010
#define TERR_NOPREAMBLE 0x000020
#define TERR_ALLOC      0x000040


/*****   Types     *****/
#define PREC_NSAMP int		/* Type for radius and wavelength
				   indices */
#define PREC_NREC long		/* Type for record indices */
#define PREC_ZREC double	/* Type for the partition info  */
#define PREC_LNDATA double	/* Type for the line data output */
#define PREC_RES double     	/* Type for every partial result */
#define PREC_ATM double		/* Type for atmospheric data */
#define PREC_CS  double		/* Type for cross-section */


#ifdef NODEBUG_TRANSIT
#define transitDEBUG(...) ((void)0)
#define transitASSERT(...) ((void)0)
#else
#define transitASSERT(a,...) if(a) transiterror(TERR_CRITICAL,__VA_ARGS__)
#define transitDEBUG(...) transitprint(__VA_ARGS__)
#endif


/***** Constants  *****/
/* units in cgs */
#if defined(AMU) || defined(EC) || defined(LS) || defined(ME) || \
defined(KB) || defined(H) || defined(PI) || defined (SIGCTE) ||  \
defined(EXPCTE) || defined(WNU_O_WLU)
#error Some of the preprocessor constants I wanted to use were already defined elsewhere!
#endif

#define AMU (1.6605402e-24)                 //Atomic Mass unit
#define EC (4.8032068e-10)                  //Electronic charge
#define LS (2.99792458e10)                  //Light Speed
#define ME (9.1093897e-28)                  //Electron mass
#define KB (1.380658e-16)                   //Boltzmann constant
#define H (6.6260755e-27)                   //Planck's constant
#define PI (3.141592653589793)              //PI!
#define RWATER (3.2e-8/2.0)                 //water molecule radius
#define HC (H*LS)                           //for lower energy conversion
#define SIGWATER (PI*RWATER*RWATER)         //water cross section
#define SIGCTE (PI*EC*EC/LS/LS/ME/AMU)      //Cross-sec constant
#define EXPCTE (H*LS/KB)                    //Exponent constant

#define WNU_O_WLU (1e7)              //Waven(cm) over wavel(nm) (units)
#define ONEOSQRT2PI (0.3989422804)   // 1.0/sqrt(2pi)
#define SQRTLN2  0.83255461115769775635 //sqrt(ln(2))

extern const int maxeisoname;
extern int transit_nowarn;
extern int verblevel;
extern int maxline;

/*****   Macros   *****/
#define stateeqnford(q,m,p,t) (AMU*(q)*(m)*(p)/(KB*(t)))

#define transitassert(a,...) if(a) transiterror(TERR_CRITICAL,__VA_ARGS__)
#define transitacceptflag(transit,hint,flag) do{                            \
        transit|=hint&flag;hint&=~(flag);}while(0)
#define transitaccepthint(transit,hint,flags,flagvalue) do{                 \
        transit=hint;                    }while(0)
#define transitprint(thislevel, verblevel, ...) if(thislevel <= verblevel)  \
        fprintf(stderr,__VA_ARGS__)
#define transitallocerror(nmb)                                              \
        transiterror(TERR_CRITICAL,                                         \
	             "transit:: %s: Allocation failed for %i allocation\n"  \
	             "units in line %i. Impossible to continue.\n"          \
	             ,__FILE__,nmb,__LINE__)

enum isodo {unclear=0,atmfile,ignore,fixed};


/***** Structures *****/

typedef struct {          	/* One item per sampling element */
  PREC_NREC n;			/* number of elements */
  PREC_RES d;			/* Spacing */
  PREC_RES i;			/* initial value */
  PREC_RES f;			/* final value */
  int o;			/* oversampling */
  PREC_RES *v;			/* values of the sampling */
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
  PREC_LNDATA wl;		//Wavelength in nm.
  PREC_LNDATA elow;		//Lower energy level in cm-1
  PREC_LNDATA gf;		//gf value
  short isoid;			//Isotope ID (Assumed to be in range)
};


struct atm_data{		/* Keeps parameters in readatminfo() */
  prop_samp rads;		/* radius sampling */
  prop_isov *isov;		/* variable isotope info
				   [isoext] */
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
  int n_aiso;			/* Number of isotopes in the atmosphere
				   file */
};


struct extinction{
  PREC_RES ***k;		/* Extinciton value [rad][iso:][wav]*/
  PREC_VOIGTP **al;		/* Lorentz width [rad][iso] */
  PREC_VOIGTP **ad;		/* Doppler width/central wavenumber
				   [rad][iso] */
  int *lw;			/* number of bins calculated per
				   wavelength [wav] */
  int **recalc;			/* Indicates in how many wavenumbers it
				   have to change profile again, if 0
				   means that it kept profile from
				   previous wavenumber */
  float maxratio;		/* Maximum Doppler width ratio between
				   current and last calculated profile.
				   If the value is greater than this,
				   then recalculate */
  int vf;			/* Fine binning of Voigt function */
  float ta;			/* number of alphas that have to be
				   contained in the profile */
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


struct transithint {		/* Structure with user hinted data that
				   should go to the 'struct transit'
				   upon approval */
  char *f_atm,*f_line,*f_out;	/* Filenames */
  prop_samp rads,wavs,wns;	/* Sampling properties of
				   radius(planetary radius),
				   wavelength(nanometers) and
				   wavenumber(cm-1) */
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
};

struct transit {		/* Main data structure */
  char *f_atm,*f_line,*f_out;	/* Filenames */
  FILE *fp_atm,*fp_out,*fp_line;/* Filepointers */
  float allowrq;		/* How much less than one is accepted,
				   and no warning is issued if
				   abundances don't ad up to that */
  PREC_RES telres;		/* Telescope resolution */
  PREC_RES m;			/* Amount of nanometers not trusted at
				   the boundaries, also how much out of
				   requested range it have to look for
				   transitions */
  PREC_RES wnm;			/* Amount of cm-1 not trusted */
  prop_samp rads,wavs,wns;	/* Sampling properties of radius,
				   wavelength and wavenumber */
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

  long fl;			/* flags */
  long pi;			/* progress indicator */

  struct line_transition *lt;	/* line transition */


  struct {			/* data structures pointers, this is
				   data that is not required for the
				   final computation */
    struct iso_noext *in;
    struct atm_data *at;
    struct transithint *th;	/*  */
    struct lineinfo *li;	/*  */
    struct extinction *ex;
  }ds;
};


/***** Prototypes *****/
#include <transit_proto.h>
#include <readlineinfo_proto.h>
#include <transitstd_proto.h>
#include <readatminfo_proto.h>
#include <makesample_proto.h>
#include <extinction_proto.h>


#endif
