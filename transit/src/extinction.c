/*
 * extinction.c   - Computes extinction coefficient. Component of the
 *                  Transit program.
 *
 * Copyright (C) 2004-2006 Patricio Rojo (pato@astro.cornell.edu)
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

#include <transit.h>

/* If you plan to run extwn() for different atmospheric files without
   rerunning 'transit', then you will need to define EXTINCTION_Q as
   blank, keep it to 'const' otherwise */
#define EXTINCTION_Q const

static _Bool extinctperiso, extwncalledonce=0, gominelow;
static PREC_VOIGT ***profile;
static PREC_RES ***kalt, /* Shape: [nrad, niso, nwn] */
                minelow; 
/* Line info content: */
static EXTINCTION_Q PREC_LNDATA *ltwl;
static EXTINCTION_Q PREC_LNDATA *ltgf, *ltelow;
static EXTINCTION_Q short *ltisoid;
static double efct, wfct;

/* Various radius independent variable: */
static int niso, nmol; //neiso
static PREC_RES *wn, iniwn, wni, wnf, dwn, wavfct;
static PREC_NREC nwn, nlines, nrad;
static struct isotopes *iso;
struct molecules *mol;

/* Following to store isotope data, these are auxiliary variables: */
static PREC_NREC *wa, *wrc, *nwnh;
static PREC_VOIGTP *alphal, *alphad;
static PREC_ATM *densiso;
static PREC_ZREC *ziso, *mass;


/* \fcnfh
   Wrapper to calculate a Voigt profile 'pr'

   Return: 1/2 of the number of points in the non-oversampled profile */
inline int
newprofile(PREC_VOIGT **pr,  /* Output 2D profile               */
           int vf,           /* Number of fine resolution bins  */
           PREC_RES dwn,     /* wavenumber spacing              */
           PREC_VOIGT dop,   /* Doppler width                   */
           PREC_VOIGT lor,   /* Lorentz width                   */
           float ta){        /* times of alpha                  */

  PREC_VOIGTP bigalpha; /* Largest width (Doppler or Lorentz) */
  PREC_VOIGTP wvgt;     /* Calculated half-width of profile   */
  int nvgt,             /* Number of points in profile        */
      j;                /* Auxiliary index                    */

  /* Get the largest width (alpha Doppler or Lorentz): */
  bigalpha = dop;
  if(bigalpha<lor)
    bigalpha = lor;

  /* Half width from line center where the profile will be computed: */
  wvgt = bigalpha*ta;
  /* Number of points in profile: */
  nvgt = 2*(long)(wvgt/dwn+0.5) + 1;

  /* Basic check that 'lor' or 'dop' are within sense: */
  if(nvgt<0)
    transiterror(TERR_CRITICAL,
                 "Number of voigt bins are not positive!.  Doppler width: %g, "
                 "Lorentz width: %g.\n", dop, lor);

  /* Allocate profile array: */
  *pr = (PREC_VOIGT *)calloc(nvgt*vf, sizeof(PREC_VOIGT));
  for(j=1; j<vf; j++)
    pr[j] = pr[0] + j*nvgt;

  /* Calculate voigt using a width that gives an integer number of 'dwn'
     spaced bins: */
  if((j=voigtn(vf, nvgt, dwn*(long)(nvgt/2), lor, dop, pr, -1, 
               nvgt>_voigt_maxelements?VOIGT_QUICK:0)) != 1)
    transiterror(TERR_CRITICAL, "voigtn() returned error code %i.\n", j);

  return nvgt/2;
}


/* \fcnfh
   Compute radius reordering extinction array to be processed by
   extradius

   Return: 0 on success, else extradius()                     */
int
computeextradius(PREC_NREC r,            /* Radius index      */
                 PREC_ATM temp,          /* Temperature       */
                 struct extinction *ex){ /* Extinction struct */
  int res;
  if((res=extradius(r, kalt[r], temp, ex->vf, ex->ta, ex->maxratio))==0)
    ex->computed[r] = 1;

  return res;
}


/*\fcnfh
  Compute extinction coefficient at one specific radius.

  Return: 0 on success, else
         -2 if no radii has been specified.
         -3 if only 2 wavelengths points
         -5 if no isotopes selected!
         -6 if fine binning less than 1
         -7 if timesalpha less than 1
         -8 if maxratio less than 0                                  */
int
extradius(PREC_NREC r,      /* Radius index                          */
          PREC_RES **kiso,  /* Extinction [iso][wn]                  */
          PREC_ATM temp,    /* Temperature                           */
          int fbinvoigt,    /* Fine-bins for voigt profile           */
          float timesalpha, /* Number of Voigt half-widths computed  */
          float maxratio){  /* Maximum allowed Doppler radius before
                                   recalculating profile             */

  PREC_NREC ln;
  int i, imol;
  long j, maxj, minj;
  PREC_VOIGT *profwn;
  PREC_NREC w, subw;
  PREC_RES *k, wavn;
  double csdiameter; /* Collision diameter */
  double propto_k;
  double maxk = 0;

  if(!extwncalledonce)
    transiterror(TERR_CRITICAL, "Trying to call extradius before extwn.\n");

  /* Constant proportionality factor of the Doppler width (in
     sqrt(AMU) units):
      \alpha_D = \frac{\nu}{\sqrt{m}}
         \underbrace{\frac{\sqrt{2k_B T \ln{2}}}{c}}_{\rm{propto\_adop}}
      \label{dopbr}                                               */
  double propto_adop = sqrt(2*KB*temp/AMU) * SQRTLN2/LS;

  /* Constant proportionality factor of the Lorenz width (FINDME: units?):
     \alpha_L = \underbrace{\Gamma_{nat}}_{\mathrm{ignored}} +
                \underbrace{\sqrt{\frac{2k_B T}{\pi c^2}}}_{\rm{propto\_alor}}
                \sum_{i}n_i\sigma_{ir}^2
                      \sqrt{\left(\frac{1}{m_r} + \frac{1}{m_i}\right)}
     \label{lorwidth}                                                 */
  double propto_alor = sqrt(2*KB*temp/PI/AMU)/(AMU*LS);

  /* Initialize a Voigt profile for every isotope as well for
     the mass, ziso, and densiso arrays:                      */
  for (i=0; i<nmol; i++)
    densiso[i] = mol->molec[i].d[r];
  for(i=0; i<niso; i++){    /* Initialize arrays:             */
    mass[i]    = iso->isof[i].m;
    ziso[i]    = iso->isov[i].z[r];

    imol = iso->imol[i];  /* Isotope's molecule index */
    /* Calculate Lorentz line width: */
    alphal[i] = 0.0;
    for(j=0; j<nmol; j++){
      csdiameter = (mol->radius[j] + mol->radius[imol]);
      transitprint(20, verblevel, "Iso's coll diameter square %.21f.\n"
                                  "Iso's cross section/pi:  %.21f.\n",
                                  csdiameter*csdiameter, iso->isov[i].c[r]/PI);
      alphal[i] += mol->molec[j].d[r]/mol->mass[j] * csdiameter * csdiameter *
                   sqrt(1/mass[i] + 1/mol->mass[j]);
    }
    alphal[i] *= propto_alor;

    /* Calculate Doppler width divided by central wavenumber: */
    alphad[i] = propto_adop/sqrt(mass[i]);

    transitprint(30, verblevel, "Lorentz: %.9f, Doppler: %.9f broadening.\n",
                               alphal[i], alphad[i]);
    /* Get a new profile: 'profile[i]' has dimensions ex->vf x (2*nwnh[i]+1).
       The latter is calculated by newprofile, though: */
    if((nwnh[i] = newprofile(profile[i], fbinvoigt, dwn,
                             wn[0]*alphad[i], alphal[i], timesalpha)) < 1)
      transiterror(TERR_CRITICAL, "newprofile() returned error code %i on its "
                                  "first try for isotope %i.\n", nwnh[i], i);

    w = nwn-1; /* Index of last wavenumber */
    /* How many other wavenumbers bins to recalculate the Voigt profile: */
    j = (int)(maxratio*wn[w]/dwn + 0.5);
    if(!j) j = 1;  /* Has to be at least 1 */
    /* Wavenumber index where the profile was last calculated:  */
    wa[i] = w;
    /* Wavenumber index where the profile will be recalculated: */
    wrc[i] = w-j;
  }

  /* Compute the spectra, proceed for every line: */
  for(ln=0; ln<nlines; ln++){
    /* Skip lines which strength is lower than minelow:          */
    if(gominelow && ltelow[ln]<minelow)
      continue;  /* FINDME: Defining gominelow seems unnecessary */

    /* Wavenumber of line transition: */
    wavn = 1.0/(ltwl[ln]*wfct);

    /* If it is beyond lower limit, skip to next line transition: */
    if(wavn < iniwn)
      continue;
    else
      /* Closest wavenumber index to transition, but no larger than: */
      w = (wavn-iniwn)/dwn;
    transitDEBUG(25, verblevel, "wavn: %g,  lgf: %g.\n", wavn, ltgf[ln]);

    /* If it is beyond the last index, skip to next line transition: */
    if(w >= nwn)
      continue;

    /* Distance from center of line to sampled wavenumber by 'w', in number
       of fine bins:  */
    subw = fbinvoigt*(wavn - w*dwn - iniwn)/dwn;

    i = ltisoid[ln]; /* Isotope ID of line     */
    k = kiso[i];     /* Extinction coefficient array[wn] for given isotope */

    /* FINDME: Comment me: */
    transitASSERT(wa[i] != -1 && wa[i]<w,
                  "Database is not ordered!, previous wavenumber was at index "
                  "%i, new one at %i (it should have been smaller).\n",
                  wa[i], w);

    /* If w <= wrc, recalculate Voigt profile: */
    if(w<=wrc[i]){
      /* Calculate next index where to recalculate the profile: */
      j = (int)(maxratio*(wn[w])/dwn+0.5);
      if(!j) j = 1;
      wrc[i] = w-j;
      transitDEBUG(22, verblevel,
                   "Recalculating Voigt for isotope %i ... current wavenumber "
                   "%li, next wavenumber %li/%li.\n", i, w, wrc[i], nwn);

      /* Free old profile and recalculate Voigt: */
      free(profile[i][0]);
      if((nwnh[i] = newprofile(profile[i], fbinvoigt, dwn,
                               wn[w]*alphad[i], alphal[i], timesalpha)) < 1)
        transiterror(TERR_CRITICAL, "newprofile() returned error code %i for "
                                    "isotope %i.\n", nwnh[i], i);
    }

    /* TD: SAHA equation is missing in level population!! Code does not
       consider ionizable levels of species. Error in thesis' equation 3.35! */

    /* Calculate opacity coefficient except the broadening factor: */
    propto_k = densiso[iso->imol[i]] * iso->isoratio[i] * /* Density     */
               SIGCTE     * ltgf[ln]             * /* Constant * gf      */
               exp(-EXPCTE*efct*ltelow[ln]/temp) * /* Level population   */
               (1-exp(-EXPCTE*wavn/temp))        / /* Induced emission   */
               mass[i]                           / /* Mass               */
               ziso[i];                            /* Partition function */

    /* Maximum line opacity coefficient at central wavenumber: */
    if (propto_k > maxk)
      maxk = propto_k;

    transitDEBUG(24, verblevel,
                 "i=%i   temp=%g   Elow=%g\n"
                 "aD=%.7g   aL=%.7g\n"
                 "wl=%.10g  wn=%.10g\n"
                 "k =%12.5g   // densiso[imol]\n"
                 "  *%12.5g   // isoratio\n"
                 "  *%12.5g   // SIGCTE\n"
                 "  *%12.5g   // ltgf[ln]\n"
                 "  *%12.5g   // exp(-EXPCTE*ltelow[ln]/temp)\n"
                 "  *%12.5g   // (1-exp(-EXPCTE*wavn/temp))\n"
                 "  /%12.5g   // mass[i]\n"
                 "  /%12.5g   // ziso[i]\n"
                 " = %12.5g   // extinction\n\n",
                 i, temp, ltelow[ln],
                 alphad[i]*wavn, alphal[i],
                 ltwl[ln], 1.0/(wfct*ltwl[ln]*wavfct),
                 densiso[iso->imol[i]],
                 iso->isoratio[i],
                 SIGCTE,
                 ltgf[ln],
                 exp(-EXPCTE*ltelow[ln]/temp),
                 (1-exp(-EXPCTE*wavn/temp)),
                 mass[i],
                 ziso[i],
                 propto_k);

    /* nwnh is (number of points in the profile)/2.                */
    /* nwnh-w-1 is the starting point of the wavenumber array with
       respect to the starting index of the profile.               */
    /* Set profwn such that the index mimic wavenumber's array:    */
    profwn = profile[i][subw] + nwnh[i]-w-1;

    /* Set the lower and upper indices of the profile to be used:  */
    minj = w-nwnh[i]+1;
    if(minj<0)
      minj = 0;

    maxj = w+nwnh[i]+1;
    if(maxj>nwn)
      maxj = nwn;

    /* Distribute the oscillator strength according to the Voigt profile: */
    for(j=minj; j<maxj; j++)
      k[j] += propto_k*profwn[j];

    if (r==43 && verblevel==21)
      printf("%-9li%-20.9g%-20.9g%-20.9g\n", ln, wavn, ltgf[ln], k[5763]);

#if 0
    if(ltwl[ln]>1696.8267){
      k=5;
    }
    else if(ln==2594618){
      k=5;
    }
    else if(ln==2595091){
      k=5;
    }
#endif

    wa[i] = w;  /* Last wavelength per isotope */
  }

  /* Free the profiles of every non-ignored isotopes: */
  for(i=0; i<niso; i++)
      free(profile[i][0]);
  transitprint(2, verblevel, "Done.\n");
  return 0;
}


/* \fcnfh
 Saving extinction for a possible next run
*/
void
savefile_extinct(char *filename,
                 PREC_RES **e,
                 _Bool *c,
                 long nrad,
                 long nwav){

  FILE *fp;

  if((fp=fopen(filename, "w")) == NULL){
    transiterror(TERR_WARNING,
                 "Extinction savefile '%s' cannot be opened for writing.\n"
                 " Continuing without saving\n"
                 ,filename);
    return;
  }

  transitprint(2, verblevel, "Saving extinction file '%s'", filename);

  const char mn[] = "@E@S@";
  fwrite(mn, sizeof(char), 5, fp);
  fwrite(e[0], sizeof(PREC_RES), nrad*nwav, fp);
  fwrite(c, sizeof(_Bool), nrad, fp);

  fclose(fp);

  int i;
  for (i=0 ; i<nrad ; i++)
    if (c[i]) break;

  transitprint(2, verblevel, " done (%li/%li radii computed)\n", nrad-i, nrad);
}


/* \fcnfh
   Restoring extinction for a possible next run
*/
void
restfile_extinct(char *filename,
                 PREC_RES **e,
                 _Bool *c,
                 long nrad,
                 long nwav){

  FILE *fp;

  if((fp=fopen(filename, "r")) == NULL){
    transiterror(TERR_WARNING,
                 "Extinction savefile '%s' cannot be opened for reading.\n"
                 "Continuing without restoring. You can safely ignore "
                 "this warning if this the first time you run for this "
                 "extinction savefile.\n"
                 ,filename);
    return;
  }

  char mn[5];
  if(fread(mn, sizeof(char), 5, fp)!=5 || strncmp(mn,"@E@S@", 5)!=0){
     transiterror(TERR_WARNING,
                  "Given filename for extinction savefile '%s' exists\n"
                  "and is not a valid extinction file. Remove it\n"
                  "before trying to use extinction savefile\n"
                  ,filename);
     fclose(fp);
     return;
  }

  transitprint(2, verblevel, "Restoring extinction file '%s'", filename);

  fread(e[0], sizeof(PREC_RES), nrad*nwav, fp);
  fread(c, sizeof(_Bool), nrad, fp);

  long i;
  for (i=0 ; i<nrad ; i++)
    if (c[i]) break;

  transitprint(2, verblevel,
               " done (From the %lith radii)\n"
               , i);

  fclose(fp);

}



/* \fcnfh
   Output info to file regarding extinction computation. Designed to be
   called from the debugger only

   Example: calling from computeextradius():
    outputinfo("out", w, 5, ln, 500, kiso, timesalpha, fbinvoigt, temp, rad); */
void
outputinfo(char *outfile,
           long w,
           long dw,
           long ln,
           long dln,
           double **kiso,
           double timesalpha,
           int fbinvoigt,
           double temp,
           double rad){
  long i, maxnwn=0;
  int j;

  FILE *out = fopen(outfile, "w");
  if(!out){
    fprintf(stderr, "Cannot write to file '%s'.\n", outfile);
    return;
  }

  dw += w;
  fprintf(out, "Debuging output:\n"
               "Radius: %g.\nTemperature: %g.\n"
               "Number_of_extinction_points: %li.\n"
               "Number_of_line_info: %li.\n"
               "timesalpha: %.9g.\nNumber of finebins: %i.\n",
               rad, temp, dw-w, dln, timesalpha, fbinvoigt);
  fprintf(out, "--------------------------------------------------\n"
               "Procesed info from index %li to %li, for each isotope",
               w, dw-1);

  for(i=w; i<dw; i++){
    fprintf(out, "\n%-15.9g", wn[i]);
    for(j=0; j<niso; j++)
      fprintf(out, "%-15.9g", kiso[j][i]);
  }

  fprintf(out, "\n--------------------------------------------------\n"
                 "Info from Doppler, next recalculation will occur  \n"
                 "at the following wavelengths for the %i different \n"
                 "isotopes", niso);
  for(i=0; i<niso; i++)
    fprintf(out, "%15.9g(%li)  ", wn[wrc[i]], wrc[i]);
  fprintf(out, "\nApprox_Doppler    Lorentz     #elem_width\n");

  for(i=0; i<niso; i++){
    fprintf(out, " %-18.9g%-15.9g%-15li\n",
            wn[w]*alphad[i], alphal[i], nwnh[i]*2+1);
    if(maxnwn<nwnh[i])
      maxnwn = nwnh[i];
  }
  fprintf(out, "Doppler profile (shown by finebinning and then by isotope):\n");
  int k;
  for(  i=0; i<=maxnwn; i++){
    for(j=0; j< niso;   j++){
      if(i<=nwnh[j])
        for(k=0; k<fbinvoigt; k++)
          fprintf(out, "%-15.9g", profile[j][k][i]);
      fprintf(out, " | ");
    }
    fprintf(out, "\n");
  }
  fprintf(out,
          "---------------------------------------------------\n"
          "Line information, showing %li lines:\n"
          "index       Wavenumber-cm    Wavelength-nm    GF             "
          "Elow         Iso\n", dln);
  for(i=ln; i<ln+dln; i++)
    fprintf(out, "%-11li%-15.9g%-15.9g%-15.9g%-15.9g%2i\n",
            i, 1e7/ltwl[i], ltwl[i], ltgf[i], ltelow[i], ltisoid[i]);

  fclose(out);
}



/*\fcnfh
  Fill up the extinction information in tr->ds.ex  
  TD: Scattering parameters should be added at some point here.

  Return: 0 on success, else
          computeextradius()  */
int
extwn(struct transit *tr){
  static struct extinction st_ex;
  tr->ds.ex = &st_ex;
  struct extinction *ex = &st_ex;
  int i, j;

  /* Check these routines have been called: */
  transitcheckcalled(tr->pi, "extwn", 4,
                             "readinfo_tli",  TRPI_READINFO,
                             "readdatarng",   TRPI_READDATA,
                             "makewnsample",  TRPI_MAKEWN,
                             "makeradsample", TRPI_MAKERAD);
  transitacceptflag(tr->fl, tr->ds.th->fl, TRU_EXTBITS);

  /* Check that voigtfine (oversampling factor of the profile binning)
     is > 1, and set it's value in transit: */
  if(tr->ds.th->voigtfine < 1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Fine binning of Voigt function has to be positive: %i.\n",
                 tr->ds.th->voigtfine);
    return -6;
  }
  ex->vf = tr->ds.th->voigtfine;

  /* Check that timesalpha (profile half width in units of Doppler/Lorentz
     broadening half widths) is > 1, and set it's value in transit: */
  if(tr->ds.th->timesalpha < 1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Times of maximum width has to be greater than one: %i\n",
                 tr->ds.th->voigtfine);
    return -7;
  }
  ex->ta = tr->ds.th->timesalpha;

  /* Copy minelow into ex if value is positive: */
  gominelow = tr->ds.th->minelow>0;
  if(gominelow)
    minelow = ex->minelow = tr->ds.th->minelow;

  /* maxratio is the maximum allowed ratio change before recalculating
     profile array: */
  if(tr->ds.th->maxratio_doppler<0){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Maximum allowed Doppler width ratio change. Has to "
                 "be 0 or positive (%g).\n", tr->ds.th->maxratio_doppler);
    return -8;
  }
  ex->maxratio = tr->ds.th->maxratio_doppler;

  /* Get line-transition info: */
  struct line_transition *lt = &(tr->ds.li->lt);
  ltgf    = lt->gf;
  ltelow  = lt->elow;
  ltwl    = lt->wl;
  ltisoid = lt->isoid;
  efct    = lt->efct;
  wfct    = lt->wfct;
  nlines  = tr->ds.li->n_l;
  nrad    = tr->rads.n;

  /* Radius-independent variables: */
  iso    = tr->ds.iso;
  mol    = tr->ds.mol;
  nmol   = mol->nmol;
  niso   = iso->n_i;
  wavfct = tr->wavs.fct;  /* FINDME: typo, should be wavs.fct? */
  wn     = tr->wns.v;
  nwn    = tr->wns.n;
  iniwn  = tr->wns.i;
  wni    = wn[    0]-tr->wnmi;
  wnf    = wn[nwn-1]+tr->wnmf;
  dwn    = tr->wns.d/tr->wns.o;
  extwncalledonce = 1;

  /* Allocate line-profile info: */
  alphal  = (PREC_VOIGTP *)calloc(niso*2, sizeof(PREC_VOIGTP));
  ziso    = (PREC_ZREC   *)calloc(niso*2, sizeof(PREC_ZREC)  );
  densiso = (PREC_ATM    *)calloc(nmol  , sizeof(PREC_ATM)   );
  wa      = (PREC_NREC   *)calloc(niso*3, sizeof(PREC_NREC)  );
  wrc    = wa     + niso;
  nwnh   = wa     + niso*2;
  alphad = alphal + niso;
  mass   = ziso   + niso;

  /* TD: enable nisoalloc, currently alloc is used  */
  /* Do not allocate memory in multidimensional arrays if we are
     ignoring that particular isotope. We are not using this in
     unidimensional arrays because of the extra hassle:           */
  //int nisoalloc = niso; /* FINDME: Number of isotopes to alloc? for what?*/
  //for(i=0; i<niso; i++)
  //  if(iso->isodo[i]==ignore)
  //    nisoalloc--;
  /* FINDME: nisoalloc unused */

  /* Check there is at least one atmospheric layer:        */
  if(nrad<1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "There are no atmospheric parameters specified. I need at "
                 "least one atmospheric point to calculate a spectra.\n");
    return -2;
  }
  /* Check there are at least two wavenumber sample points: */
  if(nwn<2){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "I need at least 2 wavenumber points to compute "
                 "anything; I need resolution.\n");
    return -3;
  }
  /* Check there is at least one isotope linelist          */
  /* FINDME: This should not be a condition, we may want
     to calculate an atmosphere with only CIA for example. */
  if(niso < 1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "You are requiring a spectra of zero isotopes!.\n");
    return -5;
  }

  /* Allocate array for the Voigt profile, also auxiliary: */
  profile  = (PREC_VOIGT ***)calloc(niso,        sizeof(PREC_VOIGT **));
  *profile = (PREC_VOIGT  **)calloc(niso*ex->vf, sizeof(PREC_VOIGT *));
  for(i=1; i<niso; i++)
    profile[i] = *profile+i*ex->vf;

  /* Set extinction-per-isotope boolean: */
  ex->periso = extinctperiso = ((tr->fl&TRU_EXTINPERISO)==TRU_EXTINPERISO);

  /* Arrange the extinctions so that the order is [iso][rad][wn] */
  int nni = extinctperiso?niso:1;
  int nnr = extinctperiso?nrad:0;

  /* FINDME: Should this be nni instead of niso? */
  ex->e           = (PREC_RES ***)calloc(niso,         sizeof(PREC_RES **));
  ex->e[0]        = (PREC_RES  **)calloc(nni*nrad,     sizeof(PREC_RES *));
  if((ex->e[0][0] = (PREC_RES   *)calloc(nni*nrad*nwn, sizeof(PREC_RES)))==NULL)
    transiterror(TERR_CRITICAL|TERR_ALLOC,
                 "Unable to allocate %li = %li*%li*%li to calculate "
                 "extinction for every radii, %stry to shorten the  "
                 "wavenumber range.\n", nrad*i*nwn, nrad, i, nwn,
                 extinctperiso?"try disabling exctinction per isotope "
                               "(option --no-per-iso), or ":"");
  for(i=0; i<niso; i++){
    ex->e[i] = ex->e[0]+i*nnr;
    if(!i || extinctperiso)
      for(j=0; j<nrad; j++)
        ex->e[i][j] = ex->e[0][0] + nwn*(j + nnr*i);
  }

  /* Have a local array which will have extinction ready to be
     used by extradius:                                         */
  kalt    = (PREC_RES ***)calloc(nrad,      sizeof(PREC_RES **));
  kalt[0] = (PREC_RES  **)calloc(nrad*niso, sizeof(PREC_RES  *));
  for(i=0; i<nrad; i++){
    kalt[i] = kalt[0]+i*niso;
    for(j=0; j<niso; j++)
      kalt[i][j] = ex->e[j][i]; /* Note, swapped indices */
  }
  
  /* Has the extinction been computed at given radius boolean: */
  ex->computed = (_Bool *)calloc(nrad, sizeof(_Bool));

  /* For each radius (index 'r'): */
  transitprint(1, verblevel, "\nThere are %li radii samples.\n", nrad);

  /* Save current status if requested: */
  //  savefile_extwn(tr); 

  /* Set progress indicator, and print and output extinction if one P,T
     was desired, otherwise return success: */
  tr->pi |= TRPI_EXTWN;
  if(tr->rads.n == 1)
    printone(tr);
  return 0;
}


/* \fcnfh
   Printout for one P,T conditions */
void
printone(struct transit *tr){
  int rn;
  FILE *out=stdout;

  /* Open file: */
  if(tr->f_out && tr->f_out[0] != '-')
    out = fopen(tr->f_out, "w");

  transitprint(1, verblevel, "\nPrinting extinction for one radius (at %gcm) "
                             "in '%s'\n", tr->rads.v[0],
                             tr->f_out?tr->f_out:"standard output");

  /* Print: */
  fprintf(out, "#wavenumber[cm-1]   wavelength[nm]   extinction[cm-1]   "
               "cross-section[cm2]\n");
  for(rn=0; rn < tr->wns.n; rn++)
    fprintf(out, "%12.6f%14.6f%17.7g%17.7g\n", tr->wns.fct*tr->wns.v[rn],
            1/(tr->wavs.fct * tr->wns.v[rn] * tr->wns.fct),
            tr->ds.ex->e[0][0][rn],
            AMU*tr->ds.ex->e[0][0][rn] * iso->isof[0].m/mol->molec[iso->imol[0]].d[0]);

  exit(EXIT_SUCCESS);
}


/* \fcnfh
   Free extinction coefficient structure arrays

   Return 0 on success */
int
freemem_extinction(struct extinction *ex, /* Extinciton struct       */
                   long *pi){             /* progress indicator flag */
  /* Free arrays: */
  free(ex->e[0][0]);
  free(ex->e[0]);
  free(ex->e);
  free(ex->computed);

  /* Update indicator and return: */
  *pi &= ~(TRPI_EXTWN);
  return 0;
}


/* \fcnfh
   Restore hints structure, the structure needs to have been allocated
   before

   @returns 0 on success
            -1 if not all the expected information is read
            -2 if info read is wrong
            -3 if cannot allocate memory
            1 if information read was suspicious
*/
int
restextinct(FILE *in,
            PREC_NREC nrad,
            short niso,
            PREC_NREC nwn,
            struct extinction *ex){
  PREC_NREC nr,nw;
  short i;
  size_t rn;

  rn=fread(&nr,sizeof(PREC_NREC),1,in);
  if(rn!=1) return -1;
  rn=fread(&i,sizeof(short),1,in);
  if(rn!=1) return -1;
  rn=fread(&nw,sizeof(PREC_NREC),1,in);
  if(rn!=1) return -1;

  //no more than 10 000 isotopes, or 10 000 000 of radii or wavelenth
  if(nr!=nrad||nw!=nwn||i!=niso||
     niso>10000||nrad>10000000||nwn>10000000)
    return -2;

  int nni=ex->periso?niso:1;

  if((ex->e      =(PREC_RES ***)calloc(niso        ,sizeof(PREC_RES **)))
     ==NULL)
    return -3;
  if((ex->e[0]   =(PREC_RES **) calloc(nni*nrad    ,sizeof(PREC_RES * )))
     ==NULL)
    return -3;
  if((ex->e[0][0]=(PREC_RES *)  calloc(nrad*nni*nwn,sizeof(PREC_RES   )))
     ==NULL)
    return -3;

  rn=fread(ex->e[0][0],sizeof(PREC_RES),nwn*nni*nrad,in);
  if(rn!=nwn*nni*nrad) return -1;
  rn=fread(ex->computed,sizeof(PREC_RES),nrad,in);
  if(rn!=nrad) return -1;

  int nnr = ex->periso?nrad:0;
  for(i=0;i<niso;i++){
    ex->e[i]=ex->e[0]+i*nnr;
    if(!i||ex->periso)
      for(nr=0;nr<nrad;nr++)
        ex->e[i][nr]=ex->e[0][0] + nwn*( nr + nni*i );
  }

  return 0;
}


#if 0
/* \fcnfh
   Save program memory after being through tau computation, whether
   ->save.tau is defined have to be checked before.
*/
void
savefile_exsofar(struct transit *tr)        /* Main structure */
{
  FILE *fp;

  if((fp=fopen(tr->save.tau,"w"))==NULL){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "\nCannot open file '%s' for saving extinction data.\n"
                 " Continuing without saving\n\n"
                 ,tr->save.tau);
    return;
  }

  const char *idchar=__TR_SAVEFILE_MN__"extinctionsofar";
  fwrite(idchar,1,strlen(idchar),fp);

  /* TD: following function is not written! */
  saveline(fp,tr->ds.li);
  savesample(fp,&tr->rads);
  savesample(fp,&tr->wns);
  saveextinct(fp,tr->rads.n,tr->ds.iso->n_i,tr->wns.n,tr->ds.ex);

  /*                     "readinfo_tli",TRPI_READINFO,
                     "readdatarng",TRPI_READDATA,*/
                       

  fclose(fp);
}


/* \fcnfh
   Restore program memory of calculated extinction, whether
   ->save.tau is defined have to be checked before.
*/
void
restorefile_exsofar(struct transit *tr)        /* Main structure */
{
  FILE *fp;

  if((fp=fopen(tr->save.tau,"r"))==NULL){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "\nCannot open file '%s' for restoring extinction data.\n"
                 " Continuing anyways\n\n"
                 ,tr->save.tau);
    return;
  }

  const char idchar[]=__TR_SAVEFILE_MN__"extinctionsofar";
  int len=strlen(idchar);
  char nidchar[len];
  fread(nidchar,1,len,fp);

  if(strncmp(nidchar,idchar,len)!=0){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "\nFile '%s' is invalid for restoring extinction data.\n"
                 " Continuing anyways\n\n"
                 ,tr->save.tau);
    return;
  }

  restextinct(fp,tr->rads.n,tr->ds.iso->n_i,tr->wns.n,tr->ds.ex);

  fclose(fp);
}


/* \fcnfh
   Saves index of refraction structure
*/
void
saveextinct(FILE *out,
            PREC_NREC nrad,
            short niso,
            PREC_NREC nwn,
            struct extinction *ex)
{
  //save structure
  fwrite(ex,sizeof(struct extinction),1,out);

  int n=ex->periso?niso:1;

  fwrite(&nrad,sizeof(PREC_NREC),1,out);
  fwrite(&niso,sizeof(short),1,out);
  fwrite(&nwn,sizeof(PREC_NREC),1,out);

  if(n<=0&&nrad<=0&&nwn<=0)
    transiterror(TERR_SERIOUS,
                 "Quantities of all of radius(%g), isotope(%i) and wavenumber(%i)\n"
                 "have to be bigger than 1\n"
                 ,nrad,niso,nwn);
  fwrite(ex->e[0][0],sizeof(PREC_RES),nrad*n*nwn,out);
  fwrite(ex->computed,sizeof(_Bool),nrad,out);
}



/* \fcnfh 
   Saves all the memory that is going to be used after running extwn()

   Return: 0 on success                                                  */
int
savefile_extwn(struct transit *tr){
  return 0;
}

#endif


/* \fcnfh
   Frees voigt profile pointer arrays. Data array should already be free  */
void
freemem_localextinction(){
  /* profile[0][0] will always be freed after a succesfull completion
     of extradius: */
  free(profile[0]);
  free(profile);

  /* kalt[0][0] will actually be tr.ds->ex->e[0][0], which is not in
     charge of freemem_extinction: */
  free(kalt[0]);
  free(kalt);

  /* Free auxiliar variables that were allocated to number of isotopes: */
  free(wa);
  free(alphal);
  free(ziso);
  free(densiso);
}
