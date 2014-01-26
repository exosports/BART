/*
 * tau.c   - Finds the optical depth. Component of the Transit program.
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

#include <transit.h>
#include <extraext.h>

#define CIA_DOFLOAT  2
#define CIA_RADFIRST 1

/* \fcnfh
   Calculates optical depth as a function of radii for a spherically
   symmetric planet.

   Return: 0 on success   */
int
tau(struct transit *tr){
  /* Different structures */
  static struct optdepth tau;        /* Optical depth struct    */
  tr->ds.tau     = &tau;             /* Optical depth           */
  prop_samp *rad = &tr->rads;        /* Radius sample           */
  prop_samp *wn  = &tr->wns;         /* Wavenumber sample       */
  prop_samp *ip  = &tr->ips;         /* Impact parameter sample */
  struct extinction *ex = tr->ds.ex; /* Extinction struct       */
  PREC_RES **e = ex->e[tr->tauiso];  /* Extinction coefficient  */
  PREC_RES (*fcn)() = tr->sol->tauperb; /* transit ray path?    */

  long wi, ri; /* Indices for wavenumber, and radius */
  int rn;

  /* Number of elements for each array:       */
  long int wnn = wn->n;   /* Wavenumber       */
  long int inn = ip->n;   /* Impact parameter */
  long int rnn = rad->n;  /* Radius           */
  double wfct  = wn->fct; /* Wavenumber units factor */

  PREC_RES *bb = ip->v;          /* Impact parameter array */
  PREC_RES *n  = tr->ds.ir->n;   /* Index of refraction    */
  PREC_RES *r  = rad->v;         /* radius array           */
  PREC_RES *tau_wn;              /* optical depth          */
  PREC_ATM *temp = tr->atm.t,    /* Temperatures           */
           tfct  = tr->atm.tfct; /* Temperature units      */ 

  PREC_RES er[rnn];        /* Array of extinction per radius           */
  int lastr = rnn-1;       /* Radius index of last computed extinction */
  int wnextout = (long)(wnn/10.0); /* (Wavenumber sample size)/10,
                                      used for progress printing       */

  /* Get a copy of the radius units factor, and  the impact
     parameter-to-radius units factor ratio:                     */
  double rfct = rad->fct;
  double riw  = ip->fct/rfct;

  double e_s[rnn],                  /* Extinction from scattering   */
         e_c[rnn];                  /* Extinction from clouds       */
  PREC_CIA **e_cia = tr->ds.cia->e; /* Extinction from CIA          */
  struct extscat *sc = tr->ds.sc;   /* Scattering extinction struct */

  /* Check idxrefrac and extwn have been called:     */
  transitcheckcalled(tr->pi, "tau", 2, "idxrefrac", TRPI_IDXREFRAC,
                                       "extwn",     TRPI_EXTWN);

  /* FINDME: what's going on here? */
  transitacceptflag(tr->fl, tr->ds.th->fl, TRU_TAUBITS);

  /* Get optical depth calculation parameters:                      */
  struct transithint *th = tr->ds.th;
  /* Filename to save/restore the extinction coefficient:           */
  tr->save.ext = th->save.ext;
  /* Line strength enhance factor:                                  */
  const double blowex = tr->blowex   = th->blowex;
  /* Use constant (taulevel=1) or variable (2) index of refraction: */
  const int taulevel  = tr->taulevel = th->taulevel;
  /* Set transit maximum optical depth to calculate:                */
  tau.toomuch = 50;  /* Default value */
  if(tr->ds.th->toomuch > 0)
    tau.toomuch = th->toomuch; /* FINDME: Set default in argum.c */

  /* Radius index where tau reached toomuch: */
  tau.last = (long      *)calloc(wnn,       sizeof(long));
  /* Optical depth per impact parameter:     */
  tau.t    = (PREC_RES **)calloc(wnn,       sizeof(PREC_RES *));
  tau.t[0] = (PREC_RES  *)calloc(wnn*ip->n, sizeof(PREC_RES));
  for(wi=1; wi<wnn; wi++)
    tau.t[wi] = tau.t[0] + wi*ip->n;

  /* Set cloud structure:                           */
  static struct extcloud cl;
  cl.maxe = tr->ds.th->cl.maxe; /* Maximum opacity  */
  cl.rini = tr->ds.th->cl.rini; /* Top layer radius */
  cl.rfin = tr->ds.th->cl.rfin; /* Radius of maxe   */
  if(tr->ds.th->cl.rfct==0)
    cl.rfct = rad->fct;
  else
    cl.rfct = tr->ds.th->cl.rfct;
  tr->ds.cl = &cl;

  /* Has the extinction coefficient been calculated boolean: */
  _Bool *comp = ex->computed;
  /* Restore extinction savefile if exists: */
  if(tr->save.ext)
    restfile_extinct(tr->save.ext, e, comp, rnn, wnn);

  /* Compute extinction at the outermost layer: */
  if(!comp[rnn-1]){
    transitprint(1, verblevel, "Computing extinction in the "
                               "outermost layer.\n");
    if((rn=computeextradius(rnn-1, tr->atm.t[rnn-1]*tr->atm.tfct, ex))!=0)
      transiterror(TERR_CRITICAL,
                   "computeextradius() returned error code %i.\n", rn);
  }

  /* Request at least four impact parameter samples to calculate
     a spline interpolation:                                     */
  if(inn<4)
    transiterror(TERR_SERIOUS,
                 "tau(): At least four impact parameter points are "
                 "required (three for spline and one for the analitical "
                 "part)!");

  transitprint(1, verblevel, "Calculating optical depth at various "
                             "radii ...\n");

  /* FINDME: If periso is True, it will calculate only one istope's
     optical depth? This is different to what was advertised in argum.c ! */
  if(ex->periso)
    transitprint(2, verblevel,
                 "Computing only for isotope '%s', others were ignored.\n",
                 tr->ds.iso->isof[tr->tauiso].n);

  /* For each wavenumber: */
  for(wi=0; wi<wnn; wi++){
    tau_wn = tau.t[wi];

    /* Print output every 10% progress: */
    if(wi>wnextout){
      transitprint(2, verblevel, "%i%%\n", (int)(100*(float)wi/wnn+0.5));
      wnextout += (long)(wnn/10.0);
    }

    /* Calculate extinction from scattering, clouds, and CIA at each level: */
    computeextscat(e_s,  rnn, sc, rad->v, rad->fct, temp, tfct, wn->v[wi]*wfct);
    computeextcloud(e_c, rnn, &cl, rad, temp, tfct, wn->v[wi]*wfct);

    /* Put the extinction values in a new array, the values may be
       temporarily overwritten by (fnc)(), but they should be restored: */
    for(ri=0; ri<rnn; ri++)
      er[ri] = e[ri][wi]*blowex + e_s[ri] + e_c[ri] + e_cia[wi][ri];

    /* FINDME: Explain this line: */
    if(wi==300)
      rn = 56;

    /* For each impact parameter: */
    for(ri=0; ri<inn; ri++){
      /* Compute extinction at new radius if the impact parameter is smaller
         than the radius of last calculated extinction:                    */
      if(bb[ri]*ip->fct < r[lastr]*rfct){
      /* FINDME: What if the ray ends up going through a lower layer because
         of the bending? */
        if(ri)
          transitprint(3, verblevel, "Last Tau (bb=%9.4g, wn=%9.4g): %10.4g.\n",
                                     bb[ri-1], wn->v[wi], tau_wn[ri-1]);
        /* While the extinction at a radius bigger than the impact
           parameter is not computed, go for it: */
        do{
          if(!comp[--lastr]){
            /* Compute a new extinction at given radius printing error if
               something happen: */
            transitprint(2, verblevel, "Radius %i: %.9g cm ... ",
                                       lastr+1, r[lastr]);
            if((rn=computeextradius(lastr, temp[lastr]*tfct, ex))!=0)
              transiterror(TERR_CRITICAL,
                           "computeextradius() return error code %i while "
                           "computing radius #%i: %g\n", rn, r[lastr]*rfct);
            /* Update the value of the extinction at the right place: */
            er[lastr] = e[lastr][wi]*blowex + e_s[lastr] +
                        e_c[lastr] + e_cia[wi][lastr];
          }
        /* FINDME: a while instead of a do-while should work better, huh? */
        }while(bb[ri]*ip->fct < r[lastr]*rfct);
      }
      /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */
      /* Calculate the optical depth and check if tau reached toomuch: */
      /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */
      if( (tau_wn[ri] = rfct * 
           fcn(bb[ri]*riw, r+lastr, n+lastr, er+lastr, rnn-lastr, taulevel))
          > tau.toomuch){
        /* Set tau.last if it reached toomuch: */
        tau.last[wi] = ri;
        if (ri<3){
          transitprint(1, verblevel,
                       "WARNING: At wavenumber %g (cm-1), the critical TAU "
                       "value (%g) was exceeded with tau=%g at the impact "
                       "parameter level %li (%g km), this should have "
                       "happened at a deeper layer (check IP sampling or ATM "
                       "file).\n", wn->v[wi], tau.toomuch, tau_wn[ri],
                       ri, bb[ri]*rfct/1e5);
        }
        /* Exit impact-parameter loop if it reached toomuch: */
        break;
      }
      transitDEBUG(22, verblevel,
                   "Tau(lambda %li=%9.07g, r=%9.4g) : %g  (toomuch: %g)\n",
                   wi, wn->v[wi], r[ri], tau_wn[ri], tau.toomuch);
    }

    if(ri==inn){
      transitprint(1, verblevel,
                   "WARNING: At wavenumber %g cm-1, the bottom of the "
                   "atmosphere was reached before obtaining the critical TAU "
                   "value of %g.\nMaximum TAU reached: %g.\n",
                   wn->v[wi], tau.toomuch, tau_wn[ri]);
      tau.last[wi] = ri-1;
    }

  }

  transitprint(1, verblevel, " Done.\nOptical depth calculated up to %g.\n",
               tr->ds.tau->toomuch);

  /* Print detailed output if requested: */
  if(tr->ds.det->tau.n)
    detailout(&tr->wns, &tr->ips,  &tr->ds.det->tau, tau.t, 0);
  if(tr->ds.det->ext.n)
    detailout(&tr->wns, &tr->rads, &tr->ds.det->ext, e, CIA_RADFIRST);
  if(tr->ds.det->cia.n)
    detailout(&tr->wns, &tr->rads, &tr->ds.det->cia, (double **)e_cia,
              CIA_DOFLOAT);

  if(tr->save.ext)
    savefile_extinct(tr->save.ext, e, comp, rnn, wnn);

  /* Print lowest impact parameter before optical depth gets too big: */
  if(tr->f_toomuch)
    printtoomuch(tr->f_toomuch, tr->ds.tau, &tr->wns, &tr->ips);

  /* Free memory that is no longer needed: */
  //freemem_lineinfotrans(tr->ds.li, &tr->pi);
  freemem_localextinction();

  /* Set progress indicator and output tau if requested, otherwise return
     success: */
  tr->pi |= TRPI_TAU;
  if(tr->fl & TRU_OUTTAU)
    printtau(tr);
  return 0;
}


/* \fcnfh
   Print to file the optical depth, cia, or extinction at the requested
   wavenumbers

   Return: 0 on success                  */
int
detailout(prop_samp *wn,         /* transit's wavenumber array */
          prop_samp *rad,        /* Radius array               */
          struct detailfld *det, /* Detail field struct        */
          PREC_RES **arr,        /* Array of values to store   */
          short flag){           /* Flags                      */
  long i,       /* Auxiliary for-loop index */
       u, d, m; /* Auxiliary binary search indices */

  /* Handle flags: */
  _Bool radfirst = (_Bool)(flag & CIA_RADFIRST); /* rad index is first in arr */
  _Bool dofloat  = (_Bool)(flag & CIA_DOFLOAT);  /* Print as float value      */

  long idx[det->n]; /* Wavenumber indices */
  double val;
  float **arrf = (float **)arr;      /* Float-casted array */
  FILE *out = fopen(det->file, "w"); /* Pointer to file    */
  if(!out)
    transiterror(TERR_SERIOUS, "Cannot open '%s' for writing fine detail.\n",
                               det->file);

  transitprint(1, verblevel, "\nPrinting in '%s'. Fine detail of %s at "
                             "selected wavenumbers.\n", det->file, det->name);

  fprintf(out, "#Radius-w=>    ");
  /* Binary search to find the index for the requested wavenumbers: */
  for(i=0; i<det->n; i++){
    val = det->ref[i];
    u = wn->n-1;
    if(val == wn->v[u])
      d = u;
    else{
      d = 0;
      while(u-d > 1){
        m = (u+d)/2;
        if(wn->v[m] > val)
          u = m;
        else
          d = m;
      }
    }
    idx[i] = d; /* Wavenumber index in transit array */
    /* Print the wavenumber */
    fprintf(out, "%-15.8g", wn->v[idx[i]]);
  }
  fprintf(out, "\n");

  /* Print radii and value: */
  if(radfirst){
    for(m=0; m<rad->n; m++){
      fprintf(out, "%-15.7g", rad->v[m]);
      for(i=0; i<det->n; i++){
        if(dofloat)
          val = arrf[m][idx[i]];
        else
          val = arr[m][idx[i]];
        fprintf(out, "%-15.7g", val);
      }
      fprintf(out, "\n");
    }
  }
  else{
    for(m=0; m<rad->n; m++){
      fprintf(out, "%-15.7g", rad->v[m]);
      for(i=0; i<det->n; i++){
        if(dofloat)
          val = arrf[idx[i]][m];
        else
          val = arr[idx[i]][m];
        fprintf(out, "%-15.7g", val);
      }
      fprintf(out, "\n");
    }
  }

  fclose(out);
  return 0;
}


/* \fcnfh
   Print (to file or stdout) the impact parameter where the optical depth
   reached toomuch (for each wavenumber)  */
void
printtoomuch(char *file,           /* Filename to save the info */
             struct optdepth *tau, /* Tau information           */
             prop_samp *wn,        /* Wavenumber sampling       */
             prop_samp *rad){      /* Radius sampling           */

  long w;             /* Auxiliary for-loop index for wavenumber */
  FILE *out = stdout; /* File pointer */

  /* Open file if it was specified: */
  if(file && file[0] != '-')
    out = fopen(file, "w");
  if(!out)
    transiterror(TERR_WARNING,
                 "Cannot open '%s' for writing maximum depth before "
                 "reaching toomuch optical depth.\n",
                 out==stdout?"STDOUT":file);

  transitprint(1, verblevel,
               "\nPrinting in '%s'.\n"
               "Maximum depth before optical depth got larger than %g, and "
               "therefore impact parameter was not calculated for deeper "
               "layers.\n\n", file, tau->toomuch);

  /* Print header: */
  fprintf(out, "#Wvn[cm-1]  Maximum_calculated_depth\n");
  fprintf(out, "#Wavenumber (cm-1)  Radius at max. calculated depth (cm)\n");
  /* Print the wavenumber and impact parameter: */
  for(w=0; w<wn->n; w++)
    fprintf(out, "%-14.10g%16.12g\n", wn->v[w]*wn->fct,
                                      rad->v[tau->last[w]]*rad->fct);
  fclose(out);
}


/* \fcnfh
   Print (to file or stdout) the optical depth at an specified radius
   (asked interactively)                                                */
void
printtau(struct transit *tr){
  int rn;
  FILE *out = stdout;
  prop_samp *rads  = &tr->ips;
  PREC_RES **t     = tr->ds.tau->t;
  long *last       = tr->ds.tau->last;
  PREC_RES toomuch = tr->ds.tau->toomuch;
  long rad;

  transitcheckcalled(tr->pi, "printtau", 1, "tau", TRPI_TAU);
  tr->ot = tr->ds.th->ot;

  /* Open file if it was specified: */
  if(tr->f_out && tr->f_out[0] != '-')
    out = fopen(tr->f_out, "w");
  if(!out)
    transiterror(TERR_WARNING,
                 "Cannot open '%s' for writing optical depth.\n",
                 out==stdout?"STDOUT":tr->f_out);

  if(tr->ot < 0){
    rad = askforposl("Radius at which you want to print the optical "
                     "depth (%li - %li): ", 1, rads->n) - 1;
    if(rad > rads->n){
      fprintf(stderr, "Value out of range, try again.\n");
      printtau(tr);
    }
  }
  else
    rad=tr->ot;

  transitprint(1, verblevel,
               "\nPrinting in '%s'.\n"
               "Optical depth for radius %li (at %g cm)\n",
               tr->f_out?tr->f_out:"standard output",
               rad+1, rads->fct*rads->v[rad]);

  transitprint(2, verblevel,
               "Optical depth calculated up to %g cm-1.\n", toomuch);

  /* Print the wavenumber, wavelength, and optical depth: */
  fprintf(out, "#Wavenumber [cm-1]\tWavelength [nm]\tOptical depth [cm-1]\n");
  for(rn=0; rn<tr->wns.n; rn++)
    fprintf(out, "%12.6f%14.6f%17.7g\n", tr->wns.fct*tr->wns.v[rn],
                 1/tr->wavs.fct/tr->wns.v[rn]/tr->wns.fct,
                 rad>last[rn]?toomuch:t[rn][rad]);

  exit(EXIT_SUCCESS);
}


/* \fcnfh
   Free tau structure

   Return: 0 on success */
int
freemem_tau(struct optdepth *tau, /* Optical depth structure */
            long *pi){            /* Progress flag           */
  /* Free arrays: */
  free(tau->t[0]);
  free(tau->t);
  free(tau->last);

  /* Update progress indicator and return: */
  *pi &= ~(TRPI_TAU);
  return 0;

}


void
outdebtauex(char *name,
            PREC_RES **e,
            prop_samp *ip,
            PREC_RES **t,
            long rn,
            long w){
  FILE *fp = fopen(name, "w");

  long j;
  for(j=0; j<rn; j++)
    fprintf(fp, "%-15.10g%-15.10g\t%-15.10g\n", ip->v[j], t[w][rn-j-1],
            e[j][w]);
  fclose(fp);
}


void
outdebex(char *name,
         PREC_RES **e,
         PREC_RES *r,
         long rn,
         long wi,
         long wf){
  FILE *fp = fopen(name, "w");

  long i, j;
  for(j=0; j<rn; j++){
    fprintf(fp, "%-15.10g\t", r[j]);
    for(i=wi; i<=wf; i++)
      fprintf(fp, "%-15.10g\t", e[j][i]);
    fprintf(fp, "\n");
  }
  fclose(fp);
}
         

void
outdebtau(char *name,
       prop_samp *ip,
       PREC_RES **t,
       long wi,
       long wf){
  FILE *fp = fopen(name, "w");

  long i, j;
  for(j=0; j<ip->n; j++){
    fprintf(fp, "%-15.10g\t", ip->v[j]);
    for(i=wi; i<=wf; i++)
      fprintf(fp, "%-15.10g\t", t[i][j]);
    fprintf(fp, "\n");
  }
  fclose(fp);
}
