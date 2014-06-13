/*
 * observable.c   - Finds the optical depth. Component of the Transit program.
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

/* \fcnfh
   Calculate the transit modulation at each wavenumber

   Return: 0 on success, else
          -1 if impact parameter sampling is not equispaced  */
int
modulation(struct transit *tr){ /* Main structure */
  static struct outputray st_out;
  tr->ds.out = &st_out;

  /* Check that impact parameter and wavenumber samples exist: */
  transitcheckcalled(tr->pi, "modulation", 3, "tau",          TRPI_TAU,
                                              "makeipsample", TRPI_MAKEIP,
                                              "makewnsample", TRPI_MAKEWN);

  /* Initial variables:  */
  long w;
  prop_samp *ip = &tr->ips;
  prop_samp *wn = &tr->wns;
  transit_ray_solution *sol = tr->sol;
  /* Check for monospaced impact parameter if required (monoip): */
  if(ip->d==0 && sol->monoip){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "To compute %s modulation, the impact parameter has to "
                 "be an equispaced array.\n", sol->name);
    return -1;
  }

  /* Allocate the modulation array: */
  PREC_RES *out = st_out.o = (PREC_RES *)calloc(wn->n, sizeof(PREC_RES));
  struct geometry *sg  = tr->ds.sg;
  struct optdepth *tau = tr->ds.tau;

  /* Set time to the user hinted default, and other user hints: */
  setgeom(sg, HUGE_VAL, &tr->pi);
  const int modlevel = tr->modlevel = tr->ds.th->modlevel;

  /* Integrate for each wavelength: */
  transitprint(1, verblevel, "\nIntegrating over wavelength.\n");

  int nextw = wn->n/10;

  /* Calculate the modulation at each wavenumber: */
  for(w=0; w<wn->n; w++){
    out[w] = sol->obsperwn(tau->t[w], tau->last[w], tau->toomuch, ip,
                           sg, modlevel);
    if(out[w]<0){
      switch(-(int)out[w]){
      case 1:
        if(modlevel == -1)
          transiterror(TERR_SERIOUS,
                       "Optical depth didn't reach limiting %g at wavenumber "
                       "%g cm-1 (only reached %g). Cannot use critical "
                       "radius technique (-1).\n", tau->toomuch,
                       tau->t[w][tau->last[w]], wn->v[w]*wn->fct);
      default:
        transiterror(TERR_SERIOUS,
                     "There was a problem while calculating modulation "
                     "at wavenumber %g cm-1. Error code %i.\n",
                     wn->v[w]*wn->fct, (int)out[w]);
        break;
      }
      exit(EXIT_FAILURE);
    }

    /* Print to screen the progress status: */
    if(w==nextw){
      nextw += wn->n/10;
      transitprint(2, verblevel, "%i%% ", (10*(int)(10*w/wn->n+0.9999999999)));
    }
  }
  transitprint(1, verblevel, "\nDone.\n");

  /* Set progress indicator, and print output: */
  tr->pi |= TRPI_MODULATION;
  printmod(tr);  
  return 0;
}


/* \fcnfh
   Print (to file or stdout) the modulation as function of wavelength */
void
printmod(struct transit *tr){
  FILE *outf = stdout;
  struct outputray *outray = tr->ds.out;
  int rn;

  /* Open file: */
  if(tr->f_out && tr->f_out[0] != '-')
    outf = fopen(tr->f_out, "w");

  transitprint(1, verblevel,
               "\nPrinting in-transit/out-transit modulation in '%s'.\n",
               tr->f_out?tr->f_out:"standard output");

  /* Print: */
  char wlu[20], /* Wavelength units name */
       wnu[20]; /* Wavenumber units name (the inverse, actually) */
  long nsd = (long)(1e6);

  /* Get wavenumber units name: */
  if((long)(nsd*tr->wns.fct)==nsd) strcpy(wnu, "cm");
  else if((long)(nsd*1e-1*tr->wns.fct)==nsd) strcpy(wnu, "mm");
  else if((long)(nsd*1e-4*tr->wns.fct)==nsd) strcpy(wnu, "um");
  else if((long)(nsd*1e-7*tr->wns.fct)==nsd) strcpy(wnu, "nm");
  else if((long)(nsd*1e-8*tr->wns.fct)==nsd) strcpy(wnu, "A ");
  else sprintf(wnu, "%6.1g cm", 1/tr->wns.fct);

  /* Get wavelength units name: */
  if((long)(nsd*tr->wavs.fct)==nsd) strcpy(wlu, "cm");
  else if((long)(nsd*1e1*tr->wavs.fct)==nsd) strcpy(wlu, "mm");
  else if((long)(nsd*1e4*tr->wavs.fct)==nsd) strcpy(wlu, "um");
  else if((long)(nsd*1e7*tr->wavs.fct)==nsd) strcpy(wlu, "nm");
  else if((long)(nsd*1e8*tr->wavs.fct)==nsd) strcpy(wlu, "A ");
  else sprintf(wlu, "%8.1g cm", tr->wavs.fct);

  /* Print header: */
  fprintf(outf, "#wvn %s-1%*s wvl %s%*s modulation\n",
          wnu, (int)(9-strlen(wnu)), "", wlu, (int)(12-strlen(wlu)),"");
  /* Print wavenumber, wavelength, and modulation at each wavenumber: */
  for(rn=0; rn<tr->wns.n; rn++)
    fprintf(outf, "%-17.9g%-17.9g%-18.9g\n",
            tr->wns.v[rn]/tr->wns.fct,
            1/tr->wavs.fct/tr->wns.v[rn]/tr->wns.fct,
            outray->o[rn]);

  fclose(outf);
  return;
}


/*\fcnfh
  Free the transit modulation array

  Return: 0 on success                          */
int
freemem_outputray(struct outputray *out,
                  long *pi){
  /* Free arrays: */
  free(out->o);

  /* Clear PI and return: */
  *pi &= ~(TRPI_MODULATION);
  return 0;
}
