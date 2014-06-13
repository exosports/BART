/*
 * transit.c   - Models the modulation produced by a Planet crossing in
 *               front of the star. Main component of the Transit program.
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

/* TBD: calloc checks */

#include <transit.h>

/* \fcnfh  */
int main(int argc,      /* Number of variables */
         char **argv){  /* Variables           */

  /* FINDME: What is s_lt? */
  /* Initialization of data structure's pointers. Note that s_lt is not
     assigned because the array that is going to point has not been
     initialized yet.  */
  struct transit transit;
  long itr=0;
  struct timeval tv;
  double t0=0.0;

  memset(&transit, 0, sizeof(struct transit));
  verblevel=2;

  /* Process the command line arguments:         */
  fw(processparameters, !=0, argc, argv, &transit);
  t0 = timecheck(verblevel, itr,  0, "processparameters", tv, t0);

  /* Accept all general hints:                   */
  fw(acceptgenhints, !=0, &transit);
  /* Presentation:                               */
  printintro();

  /* No program warnings if verblevel is 0 or 1: */
  if(verblevel<2)
    transit_nowarn = 1;

  /* Make wavenumber binning:                    */
  fw(makewnsample0, <0, &transit);
  t0 = timecheck(verblevel, itr,  1, "makewnsample0", tv, t0);
  if(fw_status>0)
    transitprint(7, verblevel,
                 "makewnsample() modified some of the hinted "
                 "parameters according to returned flag: 0x%lx.\n",
                 fw_status);

  /* Read Atmosphere information:                */
  fw(getatm, !=0, &transit);
  t0 = timecheck(verblevel, itr,  2, "getatm", tv, t0);

  /* Read line info:                             */
  fw(readlineinfo, !=0, &transit);
  t0 = timecheck(verblevel, itr,  3, "readlineinfo", tv, t0);

  /* Hack: add an index to f_out filename to get different files
     and compare them:                                               */
  int fout_len = (int)strlen(transit.f_out);  /* Length of tr.f_out  */
  char *dot    = strchr(transit.f_out, '.');  /* Search for '.' char */
  int dot_pos  = fout_len - (int)strlen(dot); /* Position of dot     */
  char fout[fout_len+1];         /* Copy of tr.f_out           */     
  char str_iter[1];              /* String of iteration number */
  strcpy(fout, transit.f_out);
  strcpy(fout+dot_pos+1, dot);
  strncpy(fout+dot_pos, "0", 1);

  /* Make radius binning and interpolate data to new value: */
  fw(makeradsample, <0, &transit);
  t0 = timecheck(verblevel, itr,  4, "makeradsample", tv, t0);
  if(fw_status>0)
    transitprint(7, verblevel,
                 "makeradsample() modified some of the hinted "
                 "parameters according to returned flag: 0x%lx.\n",
                 fw_status);

  /* Compute sampling of impact parameter: */
  fw(makeipsample, <0, &transit);
  t0 = timecheck(verblevel, itr,  5, "makeipsample", tv, t0);
  if(fw_status>0)
    transitprint(7, verblevel,
                 "makeipsample() modified some of the hinted "
                 "parameters according to returned flag: 0x%lx.\n",
                 fw_status);
 
  /* Print sampling info:                  */
  fw(outsample, !=0, &transit);
  t0 = timecheck(verblevel, itr,  6, "outsample", tv, t0);

  /* EDIT: The loop should enclose getatm (getatm will change in the future,
     but for the moment we will leave it as it is).  */
  for (itr=0; itr<1; itr++){
    t0 = timecheck(verblevel, itr,  7, "Start loop", tv, t0);

    /* Initialize CIA:                       */
    fw(interpolatecia, !=0, &transit);
    t0 = timecheck(verblevel, itr,  8, "interpolatecia", tv, t0);
 
    /* Compute index of refraction:          */
    fw(idxrefrac, !=0, &transit);
    t0 = timecheck(verblevel, itr,  9, "idxrefrac", tv, t0);
 
    /* Calculate extinction coefficient:     */
    fw(extwn, !=0, &transit);
    t0 = timecheck(verblevel, itr, 10, "extwn", tv, t0);
 
    sprintf(str_iter, "%li", itr);
    strncpy(fout+dot_pos, str_iter, 1);
    strcpy(transit.f_out, fout);

    /* Ray solutions choice:                  */
    RaySol path=transit.ds.th->path;

    /* Calculates optical depth for eclipse   */
    if(path == eclipse){
      transitprint(1,verblevel, "\nCalculating eclipse.\n\n");

      /* Angle number                        */
      struct transithint *th = transit.ds.th;
      long int an = th->ann;

      /* Sets intensity grid:                */
      fw(intens_grid, !=0, &transit);
      for(int i = 0; i < an; i++){
        /* Fills out angle index             */
        transit.angleIndex = i;

        fw(tau_eclipse, !=0, &transit);
        t0 = timecheck(verblevel, itr, 11, "tau eclipse", tv, t0);
  
        /* Calculates eclipse intensity:          */
        /* In cgs units erg/s/sr/cm               */
        fw(emergent_intens, !=0, &transit);
        t0 = timecheck(verblevel, itr, 12, "emergent intensity", tv, t0);
      }

    /* Calculates flux  erg/s/cm               */
    fw(flux, !=0, &transit);
    t0 = timecheck(verblevel, itr, 13, "flux", tv, t0);

    /* Free no longer needed memory JB added   */
    freemem_intensityGrid(transit.ds.intens, &transit.pi);
    }

    /* Calculate optical depth for transit:    */
    else{
      transitprint(1,verblevel, "\nCalculating transit.\n");
      fw(tau, !=0, &transit);
      t0 = timecheck(verblevel, itr, 11, "tau", tv, t0); 

     /* Calculates transit modulation:         */
      fw(modulation, !=0, &transit);
      t0 = timecheck(verblevel, itr, 12, "modulation", tv, t0);
   }
 
    /* Free no longer needed memory            */
    freemem_idexrefrac(transit.ds.ir,        &transit.pi);
    freemem_extinction(transit.ds.ex,        &transit.pi);
    freemem_tau(transit.ds.tau,              &transit.pi);

    free(transit.save.ext);
    freemem_cia      (transit.ds.cia, &transit.pi);
    freemem_outputray(transit.ds.out, &transit.pi);
    t0 = timecheck(verblevel, itr, 13, "THE END", tv, t0);
    transitprint(1, verblevel, "----------------------------\n");
  }
  freemem_isotopes(     transit.ds.iso, &transit.pi);
  freemem_molecules(    transit.ds.mol, &transit.pi);
  freemem_atmosphere(   transit.ds.at,  &transit.pi);
  freemem_lineinfotrans(transit.ds.li,  &transit.pi);
  freemem_transit(&transit);

  return EXIT_SUCCESS;
}


/* \fcnfh
   Frees transit structure.                 */
void
freemem_transit(struct transit *tr){
  freemem_hints(tr->ds.th);

  freemem_samp(&tr->rads);
  freemem_samp(&tr->wns);
  freemem_samp(&tr->ips);
  free_atm(&tr->atm);

  free(tr->outpret);
  /* TBD: Free saves once it is enabled
  freemem_saves();                          */
}
