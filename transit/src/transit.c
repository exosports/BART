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

/* TD: calloc checks */

#include <transit.h>

//\fcnfh
int main (int argc,		/* Number of variables */
	  char **argv)		/* Variables*/
{

  //Initialization of data structure's pointers. Note that s\_lt is not
  //assigned because the array that is going to point has not been
  //initialized yet.
  struct transit transit;
  memset(&transit, 0, sizeof(struct transit));

  verblevel=2;

  //Command line parameters' processing
  fw(processparameters, !=0, argc, argv, &transit);


  //Accept all general hints
  fw(acceptgenhints, !=0, &transit);


  //Presentation
  printintro();

  //No program warnings if verblevel is 0 or 1
  if(verblevel<2)
    transit_nowarn=1;


  //Read line info
  fw(readlineinfo, !=0, &transit);


  //Read Atmosphere information
  fw(getatm, !=0, &transit);


  //Make wavelength binning
  fw(makewavsample, <0, &transit);
  if(fw_status>0)
    transitprint(7,verblevel,
		 "makewavsample() modified some of the hinted\n"
		 " parameters according to returned flag: 0x%lx\n"
		 ,fw_status);


  //Make wavenumber binning
  fw(makewnsample, <0, &transit);
  if(fw_status>0)
    transitprint(7,verblevel,
		 "makewnsample() modified some of the hinted\n"
		 " parameters according to returned flag: 0x%lx\n"
		 ,fw_status);


  //Make radius binning and interpolate data to new value
  fw(makeradsample, <0, &transit);
  if(fw_status>0)
    transitprint(7,verblevel,
		 "makeradsample() modified some of the hinted\n"
		 " parameters according to returned flag: 0x%lx\n"
		 ,fw_status);


  //Initializes CIA
  fw(interpolatecia, !=0, &transit);


  //Computes index of refraction
  fw(idxrefrac, !=0, &transit);


  //Calculates extinction coefficient
  fw(extwn, !=0, &transit);


  //Computes sampling of impact parameter
  fw(makeipsample, <0, &transit);
  if(fw_status>0)
    transitprint(7,verblevel,
		 "makeipsample() modified some of the hinted\n"
		 " parameters according to returned flag: 0x%lx\n"
		 ,fw_status);


  //Prints sampling info
  fw(outsample, !=0, &transit);


  //Calculates optical depth
  fw(tau, !=0, &transit);


  //Calculates eclipse modulation
  fw(modulation, !=0, &transit);


  free(transit.save.ext);
  freemem_isotopes (transit.ds.iso, &transit.pi);
  freemem_cia      (transit.ds.cia, &transit.pi);
  freemem_outputray(transit.ds.out, &transit.pi);
  freemem_transit(&transit);

  return EXIT_SUCCESS;
}


/* \fcnfh
   Frees transit structure
*/
void
freemem_transit(struct transit *tr)
{
  freemem_hints(tr->ds.th);

  freemem_samp(&tr->rads);
  freemem_samp(&tr->wavs);
  freemem_samp(&tr->wns);
  freemem_samp(&tr->ips);
  free_atm(&tr->atm);

  free(tr->outpret);
  /* TD:Free saves once it is enabled
     freemem_saves(); */


}
