/*
 * transit.c
 * transit.txc - Models the modulation produced by a Planet crossing in
 *               front of the star. Main component of the Transit program.
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
  int rn;

  verblevel=2;

  //Command line parameters' processing
  if((rn=processparameters(argc,argv,&transit))!=0)
    transiterror(TERR_SERIOUS,
		 "processparameters() returned error code %i\n"
		 ,rn);

  //Accept all general hints
  if((rn=acceptgenhints(&transit))!=0)
    transiterror(TERR_SERIOUS,
		 "transitaccepthints() returned error code %i\n"
		 ,rn);

  //Presentation
  printintro();

  //No program warnings if verblevel is 0 or 1
  if(verblevel<2)
    transit_nowarn=1;

  //Read line info
  if((rn=readlineinfo(&transit))!=0)
    transiterror(TERR_SERIOUS,
		 "readlineinfo() returned error code %i\n"
		 ,rn);

  //Read Atmosphere information
  if((rn=getatm(&transit))!=0)  
    transiterror(TERR_SERIOUS,
		 "getatm() returned error code %i\n"
		 ,rn);


  //Make wavelength binning
  if((rn=makewavsample(&transit))<0)
    transiterror(TERR_SERIOUS,
		 "makewavsample() returned error code %i\n"
		 ,rn);
  if(rn>0)
    transitprint(7,verblevel,
		 "makewavsample() modified some of the hinted parameters according\n"
		 "to returned flag: 0x%x\n"
		 ,rn);

  //Make wavenumber binning
  if((rn=makewnsample(&transit))<0)
    transiterror(TERR_SERIOUS,
		 "makewnsample() returned error code %i\n"
		 ,rn);
  if(rn>0)
    transitprint(7,verblevel,
		 "makewnsample() modified some of the hinted parameters according\n"
		 "to returned flag: 0x%x\n"
		 ,rn);

  //Make radius binning and interpolate data to new value
  if((rn=makeradsample(&transit))<0)
    transiterror(TERR_SERIOUS,
		 "makeradsample() returned error code %i\n"
		 ,rn);
  if(rn>0)
    transitprint(7,verblevel,
		 "makeradsample() modified some of the hinted parameters according\n"
		 "to returned flag: 0x%x\n"
		 ,rn);

  //Initializes CIA
  if((rn=interpolatecia(&transit))!=0)
    transiterror(TERR_SERIOUS,
		 "readcia() returned error code %i\n"
		 ,rn);

  //Computes index of refraction
  if((rn=idxrefrac(&transit))!=0)
    transiterror(TERR_SERIOUS,
		 "idxrefrac() returned error code %i\n"
		 ,rn);

  //Calculates extinction coefficient
  if((rn=extwn(&transit))!=0)
    transiterror(TERR_SERIOUS,
		 "extwn() returned error code %i\n"
		 ,rn);

  //Computes sampling of impact parameter
  if((rn=makeipsample(&transit))<0)
    transiterror(TERR_SERIOUS,
		 "makeipsample() returned error code %i\n"
		 ,rn);
  if(rn>0)
    transitprint(7,verblevel,
		 "makeipsample() modified some of the hinted parameters according\n"
		 "to returned flag: 0x%x\n"
		 ,rn);

  //Prints sampling info
  if((rn=outsample(&transit))<0)
    transiterror(TERR_SERIOUS,
		 "outsample() returned error code %i\n"
		 ,rn);

  //Calculates optical depth
  if((rn=tau(&transit))!=0)
    transiterror(TERR_SERIOUS,
		 "tau() returned error code %i\n"
		 ,rn);

  //Calculates eclipse modulation
  if((rn=modulation(&transit))!=0)
    transiterror(TERR_SERIOUS,
		 "modulation() returned error code %i\n"
		 ,rn);

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
