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

/* TD: Command line parameters parsing */
/* TD: calloc checks */
/* TD: atmosphere info  retrieval */

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

#ifdef NODEBUG_TRANSIT
  verblevel=2;
#else
  verblevel=20;
#endif /* NODEBUG_TRANSIT */

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
  char revname[20];
  if(revision<0) snprintf(revname,20,"pre%i",-revision);
  else snprintf(revname,20,".%i",revision);
  transitprint(1,verblevel,
	       "-----------------------------------------------\n"
	       "                TRANSIT v%i%s\n"
	       "-----------------------------------------------\n"
	       ,version,revname);

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
    transitprint(1,verblevel,
		 "makewavsample() modified some parameters according\n"
		 "to returned flag: 0x%x\n"
		 ,rn);

  //Make wavenumber binning
  if((rn=makewnsample(&transit))<0)
    transiterror(TERR_SERIOUS,
		 "makewnsample() returned error code %i\n"
		 ,rn);
  if(rn>0)
    transitprint(1,verblevel,
		 "makewnsample() modified some parameters according\n"
		 "to returned flag: 0x%x\n"
		 ,rn);

  //Make radius binning and interpolate data to new value
  if((rn=makeradsample(&transit))<0)
    transiterror(TERR_SERIOUS,
		 "makeradsample() returned error code %i\n"
		 ,rn);
  if(rn>0)
    transitprint(1,verblevel,
		 "makeradsample() modified some parameters according\n"
		 "to returned flag: 0x%x\n"
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
    transitprint(1,verblevel,
		 "makeipsample() modified some parameters according\n"
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

  return EXIT_SUCCESS;
}
