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
  struct transithint trh;
  memset(&transit, 0, sizeof(struct transit));
  memset(&trh, 0, sizeof(struct transithint));
  transit.ds.th=&trh;

  //Initialization of wavelength sampling in nanometers.
  //Fields from 'trh.wavs'.
  //'i' and 'f' fields are for initial and final wavelength, if 0 then
  //the most extreme common minimum or maximum of the line databases is
  //chosen. 
  //'d' field is the wavelength spacing.
  //'n' is the number of sampled elements: it should always be set to 0
  //whenever 'd' is changed, so that 'makewavsample()' can update
  //correctly.
  //'o' is the oversampling.
  //'v' is the value's array.
  prop_samp *samp=&trh.wavs;
  samp->i=0;
  samp->f=0;
  samp->d=0.188;
  samp->n=0;
  samp->o=100;
  samp->v=NULL;
  samp->fct=0;
  trh.na|=TRH_WAV|TRH_WAVO;

  //Initialization of radius sampling in planetary radius.
  //Fields from 'trh.rads'.
  //'.i', '.f', '.d', '.n', '.v' are equivalent to the wavelegth
  //sampling above but for the radius. If zeroed then '.i', '.f' or '.d'
  //are taken from the atmosphere datafile.
  //'o' equal to 0 is necessary to dissable oversampling in radius.
  samp=&trh.rads;
  samp->i=0;
  samp->f=0;
  samp->d=0.6;
  samp->n=0;
  samp->v=NULL;
  samp->o=0;
  samp->fct=0;
  trh.na|=TRH_RAD;

  //Initialization of wavenumber sampling in cm-1. Note that the
  //extincitons are calculated in wavenumber and then optionally
  //interpolated to wavelength.
  //Fields from 'trh.wns'.
  //'.i', '.f', '.d', '.n', '.v' are equivalent to the wavelegth
  //sampling above but for the radius. If zeroed then '.i', '.f' are
  //taken from the wavelength transformation. '.d' is calculated so that
  //the same number of points as requested for wavelength are outputted.
  //'o' equal to 0 is necessary to dissable oversampling in radius.
  samp=&trh.wns;
  samp->i=5716.5;
  samp->f=5718;
  samp->i=0;
  samp->f=0;
  samp->d=0;
  samp->n=0;
  samp->v=NULL;
  samp->o=0;
  samp->fct=0;
  trh.na|=TRH_WN|TRH_WNO;

  //Initialization of general variables.
  //'rc' and 'rn' are the general auxiliary variables.
  //'file\_out' is the name of the output file, a '-' indicates standard
  //output.
  //'trh.verbnoise' is the noisiest possible verbose level
  //controlled by 'verbose' when in a non-debugging compilation, a value
  //of 10 or more is only used for debugging.
  //'defile\_out' is the default output filename
  int rn;
  trh.verbnoise=4;
  char defile_out[]="-";
  trh.f_out=(char *)calloc(strlen(defile_out)+1,sizeof(char));
  strcpy(trh.f_out,defile_out);
  trh.na|=TRH_FO;
#ifdef NODEBUG_TRANSIT
  verblevel=2;
#else
  verblevel=20;
#endif /* NODEBUG_TRANSIT */

  //Initialization of line database variables. 
  //'f\_line' is the name of the twii data file 
  //'trh.m' is the amount of microns that cannot be trusted at the
  //boundaries of databases range, and also how much extra out of
  //requested range it has to look for transitions.
  //'defile\_line' is the default name of the line info file.
  //(i.e. modifiable by the user)
  trh.m=0.001;
  trh.na|=TRH_WM;
  char defile_line[]="./res/lineread.twii";
  trh.f_line=(char *)calloc(strlen(defile_line)+1,sizeof(char));
  strcpy(trh.f_line,defile_line);
  trh.na|=TRH_FL;


  //Initialization of atmospheric parameters.
  //'.mass' indicates whether the abundances are by mass or number
  //'.f\_atm' is the name of the file with the atmospheric parameters,
  //a '-' indicates that a one-point solution is desired.
  trh.mass=1;
  char defile_atm[]="-";
  trh.f_atm=(char *)calloc(strlen(defile_atm)+1,sizeof(char));
  strcpy(trh.f_atm,defile_atm);
  trh.allowrq=0.01;
  trh.na|=TRH_FA;
  trh.fl|=TRU_ATMASK1P|TRU_SAMPLIN|TRH_MASS;

  //Initialization of extinction parameters.
  //'.voigtfine' is the fine binning of voigt
  trh.voigtfine=5;
  trh.timesalpha=50;
  trh.maxratio_doppler=0.001;
  trh.na|=TRH_VF|TRH_TA|TRH_DR;


  //Initialization of optical depth parameters
  char defsol[]="slant path";
  trh.solname=strdup(defsol);
  trh.tauiso=0;
  trh.toomuch=20;
  trh.na|=TRH_TOOMUCH|TRH_TAUISO|TRH_ST;
  trh.fl|=TRU_OUTTAU;


  //Command line parameters' processing
  if((rn=processparameters(argc,argv,&trh))!=0)
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
  //Print and output extinction if one P,T was desired
  if(transit.rads.n==1)
    printone(&transit);

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


  //Calculates optical depth
  if((rn=tau(&transit))!=0)
    transiterror(TERR_SERIOUS,
		 "tau() returned error code %i\n"
		 ,rn);
  //Output tau if requested
  if(transit.fl&TRU_OUTTAU)
    printtau(&transit);


  return EXIT_SUCCESS;
}


/* \fcnfh
   Printout for optical depth at a requested radius
*/
void
printtau(struct transit *tr)
{
  int rn;
  FILE *out=stdout;
  prop_samp *rads=&tr->ips;
  PREC_RES **t=tr->ds.tau->t;
  long *frst=tr->ds.tau->first;
  PREC_RES toomuch=tr->ds.tau->toomuch;

  transitcheckcalled(tr->pi,"printtau",1,
		     "tau",TRPI_TAU);

  //open file
  if(tr->f_out&&tr->f_out[0]!='-')
    out=fopen(tr->f_out,"w");

  long rad=
    askforposl("Radius at which you want to print the optical depth(%li - %li): "
	       ,1,rads->n)-1;
  if(rad>rads->v[rads->n-1]){
    fprintf(stderr,"Value out of range, try again\n");
    printtau(tr);
  }

  transitprint(1,verblevel,
	       "\nPrinting optical depth for radius %li (at %gcm) in '%s'\n"
	       "Optical depth calculated up to %g[cm-1]\n"
	       ,rad+1,rads->v[rad],tr->f_out?tr->f_out:"standard output",toomuch);

  //print!
  fprintf(out,
	  "#wavenumber[cm-1]\twavelength[nm]\toptical depth[cm-1]\n");
  for(rn=0;rn<tr->wns.n;rn++)
    fprintf(out,"%12.6f%14.6f%17.7g\n"
	    ,tr->wns.fct*tr->wns.v[rn],WNU_O_WLU/tr->wns.v[rn]/tr->wns.fct,
	    rad<frst[rn]?toomuch:t[rn][rad]);

  exit(EXIT_SUCCESS);
}


/* \fcnfh
   Printout for one P,T conditions
*/
void
printone(struct transit *tr)
{
  int rn;
  FILE *out=stdout;

  //open file
  if(tr->f_out&&tr->f_out[0]!='-')
    out=fopen(tr->f_out,"w");

  transitprint(1,verblevel,
	       "\nPrinting extinction for one radius (at %gcm) in '%s'\n"
	       ,tr->rads.v[0],tr->f_out?tr->f_out:"standard output");

  //print!
  fprintf(out,
	  "#wavenumber[cm-1]\twavelength[nm]\textinction[cm-1]\tcross-section[cm2]\n");
  for(rn=0;rn<tr->wns.n;rn++)
    fprintf(out,"%12.6f%14.6f%17.7g%17.7g\n"
	    ,tr->wns.fct*tr->wns.v[rn],WNU_O_WLU/tr->wns.v[rn]/tr->wns.fct,
	    tr->ds.ex->e[0][0][rn],
	    AMU*tr->ds.ex->e[0][0][rn]*tr->isof[0].m/tr->isov[0].d[0]);

  exit(EXIT_SUCCESS);
}

