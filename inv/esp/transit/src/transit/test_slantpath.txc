/*
 * test_slantpath.c
 * test_slantpath.txc - Tests to check correct working of slantpath.c.
 *                      Component of the transit program
 *
 * Copyright (C) 2004 Patricio Rojo (pato@astro.cornell.edu)
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

#include <transit.h>
#include "slantpath.c"
#include <test_common.h>

const transit_ray_solution *sol=&slantpath;
/* Following is number of cases to test for different extinctions:
   constant, incresing outwards, increasing inwards */
#define NWN  3
const char *ktypes[NWN]={"constant",
			 "increasing outwards",
			 "increasing inwards"};

/* \fcnfh
   Call totaltau and compare value returned with what is expected
   analytically 
*/
void
test_tau_monorad_norefr_now(PREC_RES *b, /* differential impact
					    parameter with respect to
					    maximum value */
			    PREC_RES *rad, /* radius array */
			    PREC_RES *refr, /* refractivity index */
			    PREC_RES **ex, /* extinction[rad][nwn] */
			    long nrad, /* number of radii elements */
			    long wn, /* wavenumber looked */
			    PREC_RES *res,
			    int nimp,
			    char **ipdesc,
			    int status)
{
  double acceptrelerror=.0001;
  gsl_interp_accel acc={0,0,0};
  int i;
  PREC_RES result;

  for(i=0;i<nimp;i++){
    result=totaltau(b[i],rad,refr,ex,nrad,wn,&acc);
    test_result("While %s extinction, at %s(%.7g),\n"
		"  observed %.10g, expect %.10g, error %g\n"
		,ktypes[wn],ipdesc[i],b[i],result,res[i],
		fabs(result-res[i])/res[i]);
    if(fabs(result-res[i])/res[i]>acceptrelerror)
      test_fail(status,"Error bigger than %.5g\n"
		,acceptrelerror);
    else
      test_succeed();
  }

}


/* \fcnfh
   check monospaced radius

   @returns number of errors
*/
int
test_tau_monorad_constrefr(int nrad,
			   double alpha)
{
  PREC_RES refr[nrad], rad[nrad],*ex[nrad],res[nrad];
  int i,status=0;
  ex[0]=(PREC_RES *)calloc(nrad*NWN,sizeof(PREC_RES));
  if(!ex[0]){
    fprintf(stderr,
	    "Unable to allowcate memory!. I required %i bytes in\n"
	    "line %i of file %s. ABORTING TEST...\n"
	    ,nrad*NWN*sizeof(PREC_RES),__LINE__,__FILE__);
    exit(EXIT_FAILURE);
  }

  //Print filling of arrays if necessary to debug
  if(0)
    fprintf(stderr,
	    "radius        refr          ex_cons       ex_out        ex_in\n");
  PREC_RES rm=nrad*1.0;
  for(i=0;i<nrad;i++){
    //check for equispaced radius starting from 1.0
    rad[i]=1.0+i;
    //check for nonbending rays
    refr[i]=1.0;
    //Fill different cases of extinction
    ex[i]=ex[0]+NWN*i;
    //extinction: constant
    ex[i][0]=alpha;
    //extinction: incresing outward
    ex[i][1]=alpha*rad[i];
    //extinction: incresing inward
    ex[i][2]=alpha*(rm-rad[i]);

    //Print filling of arrays if necessary to debug
    if(0)
      fprintf(stderr,"%-14.9g%-14.9g%-14.9g%-14.9g%-14.9g\n"
	      ,rad[i],refr[i],ex[i][0],ex[i][1],ex[i][2]);
  }


  //tests at impact parameters: grazing, interpolated-radius, exact-radius,
  //lowest-radius 
#define nimpact 4
  PREC_RES b[nimpact]={rad[nrad-2], (rad[(3*nrad)/4-1] + rad[(3*nrad)/4])/2.0,
		       rad[nrad/2-1], rad[0]};
  char *ipdesc[nimpact]={"grazing IP","interpolated IP","coincident IP",
			 "lowest IP"};

  //Test the atmosphere with constant extinction at different radius
  //$res=2 ex \sqrt{rm^2-b^2}$
  /*    PREC_RES res[nimpact]={8.71779788708134710448,
	13.2287565553229529525,
	17.32050807568877293528,
	19.89974874213239909468};*/
  for(i=0;i<nrad;i++)
    res[i]=2*alpha*sqrt(rm*rm-b[i]*b[i]);

  test_tau_monorad_norefr_now(b,rad,refr,ex,nrad,0
			      ,res,nimpact,ipdesc,status);

  //Test the atmosphere with increasing extinction outwards
  //$res=\alpha ( rm\sqrt{rm^2-b^2} + b^2 \log{\sqrt{(rm/b)^2-1}+rm/b})$
  // Calc function:\par
  //{\em define
  //res(b)=alpha*(rm*sqrt(rm*rm-b*b)+b*b*ln(sqrt((rm/b)\^2-1)+rm/b)/ln(10))}
  //\par $ex=\alpha r$
  /*    PREC_RES res[nimpact]={60.02215842946226671978,
			   85.57381701507597253714,
			   100.90122906677784957218,
			   100.79868387584142893220};*/
  for(i=0;i<nrad;i++)
    res[i]=alpha*(rm*sqrt(rm*rm-b[i]*b[i]) +
		  b[i]*b[i]*log(sqrt(rm*rm/b[i]/b[i]-1)+rm/b[i])/log(10.0));

  test_tau_monorad_norefr_now(b,rad,refr,ex,nrad,1
			      ,res,nimpact,ipdesc,status);

  //Test the atmosphere with increasing extinction inwards
  //Test the atmosphere with increasing extinction outwards
  //$res=\alpha ( rm\sqrt{rm^2-b^2} - b^2 \log{\sqrt{(rm/b)^2-1}+rm/b})$
  // Calc function:\par
  //{\em define
  //res(b)=alpha*(rm*sqrt(rm*rm-b*b)-b*b*ln(sqrt((rm/b)\^2-1)+rm/b)/ln(10))}
  //\par $ex=\alpha (rm-r)$
  /*    PREC_RES res[nimpact]={27.15582044135120432502,
			   46.71374853815355698786,
			   72.30385169010987978062,
			   98.19880354548256201460};*/
  for(i=0;i<nrad;i++)
    res[i]=alpha*(rm*sqrt(rm*rm-b[i]*b[i]) -
		  b[i]*b[i]*log(sqrt(rm*rm/b[i]/b[i]-1)+rm/b[i])/log(10.0));

  test_tau_monorad_norefr_now(b,rad,refr,ex,nrad,2
			      ,res,nimpact,ipdesc,status);

  return status;
}


/* \fcnfh
   check totaltau funtion

   @returns number of errors
*/
int test_tau()
{
  int status = 0;

  //First check for monospaced radius with only 10 samples, with alpha=1
  status += test_tau_monorad_constrefr(10, 1.0);
  

  /* TD: test for nonmonospaced radius */
  /* TD: test for bent rays */

  return status;
}



int 
main(int argc, char *argv[])
{
  int status = 0;

  status += test_tau( );

  if(status){
    fprintf(stderr,"\nFail: %i error%s found\n",status,status>1?"s":"");
    exit(EXIT_FAILURE);
  }
  exit(EXIT_SUCCESS);
}
