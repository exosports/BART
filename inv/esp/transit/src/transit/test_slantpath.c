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
long ntests=0;
FILE *outp=NULL;

/* \fcnfh
   Call totaltau and compare value returned with what is expected
   analytically 
*/
int
test_tau_monorad_norefr_now(PREC_RES *b, /* differential impact
					    parameter with respect to
					    maximum value */
			    PREC_RES *rad, /* radius array */
			    PREC_RES *refr, /* refractivity index */
			    PREC_RES **exi, /* extinction[wn][rad] */
			    long nrad, /* number of radii elements */
			    long wn, /* wavenumber looked */
			    PREC_RES *res, /* expected result */
			    int nimp, /* Number of impact parameters */
			    char **ipdesc, /* Name of impact parameters
					      */
			    double acceptrelerror)
{
  long i,status=0;
  PREC_RES result;
  PREC_RES *ex=exi[wn];

  test_result("** extinction %s **\n"
	      ,ktypes[wn]);
  for(i=0;i<nimp;i++){
    test_result(" => IP %s(%.3g). Expected result %g\n",ipdesc[i],b[i],res[i]);
    //test first level of precision (no ray bending)
    result=totaltau(b[i],rad,refr,ex,nrad,1);
    test_result("     RefIdx constant technique: observed %13.10g (error %g)\n"
		,result,fabs(result-res[i])/res[i]);
    if(fabs(result-res[i])/res[i]>acceptrelerror)
      status++;
    ntests++;

    //test with ray bending
    result=totaltau(b[i],rad,refr,ex,nrad,2);
    test_result("     Ray bending technique    : observed %13.10g (error %g)\n"
		,result,fabs(result-res[i])/res[i]);
    if(fabs(result-res[i])/res[i]>acceptrelerror)
      status++;
    ntests++;
  }

  //return number of fails
  return status;
}


/* \fcnfh
   check monospaced radius

   @returns number of errors
*/
int
test_tau_monorad_constrefr(int nrad,
			   double alpha,
			   double spacing)
{

  test_result("#############################\n"
	      "Testing totaltau() with alpha=%g, #elements=%i, constant refractivity index\n"
	      ,alpha,nrad);
  int ntb=ntests;
  double acceptrelerror=.0001;
  PREC_RES refr[nrad], rad[nrad],*ex[NWN];
  long i,status=0;
  ex[0]=(PREC_RES *)calloc(nrad*NWN,sizeof(PREC_RES));
  //Check that memory was really allocated and fill different cases of
  //extinction 
  if(!ex[0]){
    fprintf(stderr,
	    "Unable to allowcate memory!. I required %i bytes in\n"
	    "line %i of file %s. ABORTING TEST...\n"
	    ,nrad*NWN*sizeof(PREC_RES),__LINE__,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<NWN;i++)
    ex[i]=ex[0]+i*nrad;

  //tests at impact parameters: grazing, interpolated-radius, exact-radius,
  //lowest-radius 
#define nimpact 4
  char *ipdesc[nimpact]={"grazing","interpolated","coincident",
			 "lowest"};
  PREC_RES res[nimpact];

  //Print filling of arrays if necessary to debug
  if(outp)
    fprintf(outp,
	    "#radius        refr          ex_cons       ex_out        ex_in\n");
  PREC_RES rm=spacing*nrad*1.0;
  for(i=0;i<nrad;i++){
    //check for equispaced radius skiping 0
    rad[i]=spacing*(1.0+i);
    //check for nonbending rays
    refr[i]=1.0;
    //extinction: constant
    ex[0][i]=alpha;
    //extinction: incresing outward
    ex[1][i]=alpha*rad[i];
    //extinction: incresing inward
    ex[2][i]=alpha*(rm-rad[i]);

    //Print filling of arrays if necessary to debug
    if(outp)
      fprintf(outp,
	      "%-14.9g%-14.9g%-14.9g%-14.9g%-14.9g\n"
	      ,rad[i],refr[i],ex[0][i],ex[1][i],ex[2][i]);
  }
  PREC_RES b[nimpact]={rad[nrad-2], (rad[(3*nrad)/4-1] + rad[(3*nrad)/4])/2.0,
		       rad[nrad/2-1], rad[0]};

  //Test the atmosphere with constant extinction at different radius
  //$res=2 ex \sqrt{rm^2-b^2}$
  for(i=0;i<nimpact;i++)
    res[i]=2*alpha*sqrt(rm*rm-b[i]*b[i]);

  status+=test_tau_monorad_norefr_now(b,rad,refr,ex,nrad,0
				      ,res,nimpact,ipdesc,acceptrelerror);

  //Test the atmosphere with increasing extinction outwards
  //$res=\alpha ( rm\sqrt{rm^2-b^2} + b^2 \log{\sqrt{(rm/b)^2-1}+rm/b})$
  // Calc function:\par
  //{\em define
  //res(b)=alpha*(rm*sqrt(rm*rm-b*b)+b*b*ln(sqrt((rm/b)\^2-1)+rm/b)/ln(10))}
  //\par $ex=\alpha r$
  for(i=0;i<nimpact;i++)
    res[i]=alpha * (rm * sqrt( rm * rm - b[i] * b[i]) + b[i] * b[i] *
		    log( ( sqrt( rm * rm / b[i] / b[i] - 1) + rm ) / b[i]) );

  status+=test_tau_monorad_norefr_now(b,rad,refr,ex,nrad,1
				      ,res,nimpact,ipdesc,acceptrelerror);

  //Test the atmosphere with increasing extinction inwards
  //$res=\alpha ( rm\sqrt{rm^2-b^2} - b^2 \log{\sqrt{(rm/b)^2-1}+rm/b})$
  // Calc function:\par
  //{\em define
  //res(b)=alpha*(rm*sqrt(rm*rm-b*b)-b*b*ln(sqrt((rm/b)\^2-1)+rm/b)/ln(10))}
  //\par $ex=\alpha (rm-r)$
  for(i=0;i<nimpact;i++)
    res[i]=alpha * (rm * sqrt( rm * rm - b[i] * b[i]) - b[i] * b[i] * 
		    log( ( sqrt( rm * rm / b[i] / b[i] - 1) + rm ) / b[i] ) );

  status+=test_tau_monorad_norefr_now(b,rad,refr,ex,nrad,2
				      ,res,nimpact,ipdesc,acceptrelerror);

  //Output partial result
  if (status)
    test_result("PARTIAL FAIL: %li/%li tests with error bigger than %.5g\n"
		"################################################\n\n"
		,status,ntests-ntb,acceptrelerror);
  else
    test_result("PARTIAL SUCCEED: all %li tests with error less than %.5g\n"
		"################################################\n\n"
		,ntests-ntb,acceptrelerror);

  free(ex[0]);
  return status;
}


/* \fcnfh
   check totaltau funtion

   @returns number of errors
*/
int test_tau()
{
  long status = 0;

  //First check for monospaced radius with 10,100,1000 samples, alpha=1
  status += test_tau_monorad_constrefr(10,   1.0, 1.0);
  status += test_tau_monorad_constrefr(100,  1.0, 0.1);
  status += test_tau_monorad_constrefr(100,  1.0, 1.0);
  status += test_tau_monorad_constrefr(1000, 1.0, 0.1);
  status += test_tau_monorad_constrefr(1000, 1.0, 1.0);
  status += test_tau_monorad_constrefr(10000,1.0, 0.1);
  
  //Now with a different alpha
  status += test_tau_monorad_constrefr(10,   8.0, 1.0);
  status += test_tau_monorad_constrefr(100,  8.0, 1.0);
  status += test_tau_monorad_constrefr(1000, 8.0, 1.0);
  

  /* TD: test for nonmonospaced radius */
  /* TD: test for bent rays */

  return status;
}



int 
main(int argc, char *argv[])
{
  int status = 0;

  if(argc>1)
    outp=fopen(argv[1],"w");

  status += test_tau( );

  fclose(outp);

  if(status){
    fprintf(stdout,
	    "\nslantpath.txc result is FAILURE: %i error%s found out of %li test%s\n"
	    ,status,status>1?"s":"",ntests,ntests>1?"s":"");
    exit(EXIT_FAILURE);
  }
  fprintf(stdout,
	  "\nslantpath.txc result is SUCCESS: no errors found out of %li test%s\n"
	  ,ntests,ntests>1?"s":"");
  exit(EXIT_SUCCESS);
}
