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
#define nrad 10
#define nwn  3
#define acceptrelerror .0001

#define test_tau_monorad_norefr_now(nimp,wn,exttype)   do{      \
   for(i=0;i<nimp;i++){                                         \
      result=totaltau(b[i],rad,refr,ex,nrad,wn,&acc);           \
      test_result("While "exttype" extinction, at %s,\n"        \
		  "  observed %.10g, expect %.10g, error %g\n"  \
		  ,ipdesc[i],result,res[i],                     \
		  fabs(result-res[i])/res[i]);                  \
      if(fabs(result-res[i])/res[i]>acceptrelerror)             \
	test_fail(status,"Error bigger than %.5g\n"             \
                  ,acceptrelerror);                             \
      else                                                      \
        test_succeed();                                         \
    }                                                           \
                                               }while(0)

/* \fcnfh
   check monospaced radius

   @returns number of errors
*/
int
test_tau_monorad_constrefr(PREC_RES **ex)
{
  //check for nonbending rays
  PREC_RES refr[nrad]={1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  PREC_RES rad[nrad]={1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  gsl_interp_accel acc={0,0,0};
  int i,status=0;
  PREC_RES result;

  //tests at impact parameters: grazing, interpolated-radius, exact-radius,
  //lowest-radius 
#define nimpact 4
  PREC_RES b[nimpact]={9.0, 7.5, 5.0, 1.0};
  char *ipdesc[nimpact]={"grazing IP","interpolated IP","coincident IP",
			 "lowest IP"};

  //Test the atmosphere with constant extinction at different radius
  //$res=2 ex \sqrt{rm^2-b^2}$
  {
    PREC_RES res[nimpact]={8.71779788708134710448,
			   13.2287565553229529525,
			   17.32050807568877293528,
			   19.89974874213239909468};

    test_tau_monorad_norefr_now(nimpact,0,"constant");
  }

  //Test the atmosphere with increasing extinction outwards
  //$res=\alpha ( rm\sqrt{rm^2-b^2} + b^2 \log{\sqrt{(rm/b)^2-1}+rm/b})$
  // Calc function:\par
  //{\em define
  //res(b)=alpha*(rm*sqrt(rm*rm-b*b)+b*b*ln(sqrt((rm/b)\^2-1)+rm/b)/ln(10))}
  //\par $ex=\alpha r$
  {
    PREC_RES res[nimpact]={60.02215842946226671978,
			   85.57381701507597253714,
			   100.90122906677784957218,
			   100.79868387584142893220};

    test_tau_monorad_norefr_now(nimpact,1,"increasing outwards");
  }

  //Test the atmosphere with increasing extinction inwards
  //Test the atmosphere with increasing extinction outwards
  //$res=\alpha ( rm\sqrt{rm^2-b^2} - b^2 \log{\sqrt{(rm/b)^2-1}+rm/b})$
  // Calc function:\par
  //{\em define
  //res(b)=alpha*(rm*sqrt(rm*rm-b*b)-b*b*ln(sqrt((rm/b)\^2-1)+rm/b)/ln(10))}
  //\par $ex=\alpha (rm-r)$
  {
    PREC_RES res[nimpact]={27.15582044135120432502,
			   46.71374853815355698786,
			   72.30385169010987978062,
			   98.19880354548256201460};

    test_tau_monorad_norefr_now(nimpact,2,"increasing inwards");
  }

  return status;
}


/* \fcnfh
   check totaltau funtion

   @returns number of errors
*/
int test_tau()
{
  int status = 0;
  //Common set of extinctions, first wavelength is constant, second
  //wavelength is increasing k=1*r, third wavelength is decreasing with
  //radius k=1*(rm-r)
  PREC_RES ex[nrad][nwn]={{1.0,  1.0,  9.0},
			  {1.0,  2.0,  8.0},
			  {1.0,  3.0,  7.0},
			  {1.0,  4.0,  6.0},
			  {1.0,  5.0,  5.0},
			  {1.0,  6.0,  4.0},
			  {1.0,  7.0,  3.0},
			  {1.0,  8.0,  2.0},
			  {1.0,  9.0,  1.0},
			  {1.0, 10.0,  0.0}};

  //First check for monospaced radius
  status += test_tau_monorad_constrefr((PREC_RES **)ex);
  

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
