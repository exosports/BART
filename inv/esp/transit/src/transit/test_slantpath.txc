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

double maxerr=1e-4;
long test=0;
FILE *outp=NULL;


double calcex(double alpha, double rm, double ri)
{
  if(alpha==0)
    return 1.0;
  if(alpha<0)
    return - alpha * (rm - ri);
  if(alpha<0)
    return - alpha * (rm - ri);
  
  return alpha * ri;
}


int
tau_dens(double rm,
	 double ip,
	 long nrad,
	 double alpha,
	 double res)
{
  int status=0,i;
  double rad[nrad],ex[nrad],refr[nrad];
  double dr=rm/nrad;
  double result,err;

  test_result("      Using %li layers, an interspacing of %g\n",nrad,dr);
  if(outp)
    fprintf(outp,
	    "\n#rad        extinction  refrac\n");
  for(i=0;i<nrad;i++){
    rad[i]=(i+1)*dr;
    ex[i]=calcex(alpha, rm, rad[i]);
    refr[i]=1;
    if(outp)
      fprintf(outp,
	      "%12.9g%12.9g%12.9g\n"
	      ,rad[i],ex[i],refr[i]);
  }
  result=totaltau(ip,rad,refr,ex,nrad,1);
  err=fabs(1-result/res);
  test_result("      RefIdx constant technique: observed %13.10g (error %g)\n"
		,result,err);
  if(err<maxerr)
    status++;
  result=totaltau(ip,rad,refr,ex,nrad,2);
  err=fabs(1-result/res);
  test_result("      Ray bending technique    : observed %13.10g (error %g)\n"
		,result,err);
  if(err<maxerr)
    status++;

  test+=2;
  return status;
}


int
tau_ip(double rm,
       double ip,
       long *nrad,
       int nr,
       double alpha,
       double (*fcn)())
{
  int status=0,i;
  double res=(*fcn)(alpha,rm,ip);

  test_result("    For rays crossing at %g, expected value is %g\n"
	      ,ip,res);
  for(i=0;i<nr;i++)
    status+=tau_dens(rm,ip,nrad[i],alpha,res);

  return status;
}


int
tau_rad(double rm,
	double *ip,
	int ni,
	long *nrad,
	int nr,
	double alpha,
	double (*fcn)())
{
  int status=0,i;

  test_result("  For planets of radius %g\n",rm);
  for(i=0;i<ni;i++)
    status += tau_ip(rm,ip[i]*rm,nrad,nr,alpha,fcn);

  return status;
}



double const_ex_anal(double alpha, double rm, double ip)
{
  return 2 * sqrt( rm*rm - ip*ip );
}



double incout_ex_anal(double alpha, double rm, double ip)
{
  double rat=rm / ip;
  return alpha * (rm * ip * sqrt( rat * rat - 1 ) + ip * ip *
		  log( sqrt( rat* rat - 1) + rat ) );
}



double incin_ex_anal(double alpha, double rm, double ip)
{
  double rat=rm / ip;
  return alpha * (rm * ip * sqrt( rat * rat - 1 ) - ip * ip *
		  log( sqrt( rat* rat - 1) + rat ) );
}

int
tau_ex(double *rmax,
       int nx,
       double *ip,
       int ni,
       long *nrad,
       int nr,
       double alpha,
       char *str)
{
  int status=0,i;
  double (*analres)();

  test_result("\nFor extinction %s, alpha=%g\n",str,alpha);
  if(alpha==0.0)
    analres=&const_ex_anal;
  else if(alpha<0)
    analres=&incin_ex_anal;
  else
    analres=&incout_ex_anal;
  for(i=0;i<nr;i++)
    status += tau_rad(rmax[i],ip,ni,nrad,nr,alpha,analres);

  return status;
}


int
test_tau()
{
  int status=0;
  double rmax[]={10, 100, 1000};
  long nrad[]={10, 100, 1000};
  double ip[]  ={0.1, 0.5, 0.75, 0.9};

  status += tau_ex(rmax, 3, ip, 4, nrad, 3, 1.0, "increasing outwards");
  status += tau_ex(rmax, 3, ip, 4, nrad, 3, -1.0, "increasing inwards");
  status += tau_ex(rmax, 3, ip, 4, nrad, 3, 0, "kept constant");
  return status;
}


int 
main(int argc, char *argv[])
{
  int status = 0;

  if(argc>1)
    outp=fopen(argv[1],"w");

  status += test_tau( );

  if(outp)
    fclose(outp);

  if(status){
    fprintf(stdout,
	    "\nslantpath.txc result is FAILURE: %i error%s found out of %li test%s\n"
	    ,status,status>1?"s":"",test,test>1?"s":"");
    exit(EXIT_FAILURE);
  }
  fprintf(stdout,
	  "\nslantpath.txc result is SUCCESS: no errors found out of %li test%s\n"
	  ,test,test>1?"s":"");
  exit(EXIT_SUCCESS);
}
