/*
 * test_slantpath.c - Tests to check correct working of slantpath.c.
 *                    Component of the transit program
 *
 * Copyright (C) 2004-2006 Patricio Rojo (pato@astro.cornell.edu)
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
#include <test_common.h>
#include "slantpath.c"

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
tau_now(double ip,
	double *rad,
	double *refr,
	double *ex,
	long nrad,
	double res,
	int type,
	char *desc)
{
  double result,err;
  int status=0;

  result=totaltau(ip,rad,refr,ex,nrad,type);
  err=fabs(1-result/res);
  test_result("        %-25s: observed %13.10g (error %g)\n"
	      ,desc,result,err);
  if(err>maxerr)
    status++;
  test++;

  return status;
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

  status+=tau_now(ip,rad,refr,ex,nrad,res,1,"RefIdx constant technique");
  status+=tau_now(ip,rad,refr,ex,nrad,res,2,"Ray bending technique");

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
  return - alpha * (rm * ip * sqrt( rat * rat - 1 ) - ip * ip *
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

  test_result("\nTesting OPTICAL DEPTH for extinction %s, alpha=%g\n",str,alpha);
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


/* \fcnfh
   This function gives the expect result or the value of tau for an
   atmosphere where tau remains constant
*/
double *mod_ctau(double prm, double *res, double star, double ipmax, double first, long nip, double toomuch)
{
  //return result if that is what user wanted.
  if(res){
    double rath=ipmax/star;
    double ratl=first*ipmax/star;
    *res= - exp(-prm)*(rath*rath-ratl*ratl)
      - exp(-toomuch)*ratl*ratl
      + rath*rath;
    return NULL;
  }

  static double *tau;
  tau=(double *)calloc(nip,sizeof(double));

  while(--nip)
    tau[nip]=prm;
  tau[0]=prm;

  return tau;
}


/* \fcnfh
   This function gives the expect result or the value of tau for an
   atmosphere where tau increase outwards.
*/
double *mod_itau(double prm, double *res, double star, double ipmax, double first, long nip, double toomuch)
{
  //return result if that is what user wanted.
  if(res){
    double rath=ipmax/star;
    double ratl=first*ipmax/star;
    *res =  - 2 * ( exp(-prm*ipmax*first)*(first*ipmax+1/prm) - exp(-prm*ipmax)*(ipmax+1/prm) )
      / star / star / prm;
    *res+= - exp(-toomuch)*ratl*ratl;
    *res+= rath*rath;
    return NULL;
  }

  static double *tau;
  tau=(double *)calloc(nip,sizeof(double));

  double delt = (1-first) / (nip-1);
  while(nip--)
    tau[nip] =  prm * ipmax * (1 - nip * delt);

  return tau;
}


/* \fcnfh
   This function gives the expect result or the value of tau for an
   atmosphere where tau increasing inwards.
*/
double *mod_dtau(double prm, double *res, double star, double ipmax, double first, long nip, double toomuch)
{
  //return result if that is what user wanted.
  if(res){
    double rath=ipmax/star;
    double ratl=first*ipmax/star;
    *res = - 2 * ( (ipmax-1/prm) - exp(-prm*ipmax*(1-first))*(ipmax*first-1/prm) )
      / star / star / prm;
    *res+= - exp(-toomuch)*ratl*ratl;
    *res+= rath*rath;
    return NULL;
  }

  static double *tau;
  tau=(double *)calloc(nip,sizeof(double));

  double delt = (1-first) / (nip-1);
  while(nip--)
    tau[nip] =  prm * ipmax * nip * delt;

  return tau;
}


int
mod_now(double *tau,
	long last,
	double toomuch,
	prop_samp *ip,
	struct geometry *sg,
	double res,
	int type,
	char *desc)
{
  double result,err;
  int status=0;

  result=modulationperwn(tau,last,toomuch,ip,sg,type);
  err=fabs(1-result/res);
  test_result("          %-25s: observed %13.10g (error %g)\n"
	      ,desc,result,err);
  if(err>maxerr)
    status++;
  test++;

  return status;
}


int
mod_nip(struct geometry *sg,
	double ipmax,
	double first,
	long nip,
	double toomuch,
	double res,
	double *(*tauf)(),
	double tprm)
{
  int status=0,i;
  double delt=ipmax*(1-first)/(nip-1);

  test_result("        Using %li radial layers. Spacing is %g\n"
	      ,nip,delt);

  double *tau=tauf(tprm, NULL, sg->starrad*sg->starradfct,
		   ipmax, first, nip, toomuch);
  if(!tau){
    fprintf(stderr, 
	    "Somehow, Tau was not allocated by the (tauf)(). Aborting\n");
    exit(EXIT_FAILURE);
  }

  prop_samp ip;
  ip.v=(PREC_RES *)calloc(nip,sizeof(double));
  ip.n=nip;
  ip.o=1;
  ip.fct=1;
  ip.i=ipmax;
  ip.f=first * ipmax;
  for(i=0;i<nip;i++)
    ip.v[i] = first * ipmax + (nip - 1 - i) * delt;

  status+=mod_now(tau, nip-1, toomuch, &ip, sg, res, 1, "Full eclipse, no LD");

  free(tau);

  return status;
}


int
mod_first(struct geometry *sg,
	  double ipmax,
	  double first,
	  long *nip,
	  int nnip,
	  double toomuch,
	  double *(*tauf)(),
	  double tprm)
{
  int status=0,i;
  double result;

  tauf(tprm, &result, sg->starrad*sg->starradfct, ipmax, first, 0, toomuch);
  test_result("      Atmosphere starting at %g. Expected result is %g\n"
	      ,first*ipmax,result);
  for(i=0;i<nnip;i++)
    status+=mod_nip(sg, ipmax, first, nip[i], toomuch, result, tauf, tprm);

  return status;
}


int
mod_ip(struct geometry *sg,
       double ipmax,
       double *first,
       int nfirst,
       long *nip,
       int nnip,
       double toomuch,
       double *(*tauf)(),
       double tprm)
{
  int status=0,i;

  if(ipmax>sg->starrad)
    return 0;
  test_result("    For a planet of radius %g\n",ipmax);
  for(i=0;i<nfirst;i++)
    status+=mod_first(sg, ipmax, first[i], nip, nnip, toomuch, tauf, tprm);

  return status;
}


int
mod_star(double star,
	 double *ipmax,
	 int nipmax,
	 double *first,
	 int nfirst,
	 long *nip,
	 int nnip,
	 double toomuch,
	 double *(*tauf)(),
	 double tprm)
{
  int status=0,i;

  test_result("  For a star of radius %g\n",star);

  struct geometry sg;
  sg.smaxis=1;
  sg.smaxisfct=1;
  sg.time=0;
  sg.timefct=1;
  sg.incl=90;
  sg.inclfct=PI/180.0;
  sg.ecc=0;
  sg.eccfct=1;
  sg.lnode=0;
  sg.lnodefct=1;
  sg.aper=1;
  sg.aperfct=1;
  sg.starmass=1;
  sg.starmassfct=1;
  sg.starrad=star;
  sg.starradfct=1;
  sg.x=0;
  sg.y=0;

  for(i=0;i<nipmax;i++)
    status+=mod_ip(&sg, ipmax[i], first, nfirst, nip, nnip, toomuch, 
		   tauf, tprm);

  return status;
}


int
mod_tau(double *star,
	int nstar,
	double *ipmax,
	int nipmax,
	double *first,
	int nfirst,
	long *nip,
	int nnip,
	double toomuch,
	double *(*tauf)(),
	char *desc,
	double tprm)
{
  int status=0,i;

  test_result("\nTesting MODULATION for %s; parameter %g\n",desc,tprm);
  for(i=0;i<nstar;i++)
    status+=mod_star(star[i], ipmax, nipmax, first, nfirst, nip, nnip,
		     toomuch, tauf, tprm);

  return status;
}

int
test_mod()
{
  int status=0;
  double ipmax[]={10, 100, 1000};
  long nip[]=    {10, 100, 1000, 10000};
  double first[]={.9, .75, .5, .1};
  double star[]= {10,100,1000};
  double toomuch=30;

  int n_ipmax=sizeof(ipmax)/sizeof(double);
  int n_nip=sizeof(nip)/sizeof(long);
  int n_first=sizeof(first)/sizeof(double);
  int n_star=sizeof(star)/sizeof(double);

  status+=mod_tau(star, n_star, ipmax, n_ipmax, first, n_first, nip, n_nip,
		  toomuch, &mod_ctau, "Constant Tau", .01);
  status+=mod_tau(star, n_star, ipmax, n_ipmax, first, n_first, nip, n_nip,
		  toomuch, &mod_itau, "Increasing Tau", .01);
  status+=mod_tau(star, n_star, ipmax, n_ipmax, first, n_first, nip, n_nip,
		  toomuch, &mod_dtau, "Decreasing Tau", .01);

  return status;
}



int
test_tau()
{
  int status=0;
  double rmax[]={10, 100, 1000};
  long nrad[]={10, 100, 1000};
  double ip[]  ={0.1, 0.5, 0.75, 0.9};
  int n_rmax=sizeof(rmax)/sizeof(double);
  int n_nrad=sizeof(nrad)/sizeof(long);
  int n_ip=sizeof(ip)/sizeof(double);

  status += tau_ex(rmax, n_rmax, ip, n_ip, nrad, n_nrad, 1.0, "increasing outwards");
  status += tau_ex(rmax, n_rmax, ip, n_ip, nrad, n_nrad, -1.0, "increasing inwards");
  status += tau_ex(rmax, n_rmax, ip, n_ip, nrad, n_nrad, 0, "kept constant");
  return status;
}


int 
main(int argc, char *argv[])
{
  int status = 0;

  if(argc>1)
    outp=fopen(argv[1],"w");

  status += test_tau( );
  status += test_mod( );

  if(outp)
    fclose(outp);

  if(status){
    fprintf(stdout,
	    "\n###########################################################################\n"
	    "slantpath.c result is FAILURE: %i difference%s bigger than %g\n"
	    "                                 found out of %li test%s\n"
	    ,status,status>1?"s":"",maxerr,test,test>1?"s":"");
    exit(EXIT_FAILURE);
  }
  fprintf(stdout,
	  "\n#############################################################################\n"
	  "slantpath.c result is SUCCESS: no differences bigger than %g\n"
	  "                                 found out of %li test%s\n"
	  ,maxerr,test,test>1?"s":"");

  exit(EXIT_SUCCESS);
}
