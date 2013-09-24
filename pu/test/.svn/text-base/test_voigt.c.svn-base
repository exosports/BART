/*
 * test_voigt.c - Functions to test voigt.c file
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

#include <getopt.h>
#include <string.h>
#include <pu/iomisc.h>
#include <pu/profile.h>


/*\fcnfh
  Help output
*/
static void synhelp(double eps, double lor, double dop,int m)
{
  fprintf(stderr,"Syntax:\n\ttestvoigt <number_of_points> <line_half_range>"
	  " [...options...]\n\n"
	  "\ttestvoigt <wavenumber_center> -f <wavenumber_file> "
	  "[...options...]\n\n"
	  " Where available options and their defaults are:\n"
	  "  -p <precision>     : Convergence criterion for VOIGTXY. (%.2g)\n"
	  "  -m <#_center_shift>: Number of shifts of the center "
	  "wavenumber (%i)\n"
	  "  -d <Doppler_width> : Exactly that! (%.2g)\n"
	  "  -l <Lorentz_width> : ditto! (%.2g)\n\n"
	  ,eps,m,dop,lor);
  exit(EXIT_SUCCESS);
}

/*\fcnfh
  This testing binary is only compiled if DBGVOIGT is \#defined!
 */
int main(int argc, char **argv)
{
  double eps,dop,lor,wn0,delt;
  int rn,fname;
  FILE *fp=NULL;
  long i,nwn,j;
  PREC_VOIGT **res2,*res, *wn;
  int m;
  char rc;

  eps=-1e-6;
  dop=1;
  lor=1.5;
  fname=0;
  nwn=16;
  m=6;

  if(argc==1)
    synhelp(eps,lor,dop,m);

  while(1){
    rn=getopt(argc,argv,"p:d:l:f:hm:");
    if (rn==-1)
      break;

    switch(rn){
    case 'm':
      m=atoi(optarg);
      if(m<1){
	fprintf(stderr,
		"%s:: Value of 'm' has to be positive (%i)\n"
		,__FILE__,m);
	exit(EXIT_SUCCESS);
      }
      break;
    case 'h':
      synhelp(eps,lor,dop,m);
      break;
    case 'f':
      fname=1;
      if (strcmp("-",optarg)==0){
	fp=stdin;
      }
      else{
	if((fp=fopen(optarg,"r"))==NULL){
	  fprintf(stderr,
		  "%s:: Cannot open file '%s'\n"
		  ,__FILE__,optarg);
	  exit(EXIT_FAILURE);
	}
      }
      break;
    case 'p':
      eps=atof(optarg);
      break;
    case 'd':
      dop=atof(optarg);
      break;
    case 'l':
      lor=atof(optarg);
      break;
    default:
      synhelp(eps,lor,dop,m);
    }
  }

  if(fname){
    if(optind!=argc-1)
      synhelp(eps,lor,dop,m);
    wn0=atof(argv[optind]);
    wn=(PREC_VOIGT *)calloc(nwn,sizeof(PREC_VOIGT));
    i=0;
    while(1){
      if(i>=nwn)
	wn=(PREC_VOIGT *)realloc(wn,(nwn<<=1)*sizeof(PREC_VOIGT));
      if((wn[i]=readds(fp,&rc,NULL,0))<1){
	if(rc!=-1)
	  break;
	else{
	  fprintf(stderr,
		  "%s:: Non-numeric line ignored"
		  " while reading wavenumber file"
		  ,__FILE__);
	  continue;
	}
      }
      i++;
    }
    wn=(PREC_VOIGT *)realloc(wn,(nwn=i)*sizeof(PREC_VOIGT));
    res=(PREC_VOIGT *)calloc(nwn,sizeof(PREC_VOIGT));
    voigtf(nwn,wn,wn0,lor,dop,res,eps);
    for(i=0;i<nwn;i++)
      printf("%.6g\t%.6g\n",wn[i],res[i]);

  }
  else{
    if(argc-optind!=2)
      synhelp(eps,lor,dop,m);
    nwn=atoi(argv[optind++]);
    //\vr{wn0} in the following really stands for full width!
    wn0=atof(argv[optind]);
    res2=(PREC_VOIGT **)calloc(nwn,sizeof(PREC_VOIGT *));
    res2[0]=(PREC_VOIGT *)calloc(nwn*m,sizeof(PREC_VOIGT));
    for(i=1;i<m;i++)
      res2[i]=res2[0]+i*nwn;
    printf("#Total range:%g, Lorentz width:%g, Doppler Width %g\n",wn0,lor,dop);

    voigtn(m,nwn, wn0, lor, dop, res2, eps,0);
    printf("#wl       v(w%5.2g)"
	    ,-wn0/2/(nwn-1));
    for(j=1;j<m;j++)
      printf("  v(w%s%5.2g)"
	      ,j/(float)m<0.5?"":"+",wn0/(nwn-1)*(j/(float)m-0.5));
    printf("\n");

    delt=2*wn0/(nwn-1);
    //Note that after this output, the result is lost because the first
    //bin is a sum of all the others.
    for(i=0;i<nwn;i++){
      printf("%10.6g",-wn0+i*delt);
      for(j=0;j<m;j++){
	printf("%11.4g",res2[j][i]);
	if(i)
	  res2[j][0]+=res2[j][i];
      }
      printf("\n");
    }

    fprintf(stderr,"totals:");
    for(j=0;j<m;j++)
      fprintf(stderr,"%11.4g",res2[j][0]*delt);
    fprintf(stderr,"\n");

  }

  return EXIT_SUCCESS;
}

