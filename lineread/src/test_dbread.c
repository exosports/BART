/*
 * test_dbread.c - Test for the read drivers. Part of Transit program.
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
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */


#include <transit.h>
#include <lineread.h>
#include <math.h>

//Expect results in nanometers
double tli_fct=1e-7;
const char *tli_fct_name="nanometers";

/*\fcnfh 
  Routine for testing
*/
int main(int argc,
	 char **argv)
{

  int maxnamelen=100;
  char zfile[maxnamelen],file[maxnamelen],dbread_wl_unit[maxnamelen];
  int numisot;
  char **isonames;
  PREC_NREC (*dbread_main)(char *,  struct linedb **,
			   float, float, char *, 
			   PREC_ZREC ***, PREC_ZREC **, 
			   PREC_ZREC **, PREC_CS ***, int *, int *, 
			   char ***);

  test_fill_defaults(file, zfile, dbread_wl_unit, &numisot, 
		     &isonames, dbread_main);

  struct linedb *lines,*lp;
  PREC_NREC n;
  float wl1,wl2,szb,endb,tf1;
  int i,nbins,qb[numisot],ti1,nbinsgf;
  int ans;
  PREC_ZREC **Z,*T,*isomass;
  PREC_CS **isocs;
  int nT,nIso;
  char rc;

  wl1=1748.5;
  wl2=1828.5;
  if(argc>1){
    wl1=atof(argv[1]);
    if(argc>2)
      wl2=atof(argv[2]);
  }
  nbins=20;
  nbinsgf=10;

  printf("Welcome!\n");
  fprintf(stderr,"Lower wavelength (%s)[%g]?: ",dbread_wl_unit,wl1);
  if((tf1=readds(stdin,&rc,NULL,0))!=0||rc==-1)
    wl1=tf1;

  fprintf(stderr,"Upper wavelength (%s)[%g]?: ",dbread_wl_unit,wl2);
  if((tf1=readds(stdin,&rc,NULL,0))!=0||rc==-1)
    wl2=tf1;

  fprintf(stderr,"Number of wavelength bins[%i]?: ",nbins);
  if((tf1=readl(stdin,&rc))!=0||rc==-1)
    nbins=ti1;

  fprintf(stderr,"Number of GF bins[%i]?: ",nbinsgf);
  if((tf1=readl(stdin,&rc))!=0||rc==-1)
    nbinsgf=ti1;

  char namelow[30];
  memset(namelow,0,30);
  fprintf(stderr,"Name of low energy output file(%s): ",namelow);
  fgets(namelow,29,stdin);
  ti1=0;
  while(namelow[ti1]!='\0')
    ti1++;
  if(ti1&&namelow[ti1-1]=='\n')
    namelow[ti1-1]='\0';

  fprintf(stderr,"Reading '%s'\n--------------------------------\n",file);

  gabby_dbread=2;
  n = (*dbread_main)(file, &lines, wl1, wl2, zfile, &Z, &T,
		     &isomass, &isocs, &nT, &nIso, NULL);
  fprintf(stderr,"--------------------------------\n");

  if(n<0)
    transiterror(TERR_SERIOUS,
		 "dbread_%s() return error status %i.\n",DRIVERNAME, n);

  ti1=(int)(log10(n)+1);

  printf("\ndbread_%s() test results:\n"
	  " Chosen wavelength range was from %.2f to %.2f [%s]\n"
	  " %*li lines read\n"
	  " Choosing %i equal-sized bins, the result is\n"
	  ,DRIVERNAME, wl1, wl2, dbread_wl_unit, ti1, n, nbins);

  szb=(wl2-wl1)/nbins;
  lp=lines;
  double maxgf=0,mingf=1;
  long maxgfi,mingfi,maxgfei,mingfei,maxeli,mineli;
  double mingfe=10000,maxgfe=0;
  double minel=1e10,maxel=0;
  double gfe;
  if(!nbins)
    fprintf(stderr,"  hmmm, you chose 0 bins!\n");
  for(i=0;i<nbins;i++){
    memset(qb,0,sizeof(*qb)*4);
    endb=wl1+(i+1)*szb;
    //    fprintf(stderr,"KK %g %f\n",lp->wl,endb);
    while((lp->wl)<endb&&lp-lines<n){
      if(lp->elow>maxel){
	maxel=lp->elow;
	maxeli=lp-lines;
      }
      else if(lp->elow<minel){
	minel=lp->elow;
	mineli=lp-lines;
      }
      if(lp->gf>maxgf){
	maxgf=lp->gf;
	maxgfi=lp-lines;
      }
      else if(lp->gf<mingf){
	mingf=lp->gf;
	mingfi=lp-lines;
      }
      gfe=lp->gf*exp(-lp->elow/2500);
      if(gfe>maxgfe){
	maxgfe=gfe;
	maxgfei=lp-lines;
      }
      else if (gfe<mingfe){
	mingfe=gfe;
	mingfei=lp-lines;
      }
      qb[lp->isoid]++;
      lp++;
    }

    printf(" %*i = %i + %i + %i + %i lines shorter than %.3f\n"
	   ,ti1,qb[0]+qb[1]+qb[2]+qb[3],qb[0],qb[1],qb[2],qb[3],endb);
  }

  fprintf(stderr,
	  "\nmax GF: %.8g (position %li)\n"
	  "min GF: %.8g (position %li)\n"
	  "max GF*exp(-E/2500): %.5g (%li)\n"
	  "min GF*exp(-E/2500): %.5g (%li)\n"
	  "max Elow: %.8g (%li)\n"
	  "min Elow: %.8g (%li)\n"
	  ,maxgf,maxgfi,mingf,mingfi
	  ,maxgfe,maxgfei,mingfe,mingfei
	  ,maxel,maxeli,minel,mineli);

  //logarithmic binning of GF
  long gfb[nbinsgf];
  long bin;
  lp=lines;
  printf("     Logarithmic binning result of GF\n"
	 "Bin beggining at: ");
  for(bin=0;bin<nbinsgf;bin++)
    printf("%10.4g",mingf * ( pow(maxgf/mingf,(float)bin/nbinsgf) ));
  printf("\n");
  for(i=0;i<nbins;i++){
    memset(gfb,0,sizeof(gfb));
    endb=wl1+(i+1)*szb;
    //    fprintf(stderr,"KK %g %f\n",lp->wl,endb);
    while((lp->wl)<endb&&lp-lines<n){
      bin=nbinsgf*log(lp->gf/mingf)/log(maxgf/mingf);
      if(bin==nbinsgf) bin--;
      gfb[bin]++;
      lp++;
    }

    printf("[%6.1f,%6.1f] = "
	   ,endb-szb,endb);
    for(bin=0;bin<nbinsgf;bin++)
      printf("%10li",gfb[bin]);
    printf("\n");
  }


  //Binning low energy
  long elb[nbinsgf];
  lp=lines;
  printf("      Binning result of lower energy\n"
	 "Bin beggining at: ");
  for(bin=0;bin<nbinsgf;bin++)
    printf("%10.4g",minel + bin*(maxel-minel)/(nbinsgf-1));
  printf("\n");
  for(i=0;i<nbins;i++){
    memset(elb,0,sizeof(elb));
    endb=wl1+(i+1)*szb;
    //    fprintf(stderr,"KK %g %f\n",lp->wl,endb);
    while((lp->wl)<endb&&lp-lines<n){
      bin=(lp->elow-minel)*(nbinsgf-1)/(maxel-minel);
      if(bin==nbinsgf) bin--;
      elb[bin]++;
      lp++;
    }

    printf("[%6.1f,%6.1f] = "
	   ,endb-szb,endb);
    for(bin=0;bin<nbinsgf;bin++)
      printf("%10li",elb[bin]);
    printf("\n");
  }

  //Counting 10000 lowest energy values
  if(namelow){
    long nlow[10001],tmpl;
    memset(nlow,0,10001*sizeof(long));
    lp=lines;
    while(lp-lines<n){
      tmpl=lp->elow>10000?10000:lp->elow;
      nlow[tmpl]++;
      lp++;
    }
    FILE *out=fopen(namelow,"w");
    long acum=0;
    for(bin=0;bin<10000;bin++){
      acum+=nlow[bin];
      fprintf(out,"%-18li%-18li\n",bin,acum);
    }
    fclose(out);
  }
  /*
  printf("\nTemperature points:");
  for(ti1=0;ti1<nT;ti1++)
    printf(" %.1f",T[ti1]);
  printf("\n");

  for(i=0;i<nIso;i++){
    printf("\n%s:",isonames[lines[i].isoid]);
    for(ti1=0;ti1<nT;ti1++)
      printf(" %.3f",Z[i][ti1]);
    printf("\n");
  }
  */

  printf("\nWanna know the value of a single record?\n"
	  "If so, write record number (range 0 - %li), else "
	    "press 'enter': "
	  ,n-1);


  while((ans=readl(stdin,&rc))>=0&&rc==-1){
    if(ans<n&&ans>=0){
      lp=lines+ans;
      printf("Record Position: %li\nWavelength: %.10g [%s]\n"
	     "Wavenumber: %.10g [cm]\n"
	     ,lp->recpos,lp->wl,tli_fct_name,1/lp->wl/tli_fct);
      printf("Lower Energy Level: %.10g [cm-1]\ngf: %.10g\n"
	      , lp->elow, lp->gf);
      printf("Isotope's Name: %s\nFreq: %.10g"
	      ,isonames[lp->isoid],LS/lp->wl/tli_fct);
    }
    else
      printf("\nInvalid record number, so ...");

    printf("\nWanna know the value of another single record?\n"
	    "If so, write the record number (range 0 - %li), else just "
	    "press 'enter': "
	    ,n-1);
  }

  printf("\nHave a good day!...\n");

  return EXIT_SUCCESS;

}

