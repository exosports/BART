/*********************************************************************
 * 
 * Copyright (C) 2006 Patricio Rojo
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *    MA  02110-1301, USA.
 * 
 *
 **********************************************************************/

//\delfh
#if 0 				/* Symmetric filling of an asymetric
				   array!! voigt.txc */
    //Advance the created \vr{aint} pointer half way
    //through, later we have to bring it back before reusing it (Line
    //\ref{aint.back}).
    aint+=ob;

    //If there is an odd number of points, then evaluate mid point, and
    //set the backward filling pointer \vr{vprod} to the center of the
    //array.
    if(nint&1){                  /* odd number of points */
      transitDEBUG(20,verblevel,"ODD\n");
      vprod=aint;
      x=0;
      voigtxy(x,y,aint,eps,alphaD);
    }
    //Else, if there is an even number then set the backward filling
    //pointer to one position ahead than the forward filling pointer.
    else{
      vprod=aint--;
    }

    //Symmetric filling of answer. \vr{y} is given by $\sqrt{\ln 2}\cdot
    //|\Delta\wn| /\alpha_D$, where $\Delta\wn$ is with respect to the
    //shifted center of the spectral line according to the index \vr{j}:
    //$\Delta\wn= \wn_{bin}-\Wn =
    //\vrm{ddwn}\cdot\vrm{i}-\vrm{j}\cdot\vrm{dcshft}$, where
    //$\Wn=\wn_{central-bin}+\wn_{center-shift}$
    for(i=1;i<=ob;i++){
      x=SQRTLN2*fabs(dint*i-j*dcshft)/alphaD;
      transitDEBUG(22,verblevel,"i:%i,x:%g,d:%i\n",i,x,vprod-aint);
      voigtxy(x,y,aint+i,eps,alphaD);
      *(vprod-i)=*(aint+i);
      if(i<20) transitDEBUG(22,verblevel,"aint[%i]:%g <- %g\n",i,aint[i],x);
    }

    //Bringing back the value of \vr{aint} for its reuse
    //\linelabel{aint.back}
    aint=vprod-ob;
    transitDEBUG(20,verblevel,
		 "central values= a[%i]:%g\ta[%i]:%g\ta[%i]:%g\n"
		 ,ob-1,aint[ob-1],ob,aint[ob],ob+1,aint[ob+1]);
#endif /* Obsolete */
//\deluh

//\delfh
#if 0 				/* obsolete */
/*\fcnfh
   Computes Voigt Profile at specified wavelength points.
*/
inline int voigtwl(int J,      //Number of bins that the profile spans
			       //in both directions from the center of
		               //the bin.
		   float dl,   //wavelength separation `$\Delta\lambda$'.
		   double lc,  //central dwavelength of the bin
			       //`$\lambda_c$'.
		   PREC_VOIGTP alphaL,  //Lorentz width `$\alpha_L$'.
		   PREC_VOIGTP alphaD,	//Doppler width `$\alpha_D$'.
		   PREC_VOIGT *vpro)    //Array of length $2J+1$, where
					//the final profile will be
					//stored
{
  double y;

  //The function voigtxy(x, y, res, eps, alphaD) receives 
  //$x=\sqrt{\ln 2}(\wn-\wn_0)/\alpha_D$, $y=\sqrt{\ln
  //2}\alpha_L/\alpha_D$, the array `res' where it is going to store the
  //result, a factor for resolution `eps' and the doppler width
  //`alphaD'.\par
  //We set the value of y:
  y=SQRTLN2*alphaL/alphaD;

  //We first want to get a fine resolution Voigt profile, it has to be
  //at least of the bin span set by J
  //\begin{align}
  //\Delta\wn_{tot}=& \left|
  //\frac{1}{\lambda_c+J\Delta\lambda+\Delta\lambda/2} - 
  //\frac{1}{\lambda_c-J\Delta\lambda-\Delta\lambda/2} \right| \notag\\
  //=& \frac{1}{\Delta\lambda}\frac{2+J}
  //{\left(\frac{\lambda_c}{\Delta\lambda}\right)^2 
  //- \left(J+\frac{1}{2}\right)^2}
  //\end{align}
  // and also have a fine resolution
  //such that it is the minimum of `timesfiner' times finer than the
  //minimun separation obtained by the given $\Delta\lambda$, or
  //`timesalpha` times finer than $\alpha_D$. The minimum wavenumber
  //separation of the bin is at value +J where it is given by\par
  //\begin{align}
  //\delta\wn_{min}=&\frac{1}{\lambda_c+J\Delta\lambda-\Delta\lambda/2}
  //- \frac{1}{\lambda_c} -
  //\left(\frac{1}{\lambda_c+J\Delta\lambda+\Delta\lambda/2} -
  //\frac{1}{\lambda_c}\right)\notag\\
  //=&\frac{1}{\Delta\lambda}\frac{1}{
  //\left(\frac{\lambda_c}{\Delta\lambda}+J\right)^2 -
  //\frac{1}{4}}
  //\end{align}


  //If $\lambda_z$ is the bin wavelength that is closest to the center
  //of the voigt profile, and $\lambda_c$ is the real wavelength center
  //then the upper and lower wavenumber of the ith bin from the central
  //one ($\lambda_z$) are
  //\begin{align}
  //\wn_{\binomtn{u\\l}}=&\frac{1}{\lambda_z+i\Delta\lambda \mp
  //\Delta\lambda/2} - \frac{1}{\lambda_c}\\
  //=&\frac{2\left(\lambda_c-\lambda_z-i\Delta\lambda\right) \pm
  //\Delta\lambda}{\left(2\lambda_z+2i\Delta\lambda \mp \Delta\lambda
  //\right) \lambda_c},
  //\end{align}
  //now, if
  //\[
  //\lambda_c=\lambda_z-\frac{\Delta\lambda}{2}+\frac{m}{M}\Delta\lambda
  //\]
  //we get
  //\begin{align}
  //\wn_{\binomtn{u\\l}}=&\frac{1}{\Delta\lambda}
  //\frac{\binomsm{ +1\\-3}\frac{1}{2} + \frac{m}{M}-2i}
  //{\left(\frac{\lambda_z}{\Delta\lambda}-\left(i \mp
  //\frac{1}{2}\right) \right) \left(\frac{\lambda_z}{\Delta\lambda} -
  //\frac{1}{2} + \frac{m}{M}\right)}.
  //\end{align}
  //Is between the lower and upper limits in the previous expression
  //that we have to integrate the Voigt profile to obtain the voigt
  //modulation for the bin $i$, whence the center of the line is
  //displaced by $m$ from the real center. $m$ is given by
  //\[
  //m=\left\lfloor
  //M\frac{\lambda_c-\lambda_z+\Delta\lambda/2}{\Delta\lambda}+0.5
  //\right\rfloor
  //\]
  //{\bf \Large{Warning} function not complete!!}

  return 0;

}


inline int voigtwl2(int J,float dl,double lc,PREC_VOIGTP alphaL,
		   PREC_VOIGTP alphaD,PREC_VOIGT *vpro)
{
  double y,x;
  int i;

  y=SQRTLN2*alphaL/alphaD;

  for(i=-J;i<0;i++){
    x=SQRTLN2*((double)1.0/(lc*(1-lc/(i*dl))))/alphaD;
    voigtxy(x,y,vpro+i,-1,alphaD);
    transitDEBUG(21,verblevel,"UUU: i:%i v:%g -> v(-2):%g\n"
		 ,i,vpro[i],vpro[-2]);
    transitDEBUG(23,verblevel,
		 "XRRX: l:%i D:%g ad:%g al:%g x:%g y:%g %g\n"
	       ,-2,dl,alphaD,alphaL,25.0,y,vpro[-2]);
  }
  i=0;
  voigtxy(0.0,y,vpro+i,-1,alphaD);
  transitDEBUG(21,verblevel,"UUU: i:%i v:%g -> v(-2):%g\n"
	       ,i,vpro[i],vpro[-2]);
  for(i=1;i<=J;i++){
    x=-SQRTLN2*((double)1.0/(lc*(1-lc/(i*dl))))/alphaD;
    voigtxy(x,y,vpro+i,-1,alphaD);
    transitDEBUG(21,verblevel,"UUU: i:%i v:%g -> v(-2):%g\n"
		 ,i,vpro[i],vpro[-2]);
  }

  transitDEBUG(21,verblevel,
	       "XXXX: l:%i D:%g ad:%g al:%g x:%g y:%g %g\n"
	       ,0,(double)dl,alphaD,alphaL,25.0,y,vpro[0]);

  return 1;
}
#endif /* obsolete */
//\deluh


/* \fcnfh
   read a double value

   @returns 0 on success
            first character read if it was not appropiate
*/
char readd2(FILE *fp, 
	    double *f)
{
  char car,prim=0;
  double store=0;
  enum {ent,dec,man,uman} stage=ent;
  _Bool atleast=0,atleastman=0;
  double dlev=10.0;
  int mant=0,mans=1,ents=1;

  stage=ent;
  while(fread(&car,1,1,fp)){
    if(!prim)
      prim=car;
    if(car=='.'){
      if(stage==ent){
	stage=dec;
	continue;
      }
      *f=ents*store;
      return 0;
    }
    if(car=='+'||car=='-'){
      switch(stage){
      case uman:
	stage=man;
	mans=','-car;
	break;
      case ent:
	if(!atleast){
	  ents=','-car;
	  break;
	}
      default:
	if(!atleast)
	  return prim;
	*f=ents*store;
	return 0;
      }
    }
    if(car=='e'||car=='E'){
      if(stage!=man&&stage!=uman&&atleast){
	stage=uman;
	continue;
      }
      if(!atleast)
	return prim;
      *f=ents*store;
      return 0;
    }
    if(car<'0'||car>'9'){
      if(!atleast)
	return prim;
      *f=ents*store;
      return 0;
    }
    atleast=1;
    switch(stage){
    case ent:
      store=store*10+car-'0';
      break;
    case dec:
      store=store+(car-'0')/dlev;
      dlev*=10;
      break;
    case uman:
      stage=man;
    case man:
      mant=mant*9+car-'0';
      store=store*powi(10,mans*mant);
      break;

    }
  }
}



double powi(double x, int n)
{
  double y;

  y=1;

  for(;n>0;--n){
    while((n&1)==0){
      x*=x;
      n>>=1;
    }
    y=y*x;
  }

  return y;
}



/* \fcnfh
 Works as getopt_long except by the extra parameter 'paramfile' which
 can process long versions of options from a file in the format
 \begin{verb} long_option [=] value\end{verb}.
 Reordering of arguments is also done in the same way as getopt_long
 does.

 @return shortoption value as it would have been returned by
 getopt_long.
*/
int 
getopt_long_files_old(int argc,	/* number of arguments */
		char **argv,	/* agruments list */
		char *shortopts, /* short options accepted, format is
				    the same as getopt(). */
		struct option *getopts,	/* long options accepted, format
					   is the same as
					   getopt_long(). Although
					   optional arguments are not
					   supported */
		int *longidxp,	/* returns position in the array
				   getopts[] for the last selected
				   option */
		char *paramfilelist) /* List of files from which to
					process longoptions. Names are
					separated by commas. Value of
					this variable in the first call
					to this function, is the only
					one that matters. */
{
  static FILE *fp=NULL;
  static char **parampointer;
  int ret;

  if(freed){
    fprintf(stderr,
	    "procopt:: getopt_long_files() was called after a call to\n"
	    "getpropt_free(), which should be the very last {procopt}\n"
	    "function called\n");
    exit(EXIT_FAILURE);
  }

  //Allocate name storing array if it is first time this is called.
  if(!paramfiles)
    parampointer=paramfiles=splitnzero_alloc(paramfilelist,',');

  //process file if one is given and it was not finished
  if(*parampointer){
    //open file if it is first run.
    if(!fp)
      fp=fopen(*parampointer,"r");
    //if file is open then process next line
    if(fp){
      //skip empty lines and commented lines. If this is not the first
      //line read, free the array
      while(1){
	if(line)
	  free(line);
	//if we reach the end of line proceed with next default file or
	//go to command line values.
	if((line=fgets_alloc(fp,NULL))==NULL){
	  fclose(fp);
	  fp=NULL;
	  ++parampointer;
	  return getopt_long_files(argc,argv,shortopts,getopts,
				 longidxp,NULL);
	}
	if(*line&&*line!='#')
	  break;
      }
      //If line was read succesfully, find if option exist and set
      //optarg 
      ret=getoptfrom(line,getopts,longidxp);
    }
    //if file was not opened then is because it doesn't exists. So,
    //continue with next one.
    else{
      ++parampointer;
      return getopt_long_files(argc,argv,shortopts,getopts,
			       longidxp,NULL);
    }
  }
  //If file didn't exist, was not given or is finished processing then
  //process command line parameters
  else
    ret= getopt_long(argc,argv,shortopts,getopts,longidxp);

  //If the option was to process a new parameter file then do that. That
  //option is hence never returned out of this function.
  if(givenparamf>=0&&ret==givenparamf){
    ret=parampointer-paramfiles;
    splitnzero_add(&paramfiles,optarg,',');
    parampointer=paramfiles+ret;
    if(*parampointer==NULL&&*optarg!='\0'){
      fprintf(stderr,
	      "procopt internal error: parampointer is null despite\n"
	      "having added a new file.\n");
      exit(EXIT_FAILURE);
    }

    return getopt_long_files(argc,argv,shortopts,getopts,
			     longidxp,NULL);
  }

  return ret;
}

_Bool
fixedcmp(double d1,
	 double d2,
	 int prec)
{
  union ieee754_double du1,du2;
  _Bool d1b=1;
  _Bool d2b=1;

  du1.d=d1*2.0;
  du2.d=d2*2.0;

  if(du1.d<2){
    d1b=0;
    du1.d=2.0/d1;
  }
  if(du2.d<2){
    d2b=0;
    du2.d=2.0/d2;
  }

  int exp1=du1.ieee.exponent & 0xff;
  int exp2=du2.ieee.exponent & 0xff;

  d1*=powi(2,exp1);
  d2*=powi(2,exp2);

  double prec10=powi(10,prec);
  d1*=prec10;
  d2*=prec10;

  return (long)d1==(long)d2;
}
