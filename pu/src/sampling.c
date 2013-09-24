/*
 * sampling.c - Resampling utilities
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

#include <pu/sampling.h>
#ifdef _USE_GSL
#include <gsl/gsl_spline.h>
#endif

/* \fcnfh
   Do resampling of arrays. If getingx then process the x axis,
   otherwise processing y axis. Errors are not continued in this
   function because it is not expected to run with checks.

   @returns 0 on success
            -1 if either of the arrays has zero elements
*/
inline int
resample(const short getingx,	/* Flag indicating which axis we
				   are processing, if it set to 2 then
				   it frees the static arrays */
	 long flags,		/* flags */
	 long nrefx,		/* number of points in the
				   reference x-axis */
	 double *refx,		/* reference x axis */
	 long noutx,		/* number of points in the x
				   axis of the data to resample
				*/
	 double *outx,		/* x axis of data to resample */
	 long ny,		/* number of y axis to process
				 */ 
	 ...)			/* pairs of reference and output data
				   y-axis (in that order) */
{
  long i;
  double *soutx;
  //'indx' will store the position in the data arrays where the
  //reference value was found.
  //'D' and 'h' will be used to store info for the spline function.
  //'x' is pointing to the array's x-axis array; it is assumed that that
  //won't change between calls to resamplex() and resampley().
  //'n' and 'ndat' will keep the number of elements in the reference and
  //data array, respectively.
  //'in' and 'out' are the input and output arrays of arrays for the
  //y-axis.
  //'i', 'f' and 'm' will be used for binary search.
  static float *t=NULL;
  static long *indx=NULL;
  static double *x;
  static long ndat=0,nout=0;
  double *out,*y;
  long f,m,j;
  double val;
  va_list ap;

  //processing the x array
  switch(getingx){
  case 1:
    //if any of the arrays is not given
    if(nrefx<1||noutx<1){
      fprintf(stderr,
	      "%s:: Calling resamplex() with zero reference\n"
	      "or data points.\n"
	      ,__FILE__);
      return -1;
    }

    //If there is only one reference point, then
    if(nrefx==1){
      //warn if there was more than one requested point, but don't stop,
      //just keep setting everything as that value.
      if(noutx>1)
	fprintf(stderr,
		"%s:: In resamplex(), there was only one element\n"
		"data array, but %li elements in the resample array\n"
		,__FILE__,noutx);
      ndat=nrefx;
      nout=1;
      return 0;
    }
    //Check that we are not ask to extrapolate and that we are within
    //the range
    if((outx[0]<refx[0]||outx[noutx-1]>refx[nrefx-1])||
       (outx[noutx-1]<refx[0]||outx[0]>refx[nrefx-1])){
      fprintf(stderr,
	      "%s:: Data x array (%g to %g), spans outside the values\n"
	      "of reference x array (%g to %g). Sorry, but I cannot "
	      "extrapolate.\n"
	      ,__FILE__,outx[0],outx[noutx-1],refx[0],refx[nrefx-1]);
      exit(EXIT_FAILURE);
    }

    //If it is not first use of this function
    if(indx!=NULL){
      free(indx);
      free(t);
    }
    //'ndat' and 'nout' are assigned nrefx and noutx, and those are
    //diminished in one to speed up calculations.
    //'x' is assigned to the data x-axis.
    ndat=nrefx--;
    nout=noutx--;
    x=refx;
    soutx=outx;

    //set the static variables
    indx = (long  *)calloc(nout,sizeof(long));
    t    = (float *)calloc(nout,sizeof(float));

    //'i' is the initial index to look for in the binary search. It is
    //initialized to zero, but then we just use the value previously
    //found, unless is not in strictly ascendent order. \lin{nonascend}
    i=0;
    f=nrefx;
    for(j=0;j<nout;j++){
      val=outx[j];

      if(val<refx[0]||val>refx[nrefx]){
	fprintf(stderr,
		"%s:: resample(): Looking for a value (%g) that is outside\n"
		"reference axis (%g - %g)\n"
		,__FILE__,val,refx[0],refx[nrefx]);
	exit(EXIT_FAILURE);
      }

      //In case the next output element is not bigger than the previous,
      //then start looking from the beginning
      //\linelabel{nonascend}
      if(val<refx[i])
	i=0;
      else if(val==refx[f]){
	i=indx[j]=f++;
	t[j]=0;
	continue;
      }

      //Do binary search only if the value is bigger than the the bin
      //one position bigger than the previous find or the last element
      //was the first one. The latter happens when either this is the
      //first time the loop is tried or every time that the previous
      //result is the first element, this has a small efficiency loss,
      //but is the best alternative I could think.
      if(!i||val>refx[f]){
	f=nrefx;
	while(f-i>1){
	  m=(f+i)>>1;
	  if(refx[m]>val)
	    f=m;
	  else
	    i=m;
	}
      }
      indx[j]=i;
      t[j]=(val-refx[i])/(refx[i+1]-refx[i]);
    }
    break;
  case 0:
    //check that resamplex was called before.
    if(nout<1||ndat<1){
      fprintf(stderr,
	      "%s:: It seems that 'resamplex()' was not run before this\n"
	      "run of 'resampley()': either nout(%li) or ndat(%li)\n"
	      "are not positive.\n"
	      ,__FILE__,nout,ndat);
      return -2;
    }

    va_start(ap,ny);
    while(ny--){
      y   = va_arg(ap, double *);
      out = va_arg(ap, double *);

      if(ndat==1){
	for(j=0 ; j<ndat ; j++)
	  out[j] = y[0];
	continue;
      }

      switch(flags&SAMP_BITS){
      //use spline interpolation
      case SAMP_SPLINE:
	natcubspline(ndat,x,y,nout,indx,t,out,soutx);
	break;
      case SAMP_LINEAR:
	lineinterpol(ndat,x,y,nout,indx,t,out,NULL);
	break;
      default:
	fprintf(stderr,
		"%s:: Sampling method requested(%li) is not yet\n"
		"implemented\n"
		,__FILE__,flags&SAMP_BITS);
	exit(EXIT_FAILURE);
	break;
      }
    }
    va_end(ap);
    break;
  case 2:
    if (nout&&ndat){
      free(indx);
      free(t);
    }
    indx = NULL;
    t = NULL;
    nout = ndat = 0;
  }

  return 0;
}

/* \fcnfh
   Finds value from dataset x,y to output yout

   @returns 1 on success
*/
int 
lineinterpol(int ndat,		/* Length of x,y array */
	     double *x,		/* x-axis */
	     double *y,		/* y-axis */
	     int n,		/* length of output arrays (indx,t,yout)
				   */
	     long *indx,	/* Integer x-axis reference index */
	     float *t,		/* value from 0 to 1 to indicate how far
				   from the value at 'indx' to the next
				   value */
	     double *yout,	/* Answer array were the interpolated
				   values of x[indx]+t*dx[i] are going to
				   be stored */
	     double *dbgout)	/* Not used in linear interpolation */
{
  int i,idx;

  for(i=0;i<n;i++){
    idx=indx[i];
    if(idx>=ndat||idx<0){
      fprintf(stderr,
	      "%s:: index #%i(%i) is out of dataset\n"
	      "range (0 - %i) in lineinterpol()\n"
	      ,__FILE__,i,idx,ndat-1);
      exit(EXIT_FAILURE);
    }
    if(idx==ndat-1){
      if(t[i]!=0){
	fprintf(stderr,
		"%s:: Trying to extrapolate in last index(%i) while\n"
		"in lineinterpol(). t=%f\n"
		,__FILE__,idx,t[i]);
	exit(EXIT_FAILURE);
      }
      yout[i]=y[idx];
    }
    else
      yout[i]=y[idx]+t[i]*(y[idx+1]-y[idx]);
  }

  return 1;
}


/* \fcnfh
   Evaluates the new values of a dataset according to the indices given
   in 'indx' and fraction 'tt' in that interval.
   It doesn't check for values within the range, that have to be done
   before.
*/
void 
natcubspline(int ndat,		/* Length of (x,y) array */
	     double *x,		/* x-axis */
	     double *y,		/* y-axis */
	     int n,		/* length of output arrays
				   (indx,t,yout) */ 
	     long *indx,	/* Integer x-axis reference index */
	     float *t,		/* value from 0 to 1 to indicate how far
				   from the value at 'indx' to the next
				   value */
	     double *yout,	/* Answer array were the interpolated
				   values of x[indx]+t*dx[i] are going to
				   be stored */
	     double *xref)	/* Stores the reference X value */
{
  long i;

  if(ndat<0){
    fprintf(stderr,
	    "%s:: ndat was not positive(%i) in call to\n"
	    "evalnatcubspline()\n"
	    ,__FILE__,ndat);
    exit(EXIT_FAILURE);
  }

#ifdef _USE_GSL
  gsl_interp_accel acc={0,0,0};
  gsl_interp *spl=gsl_interp_alloc(gsl_interp_cspline,ndat);
  gsl_interp_init(spl,x,y,ndat);
  for(i=0;i<n;i++)
    yout[i]=gsl_interp_eval(spl,x,y,xref[i],&acc);
  gsl_interp_free(spl);
#else
#error non equispaced integration is not working yet without GSL
  float tt;
  double hh;
  long f;
  double h[ndat],D[ndat];
  natcubsplinecoef(ndat, x, y, h, D);

  for(f=0;f<n;f++){
    i=indx[f];
    tt=t[f];
    hh=h[f];
    if(i==ndat-1){
      if(t[i]!=0){
	fprintf(stderr,
		"%s:: Trying to extrapolate in last index(%i) while\n"
		"in natcubspline(). t=%f\n"
		,__FILE__,ndat,t[i]);
	exit(EXIT_FAILURE);
      }
      yout[f]=y[i];
    }
    else{
      yout[f]=y[i]*(1-tt) + y[i+1]*tt
	+ hh*hh*tt*(D[i]*(-2+tt*(3-tt)) + D[i+1]*(tt*tt-1))/6.0;
    }
  }

  if(dout!=NULL)
    memcpy(dout,D,sizeof(double)*ndat);
#endif /* _USE_GSL */

}

/*\fcnfh
  It calculates the 2nd derivatives of the function to interpolate
  storing them in 'D'. Requires the 'y' values and the x-axis spacing in
  'h' of the 'n' points dataset.
*/
inline void
natcubsplinecoef(long n,	/* Number of data points */
		 double *x,	/* y-axis values */
		 double *y,	/* y-axis values */
		 double *h,	/* x-axis spacing */
		 double *D)	/* 2nd-derivatives values */
{
  long i,n1;
  double u,v,w;

  //'n'-1 will be used a lot, so giving it to a new variable
  n1=n-1;
  //Initialization of array D as first derivative
  u=y[0];
  for(i=0;i<n1;i++){
    v=y[i+1];
    h[i]=x[i+1]-x[i];
    D[i+1]=(v-u)/h[i];
    u=v;
  }

  //Natural spline condition. D[0]=0.
  //U=h[i-1]/P[i-1], V=h[i-1], W=h[i], P[i] stored in h[i], where P[i]
  //denotes diagonal coefficient in the Gaussian elimination. h[i] is
  //only restored in \lin{hrest}
  D[0]=u=0;
  w=h[0];
  for(i=1;i<n1;i++){
    v=w;
    w=h[i];
    h[i]=(v+w)*2-u*v;
    D[i]=D[i+1]-D[i]-u*D[i-1];
    u=w/h[i];
  }

  D[n1]=0;
  //Back substitution and restoration of 'h'.
  for(i=n1-1;i>0;i--){
    w=x[i+1]-x[i];
    D[i]=(6*D[i]-w*D[i+1])/h[i];
    //\linelabel{hrest}
    h[i]=w;
  }
}


/* \fcnfh
   Interpolates in arrays (x,y), the value corresponding to 'refx'. 'n'
   is the arrays' length. Uses requested interpolation.

   @returns the interpolated y-value
 */
inline double
interp(double refx,		/* Reference x point */
       double *x,		/* x-axis */
       double *y,		/* y-axis */
       long n,			/* lenght of array */
       int intkind)		/* Kind of interpolation */
{

  switch(intkind){
  case INTERP_LINEAR:
    return lineinterp(refx,x,y,n);
  default:
    break;
  }
  fprintf(stderr,
	  "interp(): Requested interpolation %i, is non existent!\n"
	  ,intkind);
  exit(EXIT_FAILURE);
}



/* \fcnfh
   Linearly interpolates in arrays (x,y), the value corresponding to
   'refx'. 'n' is the arrays' length

   @returns the interpolated y-value
 */
inline double
lineinterp(double refx,		/* Reference x point */
	   double *x,		/* x-axis */
	   double *y,		/* y-axis */
	   long n)		/* lenght of array */
{
  double yout;
  _Bool ascend=(x[1]-x[0]>0);
  n--;

  //Initial range is correct?
  if((ascend && (refx<*x || refx>x[n])) || (!ascend && (refx>*x || refx<x[n]))){
    fprintf(stderr,
	    "ascend: %i, %i, %i; noascend %i, %i\n"
	    "refx:%.10g  x: %.10g - %.10g\n"
	    ,ascend,refx<*x,refx>x[n],refx>*x,refx<x[n]
	    ,refx,x[0],x[n]);
    fprintf(stderr,
	    "interpline():: refx(%.8g) is outside x range (%.8g - %.8g)\n"
	    ,refx,*x,x[n]);
    exit(EXIT_FAILURE);
  }

  //check that enough data points are given
  if(n<2){
    fprintf(stderr,
	    "interpline():: Array of length less than 2.\n");
    exit(EXIT_FAILURE);
  }

  //go through the arrays
  while(1){
    //Is final range correct?
    if(!n){
      if(*x==refx)
	return *y;
      fprintf(stderr,
	      "interp():: Last x-value (%g) is smaller than refx(%g)\n"
	      ,*x,refx);
      exit(EXIT_FAILURE);
    }
    //found!
    if((ascend && x[1]>refx) || x[1]<refx)
      break;
    n--;
    x++;
    y++;
  }

  yout=*y+(refx-*x)*(y[1]-*y)/(x[1]-*x);
  return yout;

}


//\delfh
#ifdef DBGSPLINE
#define na 5
int main(int argc, char **argv)
{
  //  double x[na]={0,1,2,3,4,5,6,7,8,9};
  //  double y[na]={0,1,2,3,4,5,6,7,8,9};
  //  double y[na]={2,1,2,3,2,1,2,3,5,9};
  double x[na]={-3,-1,0,3,4};
  double y[na]={7,11,26,56,29};
  int nn=101;
  double s,tt,hh;
  int i,ii;

  if(argc>1)
    nn=atoi(argv[1]);

  int subn=(nn-1)/(na-1)+0.5;
  nn=subn*(na-1) +1;
  long idx[nn];
  float t[nn];
  double out[nn];
  double D[na];

  for(i=0;i<nn;i++){
    idx[i]=i/subn;
    t[i]=(i%subn) / (float)subn;
  }

  natcubspline(na,x,y,nn,idx,t,out,D);
  /*
  for(i=0;i<na;i++){
    printf("%g\t%g\t%g\t%g\n",x[i],y[i],h[i],D[i]/2);
  }
  */

  printf("#x      i        t              S             S'           S\"/2          S\"'/3\n");
  for(ii=0;ii<nn-1;ii++){
    i=idx[ii];
    tt=t[ii];
    hh=x[i+1]-x[i];
    s=(y[i+1]-y[i])/hh + hh*hh/6.0*(D[i]*(-2+6*tt/hh-3*tt*tt/hh)
				    +D[i+1]*(3*tt*tt/hh-1));

    printf("%-8.3g%-4i%6.2f%15.5g%15.5g%15.5g%15.5g\n"
	   ,x[i]+tt*hh,i,tt,out[ii],
	   s,(D[i]*(1-tt)+D[i+1]*tt)/2.0
	   ,(D[i+1]-D[i])/3.0);
  }
  i=idx[ii];
  tt=t[ii];
  hh=x[i+1]-x[i];

  printf("%-8.3g%-4i%6.2f%15.5g\n"
	 ,x[i],i,tt,out[ii]);

  
}
#endif /* DBGSPLINE */

#ifdef DBGSAMPLING
int main(int argc, char **argv)
{
  //  double x[na]={0,1,2,3,4,5,6,7,8,9};

  long fl=TRU_SAMPLIN;


#define ndat 8
  double x[ndat]= {1,3,5,5.5,7,9.3,11,14};
  double y[ndat]= {2,1,0,3.2,3,4.9,-1,58};
  double y2[ndat]={1,2,1,4.2,5,4.3,-8,28};
  int i;

#define n 12
  double outx[n]={2.3,2.8,5.6,5.7,7,9.29,9.4,10,13,13.5,11.2,14};
  double outy[n],outy2[n];


  resamplex(fl,ndat,x,n,outx);
  resampley(fl,2,y,outy,y2,outy2);

  for(i=0;i<n;i++){
    printf("%12f%12f%12f\n",outx[i],outy[i],outy2[i]);
  }

}
#endif


//\deluh
