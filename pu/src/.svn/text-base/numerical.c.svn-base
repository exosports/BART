/*
 * numerical.c - Various numerical utilities
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

#include <pu/numerical.h>

/* \fcnfh
   Search for index such that val is between arr[index] inclusive and
   arr[index+1] exclusive.

   @returns index that is looked for
            -1 if value is before first element 'i' of array
	    -2 if value is after last element 'f' of array
	    -3 if only one element in array and it was not the looked
               one 
	    -4 if 'f' is equal or smaller than 'i', or 'i' is less than
	       0 
	    -5 if value is the last element of the array
*/
inline int
binsearchie(double *arr,	/* Array of length of at least 'f' */
	    long i,		/* initial search index, cannot be
				   negative */
	    long f,		/* final search index (or array
				   length minus 1) */
	    double val)		/* number to look for in the array */
{
  long m;

  if(arr[i]>val)
    return -1;
  if(arr[f]<val)
    return -2;
  if(arr[f]==val)
    return -5;
  if(i==f&&arr[i]!=val)
    return -3;
  if(f<i&&i<0)
    return -4;

  while(f-i>1){
    m=(f+i)>>1;
    if(arr[m]>val)
      f=m;
    else
      i=m;
  }

  return i;
}


/* \fcnfh
   Search for index such that val is between arr[index] exclusive and
   arr[index+1] inclusive.

   @returns index that is looked for
            -1 if value is before first element 'i' of array
	    -2 if value is after last element 'f' of array
	    -3 if only one element in array and it was not the looked
               one 
	    -4 if 'f' is equal or smaller than 'i', or 'i' is less than
	       0 
	    -6 if it is first index of the array
*/
inline int
binsearchei(double *arr,	/* Array of length of at least 'f' */
	    long i,		/* initial search index, cannot be
				   negative */
	    long f,		/* final search index (or array
				   length minus 1) */
	    double val)		/* number to look for in the array */
{
  long m;

  if(arr[i]>val)
    return -1;
  if(arr[i]==val)
    return -6;
  if(arr[f]<val)
    return -2;
  if(i==f&&arr[i]!=val)
    return -3;
  if(f<i&&i<0)
    return -4;

  while(f-i>1){
    m=(f+i)>>1;
    if(arr[m]<val)
      i=m;
    else
      f=m;
  }

  return i;
}

/*\fcnfh
  binsearch() defaults to an inclusive, exclusive search

  @return binsearchie findings.
*/
inline int 
binsearch(double *arr,
	  long i,
	  long f,
	  double val)
{
  return binsearchie(arr,i,f,val);
}


/* \fcnfh
   Integrate using simpson method and trapezoid at one interval if even
   number, it requires an equispaced x-grid

   @returns value of the integration
*/
inline double
integ_trasim(double dx,
	     double *y,
	     long n)
{
  double restrap=0,res=0;
  long i;

  if(n<2){
    fprintf(stderr,
	    "%s:: integ_trasim: At least 2 points are required to perform\n"
	    "integration\n"
	    ,__FILE__);
    exit(EXIT_FAILURE);
  }

  //use trapezoidal in the last even element
  if(!(n&1)){
    n--;
    restrap=.5*dx*(y[n]+y[n-1]);
  }

  //if there is enough elements do a Simpson integral
  if(n>2){
    //add the middle values, start with the odd elements which will be
    //multiplied by 4
    n--;
    for(i=1;i<n;i+=2)
      res+=y[i];
    res*=2;

    //now the even elements to be multiplied by 2
    for(i=2;i<n;i+=2)
      res+=y[i];
    res*=2;

    //now the borders
    res+=y[0]+y[n];
  }

  //finish'em
  return res*dx/3.0+restrap;
}


/* \fcnfh
   Interpolates a parabola in three points and return requested
   value. X-array must need equispaced. This function doesn't check for
   less than required or non-equispaced elements.

   @returns value interpolated
*/
inline double
interp_parab(double *x,		/* x-array with at least 3 equispaced
				   elements */
	     double *y,		/* y-array with at least 3 elements */
	     double xr)		/* requested x-value to interpolate */
{
  const double dx = x[1] - x[0];
  const double x0 = x[0] / dx;
  const double my = y[0] + y[2] - 2*y[1];
  const double a  = my / 2.0 / dx / dx;
  const double b  = (y[2] - y[1] - (x0 + 1.5) * my) / dx;
  const double c  = y[0] + x0 * ( y[2] - 4*y[1] + 3*y[0] + x0 * my ) / 2.0;

  return xr * xr * a  +  xr * b  +  c;
}


/* \fcnfh
   Interpolates a line in two points and return requested
   value. This function doesn't check for less than required elements.

   @returns value interpolated
*/
inline double
interp_line(double *x,		/* x-array with at least 2 elements */
	    double *y,		/* y-array with at least 2 elements */
	    double xr)		/* requested x-value to interpolate */
{
  const double dx = x[1] - x[0];
  const double m  = (y[1] - y[0]) / dx;

  return y[0] + (xr - x[0]) * m;
}


/* \fcnfh
   return $x^n$, is a faster version of pow() that only works if n is
   integer

   @returns result
*/
double 
powi(double x,
     int n)
{
  double y;
  _Bool negn=n<0;

  y=1;
  if(negn)
    n*=-1;

  for(;n>0;--n){
    while((n&1)==0){
      x*=x;
      n>>=1;
    }
    y=y*x;
  }

  if(negn)
    y=1/y;

  return y;
}

#include <math.h>

/* \fcnfh
   Compares up to the 'prec' most significative digits

   @returns true if both numbers are the same to the required accuracy
*/
_Bool
fixedcmp(double d1,
	 double d2,
	 int prec)
{
  if(prec>8){
    fprintf(stderr,
	    "fixedcmp:: Sorry, but requested precision can be 8 at"
	    " most. Not %i. STOPPING\n"
	    ,prec);
    exit(EXIT_FAILURE);
  }

  double l10=log(d1)/log(10);
  int expv;

  double prec10=powi(10,prec);
  d1*=prec10;
  d2*=prec10;

  if(l10<0)
    expv=-l10+0.999999999999;
  else
    expv=-l10;

  prec10=powi(10,expv);
  d1*=prec10;
  d2*=prec10;

  l10=d1-d2;

  return l10*l10<1.0;
}
