/*
 * readdump.c - Read dumped arrays from gdb for the Transit program
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

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/* \fcnfh
   Read dumped data from gdb
*/
int main(int argc, 
	 char **argv)
{
  if(argc<2){
    fprintf(stderr,
	    "Syntax:\n"
	    "%s [<x_data_file>] <y_data_file> [<y_data_file2> ...]\n\n"
	    "This program reads dumped data from 'gdb' into ascii columns\n"
	    ,argv[0]);
    exit(1);
  }
  _Bool yaxis=1;
  if(argc==2)
    yaxis=0;

  FILE *xf=fopen(argv[1],"r");
  char *f1=argv[1];
  if(!xf){
    fprintf(stderr,
	    "'%s' could not be opened\n\n"
	    ,argv[1]);
    exit(1);
  }

  struct stat st;

  stat(argv[1],&st);
  float t=st.st_size/8.0;
  int sz1=(int)t;
  if(sz1!=t){
    fprintf(stderr,
	    "File %s doesn't contain an exact number of double\n"
	    " variable elements(%f)\n\n"
	    ,argv[1],t);
    exit(1);
  }    

  argc-=2;
  argv+=2;
  FILE *(yf[argc]);
  int i;
  int sz2;
  if(yaxis){
    for(i=0;i<argc;i++){
      yf[i]=fopen(*argv,"r");

      if(!yf[i]){
	fprintf(stderr,
		"'%s' could not be opened\n\n"
		,argv[0]);
	exit(1);
      }

      stat(*argv,&st);
      t=st.st_size/8.0;
      sz2=(int)t;
      if(sz2!=t){
	fprintf(stderr,
		"File %s doesn't contain an exact number of double\n"
		" variable elements(%f)\n\n"
		,*argv,t);
	exit(1);
      }
      if(sz1!=sz2){
	fprintf(stderr,
		"Files %s and %s don't contain the same number of elements\n"
		" (%i vs %i)\n"
		,f1,argv[0],sz1,sz2);
	exit(1);
      }
      argv++;
    }

  }

  double x[sz1];
  fread(x,8,sz1,xf);

  if(yaxis){
    double y[argc][sz1];

    for(i=0;i<argc;i++)
      fread(y[i],8,sz1,yf[i]);
    for(sz2=0;sz2<sz1;sz2++){
      printf("%-15.9g"
	     ,x[sz2]);
      for(i=0;i<argc;i++)
	printf("%-15.9g"
	       ,y[i][sz2]);
      printf("\n");
    }
  }
  else
    for(sz2=0;sz2<sz1;sz2++)
      printf("%-15.9g\n"
	     ,x[sz2]);


  exit(EXIT_SUCCESS);
}
