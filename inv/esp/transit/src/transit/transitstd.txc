/*
 * transitstd.c
 * transitstd.txc - Common routines for the Transit program,
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

#include <transit.h>

/* keeps tracks of number of errors that where allowed to continue. */
static int terr_allown=0;
int transit_nowarn=0;
int verblevel;
int maxline=200;

inline void transitdot(int thislevel, int verblevel)
{
  if(thislevel<=verblevel)
    fwrite(".",1,1,stderr);
}


int
transiterror (int flags,
	      const char *str,
	      ...)
{
  va_list ap;

  va_start(ap,str);
  int ret=vtransiterror(flags,str,ap);
  va_end(ap);

  return ret;
}


/*\fcnfh
  transiterror: Error function for Transit package.

  @returns Number of characters wrote to the standard error file
             descriptor if PERR\_ALLOWCONT is set, otherwise, it ends
             execution of program.
	   0 if it is a warning call and 'transit\_nowarn' is 1
*/
int vtransiterror(int flags, const char *str, va_list ap)
{
  char pre_error[]="\nTransit:: ";
  char error[7][20]={"",
		     "CRITICAL:: ",         /* Produced by the code */
		     "SERIOUS:: ",          /* Produced by the user */
		     "Warning:: ",
		     "Not implemented",
		     "Not implemented",
		     "Not implemented"
  };
  char *errormessage,*out;
  int len,lenout,xtr;

  if(transit_nowarn&&(flags & TERR_NOFLAGBITS)==TERR_WARNING)
    return 0;

  len=0;
  if(!(flags&TERR_NOPREAMBLE))
    len=strlen(pre_error)+strlen(error[flags&TERR_NOFLAGBITS]);
  len+=strlen(str)+1;
  lenout=len;

  errormessage=(char *)calloc(len,sizeof(char));
  out=(char *)calloc(lenout,sizeof(char));

  if(!(flags&TERR_NOPREAMBLE)){
    strcat(errormessage,pre_error);
    strcat(errormessage,error[flags&TERR_NOFLAGBITS]);
  }
  strcat(errormessage,str);

  xtr=vsnprintf(out,lenout,errormessage,ap)+1;

  if(xtr>lenout){
    out=(char *)realloc(out,xtr+1);
    xtr=vsprintf(out,errormessage,ap)+1;
  }

  fwrite(out,sizeof(char),xtr-1,stderr);

  if (flags&TERR_ALLOWCONT||(flags&TERR_NOFLAGBITS)==TERR_WARNING){
    terr_allown++;
    return xtr;
  }

  exit(EXIT_FAILURE);
}


/*\fcnfh
  Check whether the \vr{in} file exists and is openable. If so, return
  the opened file pointer \vr{fp}. Otherwise a NULL is returned and a
  status of why the opening failed is returned.

  @returns 1 on success in opening
           0 if no file was given
           -1 File doesn't exist
           -2 File is not of a valid kind (it is a dir or device)
	   -3 File is not openable (permissions?)
	   -4 Some error happened, stat returned -1
*/
int fileexistopen(char *in,	/* Input filename */
		  FILE **fp)	/* Opened file pointer if successful */
{
  struct stat st;
  *fp=NULL;

  if(in){
    //Check whether the suggested file exists, if it doesn't, then use
    //defaults.
    if (stat(in, &st) == -1){
      if(errno == ENOENT)
	return -1;
      else
	return -4;
    }
    //Not of the valid type
    else if(!(S_ISREG(st.st_mode)||S_ISFIFO(st.st_mode)))
      return -2;
    //Not openable
    else if(((*fp)=fopen(in,"r"))==NULL)
      return -3;
    //No problem!
    return 1;
  }

  //No file was requested
  return 0;

}

/*\fcnfh
  Output for the different cases. of fileexistopen()

  @return 1 on success
          -1 on error (doesn't always returns though
*/
int
verbfileopen(char *in,		/* Input filename */
	     FILE **fp,		/* Opened file pointer if successful */
	     char *desc)	/* Comment on the kind of file */
{
  switch(fileexistopen(in,fp)){
    //Success in opening or user don't want to use atmosphere file
  case 1:
    return 1;
  case 0:
    return -1;
    //File doesn't exist
  case -1:
    transiterror(TERR_SERIOUS,
		 "%sinfo file '%s' doesn't exist."
		 ,desc,in);
    return -1;
    //Filetype not valid
  case -2:
    transiterror(TERR_SERIOUS,
		 "%sfile '%s' is not of a valid kind\n"
		 "(it is a dir or device)\n"
		 ,desc,in);
    return -1;
    //file not openable.
  case -3:
    transiterror(TERR_SERIOUS,
		 "%sfile '%s' is not openable.\n"
		 "Probably because of permissions.\n"
		 ,desc,in);
    return -1;
    //stat returned -1.
  case -4:
    transiterror(TERR_SERIOUS,
		 "Some error happened for %sfile '%s',\n"
		 "stat() returned -1, but file exists\n"
		 ,desc,in);
    return -1;
  default:
    transiterror(TERR_SERIOUS,
		 "Ooops, something weird in file %s, line %i\n"
		 __FILE__,__LINE__);
  }
  return 1;
}


void transitcheckcalled(const long pi, /* Progress indicator variable */
			const char *fcn, /* Name of function being called
					  */ 
			const int n,	/* Number of functions that have
					   to be called before */
			...)	/* Pairs of function required to be
				   called before and their appropiate
				   TRPI\_flag */
{
  /* TD: make 'mess' dynamic size */
  //'ap' is the variable argument pointer.
  //'i' is an auxiliary variable.
  //'mess' is where the output message will be stored.
  //'name' and 'flag' is the function's name and corresponding
  //flag.
  //'cfl' is the cumulative flag of functions required.
  va_list ap;
  int i;
  char mess[1000],*name[n];
  long flag[n],cfl=0;

  //read all arguments
  va_start(ap,n);
  for(i=0;i<n;i++){
    name[i]=(char *)va_arg(ap,char *);
    flag[i]=(long)va_arg(ap,long);
    cfl|=flag[i];
  }
  va_end(ap);

  //If the functions were not called
  if(!pi&cfl){
    //start building 'mess'. Just for grammar, if it is only one function
    if(n==1){
      strcpy(mess,name[i]);
      strcat(mess,"() was not called before function %s()\n");
    }
    //otherwise
    else{ 
      strcpy(mess,"Either ");
      //for each of the functions but the last
      for(i=0;i<n-1;i++){
	strcat(mess,name[i]);
	if(pi&flag[i]) strcat(mess,"(called), ");
	else strcat(mess,"(not-called), ");
      }
      //add `or' for the last one!, of course:D
      strcat(mess,"or ");
      strcat(mess,name[i]);
      if(pi&flag[i]) strcat(mess,"(called) ");
      else strcat(mess,"(not-called) ");
      strcat(mess,"were not called before function %s()\n");
    }
    //output the error!
    transiterror(TERR_CRITICAL,mess,fcn);
  }
}


void 
error (int exitstatus,
       int something, 
       const char *fmt,
       ...)
{
  va_list ap;

  va_start(ap,fmt);
  vtransiterror(TERR_CRITICAL,fmt,ap);
  va_end(ap);

  exit(exitstatus);
}

