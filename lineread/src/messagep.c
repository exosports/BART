/*
 * messages.c   - Common messaging routines. Component
 *                  of the lineread program.
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
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *
 */

#include <messagep.h>

/* Keeps tracks of number of errors that where allowed to continue. */
static int msgp_allown=0;
int msgp_nowarn=0;
int verblevel;
int maxline=1000;
static char *prgname=NULL;

/* Set the value of prgname   */
void
messagep_name(char *name){
  if(prgname) free(prgname);
  prgname = (char *)calloc(strlen(name)+1, sizeof(char *));
  strcpy(prgname, name);
}

/* Free prgname variable      */
void 
messagep_free(){
  free(prgname);
}

int
mperror_fcn(int flags,
            const char *file,
            const long line,
            const char *str,
            ...){
  va_list ap;

  va_start(ap,str);
  int ret = vmperror_fcn(flags, file, line, str, ap);
  va_end(ap);

  return ret;
}


/*\fcnfh
  mperror: Error function for Msgp package.

  @returns Number of characters wrote to the standard error file
             descriptor if PERR\_ALLOWCONT is set, otherwise, it ends
             execution of program.
             0 if it is a warning call and 'msgp\_nowarn' is 1           */
int 
vmperror_fcn(int flags, 
                 const char *file,
                 const long line,
                 const char *str,
                 va_list ap){
  char error[7][22]={"",
                     ":: SYSTEM: ",        /* Produced by the code */
                     ":: USER: ",          /* Produced by the user */
                     ":: Warning: ",
                     ":: Not implemented",
                     ":: Not implemented",
                     ":: Not implemented"
  };
  char *errormessage,*out;
  int len, lenout, xtr;

  if (((flags & MSGP_NOFLAGBITS) == MSGP_SYSTEM) && 
      !(flags & MSGP_NODBG))
    flags |= MSGP_DBG;
  if (flags & MSGP_NODBG)
    flags &= !MSGP_DBG;

  if(msgp_nowarn&&(flags & MSGP_NOFLAGBITS)==MSGP_WARNING)
    return 0;

  len = strlen(prgname)+1;
  if(!(flags&MSGP_NOPREAMBLE))
    len += strlen(error[flags&MSGP_NOFLAGBITS]);
  //symbols + digits + file
  int debugchars = 0;
  if(flags&MSGP_DBG) debugchars = 5 + 6 + strlen(file);
  len    += strlen(str) + 1 + debugchars;
  lenout  = len;

  errormessage = (char *)calloc(len,    sizeof(char));
  out          = (char *)calloc(lenout, sizeof(char));

  strcat(errormessage,"\n");
  strcat(errormessage,prgname);
  if(flags&MSGP_DBG){
    char debugprint[debugchars];
    sprintf(debugprint," (%s|%li)", file, line);
    strcat(errormessage, debugprint);
  }

  if(!(flags&MSGP_NOPREAMBLE)){
    strcat(errormessage, error[flags&MSGP_NOFLAGBITS]);
  }
  strcat(errormessage,str);


  va_list aq;
  va_copy(aq, ap);
  xtr=vsnprintf(out, lenout, errormessage, ap)+1;
  va_end(ap);

  if(xtr>lenout){
    out=(char *)realloc(out,xtr+1);
    xtr=vsnprintf(out, xtr+1, errormessage, aq)+1;
  }
  va_end(aq);
  free(errormessage);

  fwrite(out, sizeof(char), xtr-1, stderr);
  free(out);

  if (flags&MSGP_ALLOWCONT||(flags&MSGP_NOFLAGBITS)==MSGP_WARNING){
    msgp_allown++;
    return xtr;
  }

  messagep_free();
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
           -4 Some error happened, stat returned -1                 */
int fileexistopen(char *in,    /* Input filename                    */
                  FILE **fp){  /* Opened file pointer if successful */
  struct stat st;
  *fp=NULL;

  if(in){
    /* Check whether the suggested file exists, if it doesn't, 
       then use defaults                                        */
    if (stat(in, &st) == -1){
      if(errno == ENOENT)
        return -1;
      else
        return -4;
    }
    /* Invalid type:            */
    else if( !(S_ISREG(st.st_mode) || S_ISFIFO(st.st_mode)) )
      return -2;
    /* Not openable:           */
    else if(((*fp)=fopen(in,"r"))==NULL)
      return -3;
    return 1;
  }
  /* No file was requested:    */
  return 0;

}

/*\fcnfh
  Output for the different cases. of fileexistopen()

  @return fp of opened file on success
          NULL on error (doesn't always returns though   */
FILE *
verbfileopen(char *in,          /* Input filename               */
             char *desc){       /* Comment on the kind of file  */
  FILE *fp;

  switch(fileexistopen(in, &fp)){
    /* Success in opening or user don't want to use atmosphere file */
  case 1:
    return fp;
  case 0:
    mperror(MSGP_USER, "No file was given to open\n");
    return NULL;
    /* File doesn't exist */
  case -1:
    mperror(MSGP_USER, "%s info file '%s' doesn't exist.", desc, in);
    return NULL;
    /* Filetype not valid */
  case -2:
    mperror(MSGP_USER, "%sfile '%s' is not of a valid kind (it is a dir or "
                       "device)\n", desc, in);
    return NULL;
    /* File not openable  */
  case -3:
    mperror(MSGP_USER, "%sfile '%s' is not openable.  Probably because of "
                       "permissions.\n", desc, in);
    return NULL;
    /* Stat returned -1   */
  case -4:
    mperror(MSGP_USER, "Some error happened for %sfile '%s', stat() returned "
                       "-1, but file exists\n", desc, in);
    return NULL;
  default:
    mperror(MSGP_USER, "Ooops, something weird in file %s, line %i\n",
                        __FILE__, __LINE__);
  }
  return NULL;
}


/* \fcnfh
   This function is called if a line of 'file' was longer than 'max'
   characters                                                         */
void 
linetoolong(int max,     /* Maxiumum length of an accepted line */
            char *file,  /* File from which we were reading     */
            int line)    /* Line who was being read             */
{
  mperror(MSGP_USER|MSGP_ALLOWCONT,
          "Line %i of file '%s' has more than %i characters, that is not "
          "allowed\n", file, max);
  exit(EXIT_FAILURE);
}



