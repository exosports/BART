/*
 * transitstd.c   - Common routines for the Transit program. Component
 *                  of the Transit program.
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
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

#include <transit.h>

/* keeps tracks of number of errors that where allowed to continue. */
static int terr_allown=0;
int transit_nowarn=0;
int verblevel;
int maxline=1000;

inline void transitdot(int thislevel,
                       int verblevel,
                       ...){
  if(thislevel <= verblevel)
    fwrite(".", 1, 1, stderr);
}

int
transiterror_fcn(int flags,
                 const char *file,
                 const long line,
                 const char *str,
                  ...){
  va_list ap;

  va_start(ap, str);
  int ret = vtransiterror_fcn(flags, file, line, str, ap);
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
int vtransiterror_fcn(int flags, 
                      const char *file,
                      const long line,
                      const char *str,
                      va_list ap){
  char prepre_error[] = 
                     "\n******************************************************";
  char pre_error[] = "\n*** Transit";
  char error[7][22] = {"",
                       " :: SYSTEM ERROR ***\n",  /* Produced by the code */
                       " :: USER ERROR ***\n",    /* Produced by the user */
                       " :: Warning ***\n",
                       " :: Not implemented",
                       " :: Not implemented",
                       " :: Not implemented"
  };
  char *errormessage, *out;
  int len, lenout, xtr;
  char post_error[] = 
                     "******************************************************\n";


  if(transit_nowarn && (flags & TERR_NOFLAGBITS)==TERR_WARNING)
    return 0;

  len = strlen(pre_error);
  if(!(flags & TERR_NOPREAMBLE))
    len += strlen(error[flags & TERR_NOFLAGBITS]) + strlen(post_error) + 
           strlen(prepre_error);
  /* Symbols + digits + file: */
  int debugchars = 0;
  if(flags&TERR_DBG) debugchars = 5 + 6 + strlen(file);
  len    += strlen(str) + 1 + debugchars;
  lenout  = len;

  errormessage = (char *)calloc(len   +10, sizeof(char));
  out          = (char *)calloc(lenout+10, sizeof(char));

  if(!(flags & TERR_NOPREAMBLE))
    strcat(errormessage, prepre_error);
  strcat(errormessage, pre_error);
  if(flags & TERR_DBG){
    char debugprint[debugchars];
    sprintf(debugprint," (%s|%li)", file, line);
    strcat(errormessage, debugprint);
  }

  if(!(flags&TERR_NOPREAMBLE))
    strcat(errormessage, error[flags&TERR_NOFLAGBITS]);
  strcat(errormessage, str);
  if(!(flags&TERR_NOPREAMBLE))
    strcat(errormessage, post_error);

  va_list aq;
  va_copy(aq, ap);
  xtr = vsnprintf(out, lenout, errormessage, ap)+1;
  va_end(ap);

  if(xtr > lenout){
    out = (char *)realloc(out, xtr+1);
    xtr = vsnprintf(out, xtr+1, errormessage, aq)+1;
  }
  va_end(aq);
  free(errormessage);

  fwrite(out, sizeof(char), xtr-1, stderr);
  free(out);

  if (flags&TERR_ALLOWCONT || (flags&TERR_NOFLAGBITS)==TERR_WARNING){
    terr_allown++;
    return xtr;
  }

  exit(EXIT_FAILURE);
}


/*\fcnfh
  Check whether the file 'in' exists and is openable. If so, return
  the opened file pointer 'fp'. Otherwise a NULL is returned and a
  status of why the opening failed is returned.

  Return: 1 on success in opening
          0 if no file was given
         -1 File doesn't exist
         -2 File is not of a valid kind (it is a dir or device)
         -3 File is not openable (permissions?)
         -4 Some error happened, stat returned -1                       */
int
fileexistopen(char *in,    /* Input filename                    */
              FILE **fp){  /* Opened file pointer if successful */
  struct stat st;
  *fp = NULL;

  if(in){
    /* Check if the suggested file exists. If it doesn't, use defaults: */
    if (stat(in, &st) == -1){
      if(errno == ENOENT)
        return -1;
      else
        return -4;
    }
    /* Not of the valid type: */
    else if(!(S_ISREG(st.st_mode) || S_ISFIFO(st.st_mode)))
      return -2;
    /* Not openable: */
    else if(((*fp)=fopen(in,"r")) == NULL)
      return -3;
    /* No problem!: */
    return 1;
  }
  /* No file was requested: */
  return 0;
}


/*\fcnfh
  Output for the different cases. of fileexistopen()

  @return fp of opened file on success
          NULL on error (doesn't always returns though */
FILE *
verbfileopen(char *in,     /* Input filename              */
             char *desc){  /* Comment on the kind of file */
  FILE *fp;

  switch(fileexistopen(in, &fp)){
  /* Success in opening or user don't want to use atmosphere file: */
  case 1:
    return fp;
  case 0:
    transiterror(TERR_SERIOUS, "No file was given to open.\n");
    return NULL;
  /* File doesn't exist: */
  case -1:
    transiterror(TERR_SERIOUS, "%s file '%s' doesn't exist.\n", desc, in);
    return NULL;
  /* Filetype not valid: */
  case -2:
    transiterror(TERR_SERIOUS, "%s file '%s' is not of a valid kind "
                               "(it is a dir or device)\n", desc, in);
    return NULL;
  /* File not openable: */
  case -3:
    transiterror(TERR_SERIOUS, "%s file '%s' is not openable. Probably "
                               "because of permissions.\n", desc, in);
    return NULL;
  /* stat returned -1: */
  case -4:
    transiterror(TERR_SERIOUS,
                 "Error happened for %s file '%s', stat() returned -1, "
                 "but file exists.\n", desc, in);
    return NULL;
  default:
    transiterror(TERR_SERIOUS,
                 "Something weird in file %s, line %i.\n", __FILE__, __LINE__);
  }
  return NULL;
}

/*
    Check that the 'n' functions 'variable_argument' have been called
    before 'fcn'. If any of them have not been yet called, throw an error
    message detailed all functions not yet called.                          */
void
transitcheckcalled(const long pi,   /* Progress indicator variable          */
                   const char *fcn, /* Name of function being checked       */
                   const int n,     /* Number of functions that have to be
                                       called before fcn                    */
                   ...){            /* Function's name and flag pairs       */
  /* TD: make 'mess' dynamic size */
  va_list ap;      /* Variable-argument pointer        */
  int i;           /* Auxiliary for index              */
  char mess[1000], /* Output message                   */
       *name;      /* Function's name                  */
  long flag;       /* Function's flag                  */
  _Bool stop=0;    /* Any function has not been called */

  *mess = '\0';
  va_start(ap, n);
  /* Check the pair of arguments: */
  for(i=0; i<n; i++){
    name = (char *)va_arg(ap, char *);
    flag = (long)va_arg(ap, long);
    /* Check if it was called before: */
    if(!(pi&flag)){
      /* Append text to mess: */
      if(!stop){ /* Append only once at first function not called: */
        strcpy(mess, "The following function(s) were not executed before "
                     "this execution of '%s()':\n");
      }
      strcat(mess, "  ");
      strcat(mess, name);
      strcat(mess, "()\n");
      stop = 1;
    }
  }
  va_end(ap);
  /* Print out the error: */
  if(stop)
    transiterror(TERR_CRITICAL, mess, fcn);
}


/* \fcnfh
   Called by a gsl_error */
void 
error(int exitstatus,
      int something, 
      const char *fmt,
      ...){
  va_list ap;
  int len = strlen(fmt);
  char out[len+2];
  strcpy(out, fmt);
  out[len  ] = '\n';
  out[len+1] = '\0';

  va_start(ap, fmt);
  vtransiterror(TERR_CRITICAL, out, ap);
  va_end(ap);

  exit(exitstatus);
}

void
freemem_molecules(struct molecules *mol, long *pi){
  /* Free structures: */
  for(int i=0; i<mol->nmol;  i++)
    free_mol(mol->molec+i);
  /* Free arrays:     */
  free(mol->name[0]);
  free(mol->name);
  free(mol->molec);
  free(mol->mass);
  free(mol->radius);
  /* FINDME: Define a pi for molec */
}

/* \fcnfh
   Frees array in prop_isov, this should be called only once for all the
   isotopes.  */
void
free_isov(prop_isov *isov){
  free(isov->z);
  free(isov->c);
}


/* \fcnfh
   Frees array in prop_isof, this should be called for each of the
   isotopes.                                                       */
void
free_isof(prop_isof *isof){
  free(isof->n);
}


void
free_mol(prop_mol *molec){
  free(molec->d);
  free(molec->q);
}


/* \fcnfh
   Frees array in prop_db, this should be called for each of the
   isotopes.                                                       */
void
free_db(prop_db *db){
  free(db->n);
}


/* \fcnfh
   Frees array in prop_dbnoext, this should be called once per each
   database.                                                       */
void
free_dbnoext(prop_dbnoext *db){
  free(db->T);
}

/* \fcnfh
   Frees array in prop_samp */
void
free_samp(prop_samp *samp){
  free(samp->v);
}


/* \fcnfh
   Frees array in prop_atm */
void
free_atm(prop_atm *atm){
  free(atm->p);
  free(atm->t);
  free(atm->mm);
}

/* \fcnfh
   Saves a string in binary in open file */
void
savestr(FILE *out,
        char *str){
  long len = strlen(str)+1;

  fwrite(&len, sizeof(long), 1, out);
  fwrite(str, 1, len, out);
}

/* \fcnfh
   Restores a string from a binary open file */
int
reststr(FILE *in,
        char **str){
  long len;

  if(fwrite(&len, sizeof(long), 1, in) != 1)
    return -1;
  if(len<0)
    return -2;
  if(len>1000)
    return 1;
  if((*str=(char *)calloc(len, sizeof(char)))==NULL)
    return -3;
  if(fread(*str, 1, len, in) != len)
    return -1;

  return 0;
}


/* \fcnfh
   This function is called if a line of 'file' was longer than 'max'
   characters
*/
void 
linetoolong(int max,     /* Maxiumum length of an accepted line */ 
            char *file,  /* File from which we were reading     */
            int line){   /* Line who was being read             */
  transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
               "Line %i of file '%s' has more than %i characters, "
               "that is not allowed.\n", file, max);
  exit(EXIT_FAILURE);
}


/* Print to screen the time elapsed since time t0
   Return: current time in seconds                       */
double
timecheck(int verblevel,      /* Verbosity level         */
          long iter,          /* Iteration index         */
          long index,         /* Sequencial index        */
          char *str,          /* Time stamp description  */
          struct timeval tv,  /* timeval structure       */
          double t0){         /* Time in seconds         */
  /* Get current time:          */
  gettimeofday(&tv, NULL);
  /* Calculate time in seconds: */
  double sec = tv.tv_sec + 1e-6*tv.tv_usec;
  /* Print time stamp:          */
  transitprint(1, verblevel, "Check point: %02li - %02li %s:  dt = %.4f "
                             "sec.\n\n", iter, index, str, sec-t0);
  return sec;
}
