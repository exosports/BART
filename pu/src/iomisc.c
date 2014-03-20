/*
 * iomisc.c - Miscellaneous input/output utilities.
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


/* List of functions defined in this file:
int    ncharchg (char *str, char car, char chg)
int    nchar (char *str, char car)
char  *readstr_sp_alloc (char *line, char **next, char fspace)
void   freetoolongerr ()
void   settoolongerr (void (*errfcn)(int, char *, int), char *filename,
                      long *currline)
inline
 char fgetupto_err (char *line, int max, FILE *fp,
                    void (*errfcn)(int, char *, int), char *name, long curr)
inline
  char fgetupto (char *line, int max, FILE *fp)
int    getad (int n, char sep, char *str, double **array)
int    getnd (int n, char sep, char *str, ...)
int    getnl (int n, char sep, char *str, ...)
void   fprintpad (FILE *fp, int indent, char *fmt, ...)
double readds (FILE *fp, char *c, char *string, int maxstring)
double getds (char *in, char *c, char *string, int maxstring) 
long   readl (FILE *fp, char *c)
char  *linepad (char *out, int nc, char *in)
double askforposd (char *fmt, ...)
long   askforposl (char *fmt, ...)
char  *fgets_alloc (FILE *fp, int *max)
void   splitnzero_add (char ***array, char *string, char sep)
char **splitnzero_alloc (char *string, char sep) 
void   splitnzero_free (char **multi)
long   countfields (char *l, char sep)
*/

#include <pu/iomisc.h>

char *linepad_break="-";
int printpad_columns=__PADTO_COLUMNS;
char printpad_endpadchar=' ';


/* \fcnfh
   Count ocurrences of 'car' in 'str', replacing each of them by 'chg' */
int
ncharchg(char *str,
         char car,
         char chg){
  int n=0;
  while(*str)
    if(*str++==car){
      str[-1]=chg;
      n++;
    }
  return n;
}


/* \fcnfh
   Count ocurrences of 'car' in 'str'  */
int
nchar(char *str,
      char car){
  int n=0;
  while(*str)
    if(*str++==car)
      n++;
  return n;
}


/* \fcnfh
   Returns next string in a newly allocated array. 'fspace', if any,
   indicates a false space character, which is going to be replaced
   before returning.

   Return: string in a newly allocated array, else
           NULL if allocation was not possible                               */
char *
readstr_sp_alloc(char *line,   /* Input array                                */
                 char **next,  /* pointer to char after the last one used    */
                 char fspace){ /* Character to be replaced by the space char */
  char *lp;
  static char *ret;
  int rn;

  while(*line==' ')  /* Skip empty chars                          */
    line++;
  /* Use lp to keep reading, until NULL, empty or tab char found: */
  lp = line;
  while(*lp && *lp!=' ' && *lp!='\t')
    lp++;

  rn = lp-line; /* Length of string    */
  if((ret=(char *)calloc(rn+1, sizeof(char)))==NULL)
    return NULL;
  /* Set next if not NULL input:       */
  if(next)
    *next = lp;
  /* Copy result into ret and lp:      */
  strncpy(lp=ret, line, rn);
  ret[rn]='\0';
  /* Replace 'fspace' by blank spaces: */
  if(fspace)
    while(*lp)
      if(*lp++ == fspace)
        lp[-1]=' ';

  return ret;
}

static void (*fgut_errfcn)(int, char *, int) = NULL;
static char *fgut_file = NULL;
static long *fgut_currline = NULL;


/* \fcnfh
   Free memory for fgetupto_se() */
void
freetoolongerr(){
  if(fgut_file)
    free(fgut_file);
}


/* \fcnfh
   Store error call for fgetupto_se() so that they don't need to be
   specified each time                                               */
void
settoolongerr(void (*errfcn)(int, char *, int),
              char *filename,
              long *currline){
  fgut_errfcn = errfcn;
  if(fgut_file)
    free(fgut_file);
  fgut_file = calloc(strlen(filename)+1,sizeof(char));
  strcpy(fgut_file, filename);
  fgut_currline = currline;
}


/* \fcnfh
   Reads characters from 'fp' into 'line' until 'max' characters have been 
   read, or newline or end-of-file is reached.  Strips off the newline
   character in 'line' if exists.  If no newline nor end-of-file
   is reached before 'max' chars have been read, call errfcn.

   @returns: 0, if end of file is reached while no characters were read.
             1, if no error function was given.
             The first character read, otherwise.                          */
inline char
fgetupto_err(char *line,  /* Where input is stored                         */
             int max,     /* Maximum number of characters to read          */
             FILE *fp,    /* File pointer to read from                     */
             void (*errfcn)(int, char *, int),
                          /* Error funtion to call in case a newline is not
                             the last character read and we reached 'max'  */
             char *name,  /* Second argument for errfcn                    */
             long curr){  /* Third argument for errfcn                     */
  char *lp;
  int n=0;      /* Size of output (does not include `\\0')  */

  if((lp=fgets(line, max, fp))==NULL)
    return 0;   /* End-of-line case   */
  while(*lp++)  /* Get length of char */
    n++;

  if(n==0){
    fprintf(stderr,
            "This should had never happen, maybe fgets() is not working?. "
            "File %s, line %i.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  /* Strip off \n character:        */
  if(line[n-1] == '\n'){
    line[--n] = '\0';
    if(!n) /* Zero-length output    */
      return '\n';
    return *line;
  }

  /* End of line without a newline: */
  else if(n<max-1)
    return *line;

  /* Line too long:                 */
  if(errfcn)
    errfcn(max, name, curr);
  return -1;
}


/* \fcnfh
   Wrapper for fgetupto() with error parameters previously stored */
inline char
fgetupto(char *line,
         int max,
         FILE *fp){
  long zero = 0;
  if (!fgut_currline) fgut_currline = &zero;
  return fgetupto_err(line, max, fp, fgut_errfcn, fgut_file,
                      *fgut_currline);
}


/* \fcnfh
   get 'n' (all of them if 0) double format numbers from 'str',
   they are separated by 'sep'. Allocate and store them in 'array'

   @returns number of converted values
            -x if there were 'x' fields separated by car instead of 'n'. */
int
getad(int n,          /* Number of float field to accept */
      char sep,       /* Separation character            */
      char *str,      /* String where to look for values */
      double **array) /* Pointers to double value        */
{
  char *ptr=str;
  int rn=1;

  switch(sep){
  default:
    if(!isdigit(sep)) break;
  case 'e': case 'E': case '.': case '-': case '+':
    fprintf(stderr,
            "iomisc:: getnd: invalid separator '%c'. It cannot be a\n"
            "digit, '+', '-', 'e', 'E', nor '.'\n"
            ,sep);
    exit(EXIT_FAILURE);
  }

  do{
    if(*ptr==sep)
      rn++;
  }while(*ptr++);

  if(ptr-1==str)
    return 0;
  if (n&&rn!=n)
    return -rn;
  n = rn;
  *array = (double *)calloc(n, sizeof(double));

  for(rn=0; rn<n; rn++){
    (*array)[rn] = strtod(str, &ptr);

    while(*ptr != sep)
      if(!*ptr++){
        if(rn != n-1)
          fprintf(stderr,
                  "iomisc:: getnd: the one with %f was not supposed to\n"
                  "be the last(%ith) field but it is in: %s\n",
                  (*array)[rn], n, str);
        else
          break;
      }
    str = ptr+1;
  }

  return rn;
}


/* \fcnfh
   Get 'n' double-format numbers from 'str' delimited  by 'sep'. 
   There is no way to know which of the 'n' were really converted. 
   The pointers to the non-converted doubles are returned untouched.

   @Returns number of converted values. else -x if there were only
    x (<n) values in 'str'.                                          */
int
getnd(int n,     /* Number of float field to accept                */
      char sep,  /* Separation character.  If it is a space, also
                    accept tabs as separation character.           */
      char *str, /* String where to look for values                */
      ...){      /* Pointers to double values returned             */
  char *ptr = str; /* Pointer to first character in str       */
  int rn=1,        /* Number of 'sep' separated values in str */ 
      nc=0;        /* Number of read values                   */
  double out, *in;
  va_list ap;

  switch(sep){
  default:
    if(!isdigit(sep)) break;
  case 'e': case 'E': case '.': case '-': case '+':
    fprintf(stderr, "iomisc:: getnd: Invalid separator '%c'. It cannot be a"
            "digit, '+', '-', 'e', 'E', nor '.'\n", sep);
    exit(EXIT_FAILURE);
  }

  /* Count the number of 'sep' separated strings in 'str':  */
  do{
    if(*ptr==sep || (sep==' ' && *ptr=='\t'))
      rn++;
  }while(*ptr++);
  /* Empty str string:                                      */
  if(ptr-1 == str)
    return 0;
  /* There were less 'sep'-separated values than requested: */
  if (rn<n)
    return -rn;

  /* Read the values:                               */
  va_start(ap, str);
  for(rn=0; rn<n; rn++){
    /* Get pointer to store a value in:             */
    in  = (double *)va_arg(ap, double *);
    /* Read value from string and conver to double: */
    out = strtod(str, &ptr);
    if(str != ptr){  /* If it was correctly read,   */
      nc++;          /* count that a value was read */
      *in = out;     /* Assign it to the pointer    */
    }
    /* FINDME: I don't understand what's going on here */
    while(!(*ptr==sep || (sep==' ' && *ptr=='\t')))
      if(!*ptr++)
        if(rn == n-1) /* If this is the last one,  */
          break;      /* break out of the for loop */
    str = ptr+1;
  }
  va_end(ap);

  return nc;
}


/* \fcnfh
   Split the first 'n' strings of 'str' using 'sep' as delimiter.  Format the
   substrings as long integers if possible.

   @returns: number of converted values
            -x if there were 'x' fields separated by car instead of 'n'.     */
int
getnl(int n,     /* Number of long fields to read                            */
      char sep,  /* Separator character. If it is a space, also accept tabs  */
      char *str, /* Input string                                             */
      ...){      /* Pointers to long values                                  */

  char *ptr=str;   /* Pointer to the input string         */
  int rn=1,        /* Number of sub-strings               */
      nc=0;        /* Number or stored values             */
  long out, *in;   /* */
  va_list ap;      /* Object to handle variable arguments */

  /* Complain if sep is a digit, '+', or '-' sign.        */
  switch(sep){
  default:
    if(!isdigit(sep)) break;
  case '-':
  case '+':
    fprintf(stderr,
            "iomisc:: getnd: invalid separator '%c'. It cannot be a "
            "digit, '+', nor '-'.\n", sep);
    exit(EXIT_FAILURE);
  }

  /* Count the number of sub-strings delimited by sep:    */
  do{
    if(*ptr==sep || (sep==' ' && *ptr=='\t')){
      rn++;
      while(*ptr==sep || (sep==' ' && *ptr=='\t'))
        ptr++;
    }
  }while(*ptr++);

  /* Empty imput string:                                  */
  if(ptr-1==str)
    return 0;
  /* rn is not the number of requested numbers:           */
  if (rn!=n)
    return -rn;

  /* Store values in input pointers:                      */
  va_start(ap, str);
  for(rn=0; rn<n; rn++){
    in = (long *)va_arg(ap, long *); /* Get pointer       */
    out = strtol(str, &ptr, 0);      /* Read a number     */
    if(str != ptr){
      nc++;               /* Number of saved longs        */
      *in = out;          /* Set value of the pointer     */
    }
    /* Skip separator chars:                              */
    while(!(*ptr==sep || (sep==' ' && *ptr=='\t')))
      if(!*ptr++){        /* End of string reached:       */
        if(rn != n-1)     /* Stored/requested don't match */
          fprintf(stderr,
                  "iomisc:: getnl: the one with %li was not supposed to "
                  "be the last(%ith) field but it is in: %s.\n", out, n, str);
        else
          break;
      }
    str = ptr+1; /* Move str pointer to next value       */
  }
  va_end(ap);
  return nc;
}

/* \fcnfh
   Prints a parragraph with justification using linepad, it takes a
   indent value                                                      */
void
fprintpad(FILE *fp,    /* File pointer to print                      */
          int indent,  /* Number of blanks to leave at the beginning */ 
          char *fmt,   /* Format                                     */
          ...){        /* Parameters                                 */
  char *out, *ptr;
  int cols;
  int length=80, reall;
  va_list ap;

  out = getenv("COLUMNS");
  if(!out || !(cols=atoi(out)))
    cols = printpad_columns;
  if(indent >= cols)
    fprintf(stderr,
            "iomisc::fprintpad(): Indent(%i) is bigger than columns(%i).\n",
            indent, cols);

  char top[cols+1];

  out = (char *)calloc(length, sizeof(char));
  va_start(ap, fmt);
  reall = vsnprintf(out, length, fmt, ap);
  if(reall >= length){
    out = (char *)realloc(out, (reall+1)*sizeof(char));
    reall = vsnprintf(out, reall+1, fmt, ap);
  }
  va_end(ap);

  ptr = out;
  fprintf(stderr,   "%*s",  indent, "");
  while((ptr=linepad(top, cols-indent, ptr))){
    fprintf(stderr, "%s\n", top);
    fprintf(stderr, "%*s",  indent, "");
  }
  fprintf(stderr, "%s", top);

}


/* \fcnfh
   read a double value and optionally return the rest of the string

   Return: on-success double value read, and -1 on 'c'
            on-failure 0 returned and first character read on 'c' */
double readds(FILE *fp,
              char *c,
              char *string,   /* NULL if no string return is desired */ 
              int maxstring){ /* includes '\\0'                      */
  char *str, *cur, car;
  int salloc=8, ns;
  _Bool atleast=0;
  enum {ent=0x1, dec=0x2, man=0x4, unsig=0x8} stage=ent|unsig;

  cur = str = (char *)calloc(salloc, sizeof(char));

  fflush(fp);
  while(fread(&car, 1, 1, fp)){
    if(cur-str == salloc){
      cur = (str=(char *)realloc(str,( salloc<<1)*sizeof(char))) + salloc;
      salloc <<= 1;
    }
    *cur = car;
    if(car == '\n')
      break;
    if(car=='e' || car=='E'){
      if(!atleast || stage&man)
        break;
      stage = man|unsig;
    }
    else if(car=='+' || car=='-'){
      if(!(stage & unsig))
        break;
      stage &= ~unsig;
    }
    else if(car == '.'){
      if(!stage & ent)
        break;
      stage = dec;
    }
    else if(isdigit(car)){
      atleast = 1;
      stage &= ~unsig;
    }
    else
      break;

    cur++;
  }
  if(string){
    ns = 0;
    /* An optional dash can separate number from string: */
    if(car=='-')
      fread(&car,1,1,fp);
    while(car != '\n'){
      /* If line is bigger than wanted, stop and discard rest of line: */
      if(ns++ >= maxstring)
        break;
      *string++ = car;
      fread(&car, 1, 1, fp);
    }
    *string = '\0';
  }
  while(car!='\n' && fread(&car, 1, 1, fp));

  if(!atleast){
    if(c) *c = *cur;
    return 0;
  }
  *cur = '\0';
  if(c) *c = -1;
  return atof(str);
}


/* \fcnfh
   Read a double value from a field in a string and optionally return the
   rest of the field in 'string'. A field ends in space and can be 
   separated from the number by a dash.
   Return: the read value, on success;
           0, else.                                                           */
double 
getds(char *in,       /* Input field                                          */
      char *c,        /* If not NULL, store a -1 for a valid number field; 
                         else, store the first character of the field         */
      char *string,   /* If not NULL, store the string trailing the value     */
      int maxstring){ /* Maximum number of characters to read                 */

  char *cur,       /* Copy of *in                 */
        car;       /* Pointer to chars in *cur    */
  _Bool atleast=0; /* There is at least one digit */
  enum {ent=0x1, dec=0x2, man=0x4, unsig=0x8} stage = ent|unsig;

  /* Read until *cur is no longer a vaild numeric value: */
  cur = in;
  /* Go over each character: */
  while((car = *cur)){
    if(car == '\n')  /* Stop if we reached the end of line */
      break;
    if(car=='e' || car=='E'){
      if(!atleast || stage&man) /* Stop if no digit read before or */
        break;                  /* if stage is already man         */
      stage = man|unsig;
    }
    else if(car=='+' || car=='-'){
      if(!(stage & unsig)) /* Stop if we already had read a digit */
        break;
      stage &= ~unsig;
    }
    else if(car == '.'){
      if(!stage & ent)     /* */
        break;
      stage=dec;
    }
    else if(isdigit(car)){
      atleast=1;
      stage &= ~unsig;
    }
    else
      break;
    cur++;
  }

  /* Set the return string:  */
  if(string){
    /* Skip dash if present: */
    if(*cur == '-')
      cur++;
    /* read until blank, tab, or termination character: */
    while(*cur!=' ' && *cur!='\t' && *cur!='\0'){
      /* If line is larger than wanted, stop and discard the rest: */
      if(!--maxstring)
        break;
      *string++ = *cur++;
    }
    *string = '\0';
  }

  /* Not a valid number: */
  if(!atleast){
    if(c) *c = *in;
    return 0;
  }
  if(c)
    *c = -1;
  return atof(in);
}


/* \fcnfh
   read a long value

   @returns on-success long value read and -1 on 'c'
            on-failure 0 returned and first character read on 'c'
*/
long readl(FILE *fp,
           char *c)
{
  char *str,*cur,car;
  int salloc=8;
  _Bool atleast=0;
  _Bool unsig=1;

  cur=str=(char *)calloc(salloc,sizeof(char));

  fflush(fp);
  while(fread(&car,1,1,fp)){
    if(cur-str==salloc){
      cur=(str=(char *)realloc(str,(salloc<<1)*sizeof(char)))+salloc;
      salloc<<=1;
    }
    *cur=car;
    if(car=='\n')
      break;
    if(car=='+'||car=='-'){
      if(!unsig)
        break;
      unsig=0;
    }
    else if(isdigit(car)){
      atleast=1;
      unsig=0;
    }
    else
      break;

    cur++;
  }
  while(car!='\n' && fread(&car,1,1,fp));

  if(!atleast){
    if(c) *c=*cur;
    return 0;
  }
  *cur='\0';
  if(c) *c=-1;
  return atol(str);
}


/* \fcnfh
   Align a line by adding spaces in between words, returns a pointer to
   the position right after the the last word. All newlines and tabs are
   cosidered spaces. It doesn't add a newline at the end of lines.

   @returns pointer to the beggining of new sentence
            NULL if all is processed, in this case it doesn't do alignement.
*/
char *
linepad(char *out,                /* output, it has to have a length of,
                                   at least, nc+1 characters
                                   allocated. */ 
        int nc,                        /* number of columns */
        char *in)                /* input array */
{
  //word count 'wc' initializes to zero, it will be counting
  //interspaces.'out' index 'o' also initialize to zero.
  int wc=0,o=0,lasto=0;
  _Bool prev=0;
  out[nc]='\0';

  while(o<nc){
    //\linelabel{lastline}
    if(!*in){
      out[o]='\0';
      return NULL;
    }
    if(*in==' '||*in=='\n'||*in=='\t'){
      if(!prev){
        wc++;
        out[o]=' ';
        lasto=o;
        prev=1;
      }
    }
    else{
      prev=0;
      out[o]=*in;
    }
    in++;
    o++;
  }
  //if there is only one word, which is larger than 'nc' then break with
  //the string 'linepad\_break'.
  if(!lasto){
    o=strlen(linepad_break);
    strcpy(out+nc-o,linepad_break);
    return in-o;
  }
  //'out' is now in the last position before padding.
  out+=lasto-1;
  out[1]='\0';
  //'lasto' now has the number of spaces that need to be added. 'in'
  //points to the beggining of the first non-processed word in the input
  //array
  lasto=nc-lasto;
  in-=lasto-1;
  wc--;

  while(lasto){
    //if two words don't fit in the line, then right align the first
    //word
    if(!wc){
      nc=nc-lasto;
      while(nc--){
        out[lasto]=*out;
        out--;
      }
      while(lasto)
        out[lasto--]=printpad_endpadchar;
      break;
    }
    //if space found, then add 'lasto'/'wc' (+1 if they were not an
    //exact ratio) extra spaces before this word. Diminish 'lasto' by
    //those many spaces and 'wc' by 1.
    if(*out==' '){
      o=(int)((float)lasto/wc--+0.999999999999);
      while(o--)
        out[lasto--]=' ';
    }
    out[lasto]=*out;
    out--;
  }

  return in;
}




/* \fcnfh
   ask for a positive double value from stdin with the question 'fmt'
   from stdin          */
double
askforposd(char *fmt,  /* Question asked                     */
           ...){       /* Parameters for the question if any */

  va_list ap;         /* Handle variable arguments list */
  va_start(ap, fmt);  
  double val;         /* Returned value                 */
  char rc;            /* Returned value from stdin      */

  while(1){
    vfprintf(stderr, fmt, ap); /* Print question */
    val = readd(stdin, &rc);   /* Request answer */
    /* Quit program:                   */
    if(rc=='q'){
      fprintf(stderr, "User interrupt!\n");
      exit(EXIT_SUCCESS);
    }
    /* Ask again if value is negative: */
    if(val<=0)
      fprintf(stderr," Invalid value %g, must be positive.\n", val);
    /* Break if NULL answer:           */
    else if(!rc)
      break;
    fprintf(stderr, "Try again.\n");
  }
  va_end(ap);

  return val;
}


/* \fcnfh
   Ask for a positive long value from stdin with the question 'fmt'
   from stdin  */
long
askforposl(char *fmt, /* Question asked                     */
           ...){      /* Parameters for the question if any */

  va_list ap;        /* Handle variable arguments list */
  va_start(ap, fmt);
  long val;          /* Returned value                 */
  char rc;           /* Returned value from stdin      */

  while(1){
    vfprintf(stderr, fmt, ap); /* Print question */
    val = readl(stdin, &rc);   /* Request answer */
    /* Quit program:                   */
    if(rc == 'q'){
      fprintf(stderr, "User interrupt!\n");
      exit(EXIT_SUCCESS);
    }
    /* Ask again if value is negative: */
    if(val <= 0)
      fprintf(stderr, " Invalid value %li, must be positive.\n", val);
    /* Break if NULL answer:           */
    else if(!rc)
      break;
    fprintf(stderr, "Try again!\n");
  }
  va_end(ap);

  return val;
}


/* \fcnfh
   Read a new line and return it in a newly allocated 

   @returns string newly allocated with with data read from fp
            NULL if end of file was found
 */
char *
fgets_alloc(FILE *fp,                 /* where to read from */
            int *max)                /* Number of characters read if pointer was
                                   given */
{
  static char *str;
  char *rp;
  int oldalloc,alloc=8;

  //initial allocation
  str=(char *)calloc(alloc,sizeof(char));
  rp=str;
  oldalloc=7;

  //go until end of file or end of line.
  while(1){
    //Read the mext part of the line
    if(fgets(rp,oldalloc+1,fp)==NULL&&rp==str){
      free(str);
      return NULL;
    }
    while(*rp) rp++;

    //Check if we reached an end of line
    if(rp[-1]=='\n'){
      rp[-1]='\0';
      if(max)
        *max=rp-str;
      str=(char *)realloc(str,(rp-str)*sizeof(char));
      return str;
    }

    //Now see if the whole line was not read and there was no '\\n' (as
    //checked by the previous conditional, in this case an end of file
    //without a final newline should have ocurred
    if(rp-str+1<alloc){
      if(max)
        *max=rp-str+1;
      str=(char *)realloc(str,(rp-str+1)*sizeof(char));
      return str;
    }

    //reallocate for new lecture. The whole allocated string wasa ready,
    //but the last character is '\\0'
    oldalloc=alloc;
    str=(char *)realloc(str,sizeof(char)*(alloc<<=1));
    rp=str+oldalloc-1;
  }

  //this never happens.
  return NULL;
}


/* \fcnfh
   Add elements to an array already allocated by splitnzero\_alloc. It
   can also initialize if the a pointer to NULL is given.

*/
void
splitnzero_add(char ***array,        /* Pointer to the already allocated
                                   array. If it is null, then allocates
                                   a new array */
               char *string,        /* string to parse */
               char sep)        /* separator */
{
  if(!array){
    fprintf(stderr,
            "iomisc:: there was no pointer address given to"
            " splitnzero_add!\n"
            );
    exit(EXIT_FAILURE);
  }

  if(!*array){
    *array=splitnzero_alloc(string,sep);
    return;
  }

  long n=0,ns=0;
  long len=0,lenold=0;
  char **point=*array;
  char *sp=string;

  //find length and fields of the old array
  while(*point){
    n++;
    lenold+=strlen(*point++)+1;
  }
  //find how many and length of the new strings.
  while(*sp){
    len++;
    if(*sp++==sep)
      ns++;
  }

  //reallocate arrays
  *array=(char **)realloc(*array,(n+ns+1)*sizeof(char *));
  sp=**array=(char *)realloc(**array,(len+lenold+1)*sizeof(char));

  //reassign old string pointers, at the end of the loop, sp points to
  //the first space where the new data is going to go.
  long na=1;
  while(sp-(**array)<lenold){
    if(!*sp)
      array[0][na++]=sp+1;
    sp++;
  }

  //now store the new strings
  while(*string){
    if(*string==sep){
      array[0][na++]=sp+1;
      *sp='\0';
    }
    else
      *sp=*string;
    string++;
    sp++;
  }
  *sp='\0';

  array[0][na]=NULL;

}


/* \fcnfh
   Split strings according to separators, allocating required space fot
   the 2D array. Last element is a null pointer.

   @returns 2D-array with the splitted data. The last element of the
                     returned array is null.
*/
char **
splitnzero_alloc(char *string,        /* String to split */
                 char sep)        /* Separator */
{
  //'nf' starts at 1 because there at least the starting field. 'len' is
  //the string's length until the loop in \lin{looplen}
  static char **multi;
  char *sp=string;
  int nf=1,len=0;

  //if string is null, then return an array of 1 NULL element
  if(!sp){
    multi=(char **)calloc(1,sizeof(char *));
    multi[0]=NULL;
    return multi;
  }

  //count number of fields, add one each time we find a separator.
  while(*sp){
    len++;
    if(*sp++==sep)
      nf++;
  }

  //allocate array to return, allocate nf+1 because of the last NULL
  //element...
  multi=(char **)calloc(nf+1,sizeof(char *));
  *multi=(char *)calloc(len+1,sizeof(char));
  //\linelabel{looplen}
  len=nf=0;
  sp=string;
  while(sp[len]){
    if(sp[len]==sep){
      strncpy(multi[nf],sp,len);
      multi[nf++][len]='\0';
      multi[nf]=multi[nf-1]+len+1;
      sp+=len+1;
      len=-1;
    }
    len++;
  }
  strncpy(multi[nf],sp,len);
  multi[nf][len]='\0';
  multi[++nf]=NULL;

  return multi;
}


/* \fcnfh
   Free space allocated by splitnzero\_alloc()
 */
void
splitnzero_free(char **multi)
{
  if(*multi)
    free(*multi);
  free(multi);
}


/* Count the number of fields in a string lp separated by the character sep
   Return: number of fields  */
long
countfields(char *lp,  /* Input string    */
            char sep){ /* Field separator */
  long nfields = 0;  /* Number of fields in lp */

  /* Skip sep chars at beginning: */
  while(*lp++ == sep){}

  /* Count fields:                      */
  while(*lp){
    /* Go through a field:              */
    while(*lp++ != sep  &&  *lp){}
    nfields++;
    /* Skip sep chars until next field: */
    while(*lp++ == sep){}
  }
  return nfields;
}


/* Reads characters from current position of line until it reaches a blank
   space or line-break, store string in name.                               */
void
getname(char *line,
        char *name){
  /* Copy characters from line into name: */
  while (*line != ' ' && *line != '\n')
    *name++ = *line++;
}

/* Search for string atom in list.
   Return: the index in list if found,
           -1, else                       */
int
findstring(char *string,  /* String being searched      */
           char **list,   /* List of strings            */
           int size){     /* Number of elements in list */
  int i;
  for (i=0; i<size; i++)
    if (strcmp(string, list[i]) == 0)
      return i;
  return -1;
}


/* \fcnfh
   Move pointer 'lp' to next field.
   @returns pointer to the beginning of next non-space                  */
char *
nextfield(char *lp){
  /* Skip leading blank spaces and tabs: */
  while(*lp == ' ' || *lp == '\t') lp++;
  /* Skip field                          */
  while(*lp != ' ' && *lp) lp++;
  /* Skip trailing blank spaces and tabs */
  while(*lp == ' ' || *lp == '\t') lp++;

  return lp;
}
