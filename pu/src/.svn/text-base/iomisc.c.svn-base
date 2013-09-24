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

#include <pu/iomisc.h>

char *linepad_break="-";
int printpad_columns=__PADTO_COLUMNS;
char printpad_endpadchar=' ';


/* \fcnfh
   Count ocurrences of 'car' in 'str', replacing each of them by 'chg'.
*/
int
ncharchg(char *str,
	 char car,
	 char chg)
{
  int n=0;
  while(*str)
    if(*str++==car){
      str[-1]=chg;
      n++;
    }
  return n;
}


/* \fcnfh
   Count ocurrences of 'car' in 'str'.
*/
int
nchar(char *str,
      char car)
{
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

   @returns string in a newly allocated array
            NULL if allocation was not possible
*/
char *
readstr_sp_alloc(char *line,	/* Input array */
		 char **next,	/* pointer to the character after the
				   last one used */
		 char fspace)	/* False space character, any ocurrence
				   of this will be replaced by space
				   before returning */
{
  char *lp;
  static char *ret;
  int rn;

  while(*line==' ')
    line++;
  lp=line;
  while(*lp&&*lp!=' '&&*lp!='\t')
    lp++;

  rn=lp-line;
  if((ret=(char *)calloc(rn+1,sizeof(char)))==NULL)
    return NULL;
  if(next)
    *next=lp;
  strncpy(lp=ret,line,rn);
  ret[rn]='\0';
  //change 'fspace' per spaces
  if(fspace)
    while(*lp)
      if(*lp++==fspace)
	lp[-1]=' ';

  return ret;
}

static void (*fgut_errfcn)(int, char *, int) = NULL;
static char *fgut_file = NULL;
static long *fgut_currline = NULL;

/* \fcnfh
   Free memory for fgetupto_se()
*/
void
freetoolongerr()
{
  if(fgut_file)
    free(fgut_file);
}

/* \fcnfh
   Store error call for fgetupto_se() so that they don't need to be
   specified each time
*/
void
settoolongerr(void (*errfcn)(int, char *, int),
	      char *filename,
	      long *currline)
{
  fgut_errfcn = errfcn;
  if(fgut_file)
    free(fgut_file);
  fgut_file = calloc(strlen(filename)+1,sizeof(char));
  strcpy(fgut_file, filename);
  fgut_currline = currline;
}

/* \fcnfh
   reads into 'line' from 'fp' up to 'max' characters, it checks that the
   last character read is '\\n', if not and it is not the end of file,
   then it assumes error as max was too small and calls errfcn with
   parameters 'max' and 'name'. If '\\n' is found then it takes it out
   from the returned value. the idea behind returning the first
   character is to quickly identify comment lines

   @returns 0 if end of file is reached while no characters were read
            first character read otherwise.
	    -1 if no error function was given.
*/
inline char
fgetupto_err(char *line,		/* Where input is stored, it
					   has to have dimension 'max'
					*/ 
	     int max, 		/* dimension of 'line' */
	     FILE *fp,		/* File pointer to read from */
	     void (*errfcn)(int,char *, int), /* Error funtion to call
						 in case a newline is
						 not the last
						 character read and we
						 reached 'max' */
	     char *name,	/* Second parameter with which errfcn
				   is going to be called */
	     long curr)
{
  char *lp;
  int n=0;

  if((lp=fgets(line,max,fp))==NULL)
    return 0;
  //'n' does not include `\\0'
  while(*lp++)
    n++;

  if(n==0){
    fprintf(stderr,
	    "UYYUYUY!, this should had not happened ever, maybe\n"
	    "fgets() is not working?. File %s, line %i\n"
	    ,__FILE__,__LINE__);
    exit(EXIT_FAILURE);
  }

  //regular line
  if(line[n-1]=='\n'){
    line[--n]='\0';
    if(!n)
      return '\n';
    return *line;
  }
  //End of file without a newline
  else if(n<max-1)
    return *line;

  //Line too large
  if(errfcn)
    errfcn(max,name,curr);
  return -1;
}


/* \fcnfh
   Wrapper for fgetupto() with error parameters previously stored
*/
inline char
fgetupto(char *line,
	 int max,
	 FILE *fp)
{
  long zero=0;
  if (!fgut_currline) fgut_currline = &zero;
  return fgetupto_err(line, max, fp, fgut_errfcn, fgut_file,
		      *fgut_currline);
}


/* \fcnfh
   get 'n' (all of them if 0) double format numbers from 'str',
   they are separated by 'sep'. Allocate and store them in 'array'

   @returns number of converted values
            -x if there were 'x' fields separated by car instead of 'n'.
*/
int
getad(int n,			/* Number of float field to accept */
      char sep,			/* separation character */
      char *str,		/* String where to look for values */
      double **array)		/* Pointers to double value */
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
  n=rn;
  *array=(double *)calloc(n,sizeof(double));

  for(rn=0;rn<n;rn++){
    (*array)[rn]=strtod(str,&ptr);

    while(*ptr!=sep)
      if(!*ptr++){
	if(rn!=n-1)
	  fprintf(stderr,
		  "iomisc:: getnd: the one with %f was not supposed to\n"
		  "be the last(%ith) field but it is in: %s\n"
		  ,(*array)[rn],n,str);
	else
	  break;
      }
    str=ptr+1;
  }

  return rn;
}


/* \fcnfh
   get 'n' double format numbers from 'str', they are separated by
   'sep'. There is no way to know which of the 'n' were really
   converted. The pointers to the non-converted doubles are returned
   untouched.

   @returns number of converted values
            -x if there were only 'x' fields separated by car instead 
	       of 'n'.
*/
int
getnd(int n,			/* Number of float field to accept */
      char sep,			/* separation character. If it is a
				   space, then it also accept tabs as
				   separation character */
      char *str,		/* String where to look for values */
      ...)			/* Pointers to double value */
{
  char *ptr=str;
  int rn=1,nc=0;
  double out,*in;
  va_list ap;

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
    if(*ptr==sep||(sep==' '&&*ptr=='\t'))
      rn++;
  }while(*ptr++);
  if(ptr-1==str)
    return 0;
  if (rn<n)
    return -rn;

  va_start(ap,str);
  for(rn=0;rn<n;rn++){
    in=(double *)va_arg(ap,double *);
    out=strtod(str,&ptr);
    if(str!=ptr){
      nc++;
      *in=out;
    }

    while(!(*ptr==sep||(sep==' '&&*ptr=='\t')))
      if(!*ptr++)
	if(rn==n-1)
	  break;
    str=ptr+1;
  }
  va_end(ap);

  return nc;
}


/* \fcnfh
   get 'n' long format numbers from 'str', they are separated by
   'sep'. There is no way to know which of the 'n' were really
   converted. The non-converted doubles are returned untouched.

   @returns number of converted values
            -x if there were 'x' fields separated by car instead of 'n'.
*/
int
getnl(int n,			/* Number of long fields to accept */
      char sep,			/* separation character. If it is a
				   space, then it also accept tabs as
				   separation character */
      char *str,		/* String where to look for values */
      ...)			/* Pointers to long value */
{
  char *ptr=str;
  int rn=1,nc=0;
  long out,*in;
  va_list ap;

  switch(sep){
  default:
    if(!isdigit(sep)) break;
  case '-': case '+':
    fprintf(stderr,
	    "iomisc:: getnd: invalid separator '%c'. It cannot be a\n"
	    "digit, '+', nor '-'\n"
	    ,sep);
    exit(EXIT_FAILURE);
  }

  do{
    if(*ptr==sep||(sep==' '&&*ptr=='\t')){
      rn++;
      while(*ptr==sep||(sep==' '&&*ptr=='\t'))
	ptr++;
    }
  }while(*ptr++);
  if(ptr-1==str)
    return 0;
  if (rn!=n)
    return -rn;
  va_start(ap,str);
  for(rn=0;rn<n;rn++){
    in=(long *)va_arg(ap,long *);
    out=strtol(str,&ptr,0);
    if(str!=ptr){
      nc++;
      *in=out;
    }

    while(!(*ptr==sep||(sep==' '&&*ptr=='\t')))
      if(!*ptr++){
	if(rn!=n-1)
	  fprintf(stderr,
		  "iomisc:: getnl: the one with %li was not supposed to\n"
		  "be the last(%ith) field but it is in: %s\n"
		  ,out,n,str);
	else
	  break;
      }
    str=ptr+1;
  }
  va_end(ap);
  return nc;
}

/* \fcnfh
   Prints a parragraph with justification using linepad, it takes a indent
   value
*/
void
fprintpad(FILE *fp,		/* File pointer to print */
	  int indent,		/* Number of blanks to leave at the
				   beggining */ 
	  char *fmt,		/* format */
	  ...)			/* parameters */
{
  char *out,*ptr;
  int cols;
  int length=80,reall;
  va_list ap;

  out=getenv("COLUMNS");
  if(!out||!(cols=atoi(out)))
    cols=printpad_columns;
  if(indent>=cols)
    fprintf(stderr,
	    "iomisc::fprintpad(): Indent(%i) is bigger than columns(%i)\n"
	    ,indent,cols);

  char top[cols+1];

  out=(char *)calloc(length,sizeof(char));
  va_start(ap,fmt);
  reall=vsnprintf(out,length,fmt,ap);
  if(reall>=length){
    out=(char *)realloc(out,(reall+1)*sizeof(char));
    reall=vsnprintf(out,reall+1,fmt,ap);
  }
  va_end(ap);

  ptr=out;
  fprintf(stderr,"%*s",indent,"");
  while((ptr=linepad(top,cols-indent,ptr))){
    fprintf(stderr,"%s\n",top);
    fprintf(stderr,"%*s",indent,"");
  }
  fprintf(stderr,"%s",top);

}


/* \fcnfh
   read a double value and optionally return the rest of the string

   @returns on-success double value read, and -1 on 'c'
            on-failure 0 returned and first character read on 'c'
*/
double readds(FILE *fp,
	      char *c,
	      char *string,	/* NULL if no string return is desired
				   */ 
	      int maxstring)	/* includes '\\0' */
{
  char *str,*cur,car;
  int salloc=8,ns;
  _Bool atleast=0;
  enum {ent=0x1,dec=0x2,man=0x4,unsig=0x8} stage=ent|unsig;

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
    if(car=='e'||car=='E'){
      if(!atleast||stage&man)
	break;
      stage=man|unsig;
    }
    else if(car=='+'||car=='-'){
      if(!(stage&unsig))
	break;
      stage&=~unsig;
    }
    else if(car=='.'){
      if(!stage&ent)
	break;
      stage=dec;
    }
    else if(isdigit(car)){
      atleast=1;
      stage&=~unsig;
    }
    else
      break;

    cur++;
  }
  if(string){
    ns=0;
    //an optional dash can separate number from string
    if(car=='-')
      fread(&car,1,1,fp);
    while(car!='\n'){
      //if line is bigger than wanted, stop storing and discard rest of
      //line
      if(ns++>=maxstring)
	break;
      *string++=car;
      fread(&car,1,1,fp);
    }
    *string='\0';
  }
  while(car!='\n' && fread(&car,1,1,fp));

  if(!atleast){
    if(c) *c=*cur;
    return 0;
  }
  *cur='\0';
  if(c) *c=-1;
  return atof(str);
}


/* \fcnfh
   read a double value from a string and optionally return the rest of
   the field in 'string', a field ends in space and can be separated
   from the number by a dash.

   @returns on-success double value read, and -1 on 'c'
            on-failure 0 returned and first character read on 'c'
*/
double getds(char *in,		/* input field */
	     char *c,		/* 0 if there is a valid number at the
				   beggining of the field, otherwie
				   return here the first character of
				   the field. Optional  */
	     char *string,	/* return string here. Optional: NULL if
				   no return is desired */
	     int maxstring)	/* includes '\\0' */
{
  char *cur,car;
  _Bool atleast=0;
  enum {ent=0x1,dec=0x2,man=0x4,unsig=0x8} stage=ent|unsig;

  cur=in;
  while((car=*cur)){
    if(car=='\n')
      break;
    if(car=='e'||car=='E'){
      if(!atleast||stage&man)
	break;
      stage=man|unsig;
    }
    else if(car=='+'||car=='-'){
      if(!(stage&unsig))
	break;
      stage&=~unsig;
    }
    else if(car=='.'){
      if(!stage&ent)
	break;
      stage=dec;
    }
    else if(isdigit(car)){
      atleast=1;
      stage&=~unsig;
    }
    else
      break;

    cur++;
  }
  if(string){
    //an optional dash can separate number from string
    if(*cur=='-')
      cur++;
    while(*cur!=' '&&*cur!='\t'&&*cur!='\0'){
      //if line is bigger than wanted, stop storing and discard rest of
      //line
      if(!--maxstring)
	break;
      *string++=*cur++;
    }
    *string='\0';
  }

  if(!atleast){
    if(c) *c=*in;
    return 0;
  }
  if(c) *c=-1;
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
linepad(char *out,		/* output, it has to have a length of,
				   at least, nc+1 characters
				   allocated. */ 
	int nc,			/* number of columns */
	char *in)		/* input array */
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
   ask for a positive double value from stdin with the question 'fmt';
   to stderr
*/
double
askforposd(char *fmt,		/* Question asked */
	   ...)			/* parameters for the question if any */
{
  va_list ap;
  va_start(ap,fmt);
  double val;
  char rc;

  while(1){
    vfprintf(stderr,fmt,ap);
    val=readd(stdin,&rc);
    if(rc=='q'){
      fprintf(stderr,"User interrupt!\n");
      exit(EXIT_SUCCESS);
    }
    if(val<=0)
      fprintf(stderr," Invalid value %g, it has to be positive\n",val);
    else if(!rc)
      break;
    fprintf(stderr,"Try again!\n");
  }
  va_end(ap);

  return val;
}


/* \fcnfh
   ask for a positive long value from stdin with the question 'fmt';
   to stderr
*/
long
askforposl(char *fmt,		/* Question asked */
	   ...)			/* parameters for the question if any */
{
  va_list ap;
  va_start(ap,fmt);
  long val;
  char rc;

  while(1){
    vfprintf(stderr,fmt,ap);
    val=readl(stdin,&rc);
    if(rc=='q'){
      fprintf(stderr,"User interrupt!\n");
      exit(EXIT_SUCCESS);
    }
    if(val<=0)
      fprintf(stderr," Invalid value %li, it has to be positive\n",val);
    else if(!rc)
      break;
    fprintf(stderr,"Try again!\n");
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
fgets_alloc(FILE *fp, 		/* where to read from */
	    int *max)		/* Number of characters read if pointer was
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
splitnzero_add(char ***array,	/* Pointer to the already allocated
				   array. If it is null, then allocates
				   a new array */
	       char *string,	/* string to parse */
	       char sep)	/* separator */
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
splitnzero_alloc(char *string,	/* String to split */
		 char sep)	/* Separator */
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

/* \fcnfh
   Count the number of fields, if no separator is given, then it will
   use any number of adjacents blank characters (as given by isblank)

   @returns number of fields
*/
long
countfields(char *l,
	    char sep)
{
  //If there is an empty string return 0, otherwise there will always be
  //at least one field.
  if(!*l)
    return 0;
  long n=1;

  //Are we skiping a specific character
  if(sep){
    while(*l)
      if(*l++==sep)
	n++;
  }
  //or blanks?
  else{
    while(*l)
      if(isblank(*l++)){
	while(isblank(*l++));
	n++;
      }
  }
  return n;
}


