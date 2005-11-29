/*
 * dbread_pands.c - Driver to read Partridge & Schwenke for Transit.
 *              Part of Transit program.
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
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include <transit.h>
#include <lineread.h>
#include <math.h>

/* TD: Is not safe to use structre to read data, because it may
   contain padding */
/* TD: BS can be fooled by lines with identical line centers */

#define PANDS_RECLENGTH 8
#define PANDS_NCODIDX 32786
#define PREC_BREC int		/* Type for the Big record  */
#define NUM_ISOT 4
#define DRIVERNAME "pands"

static char *pands_name="Partridge & Schwenke(1997)";
static char *isotope[NUM_ISOT]={"1H1H16O","1H1H17O","1H1H18O","1H2H16O"};
/* Wavelength is stored in nanometers */
static double pands_fct=1e-7;
static const char *pands_fct_ac="nm";
short gabby_dbread;

static int isoname(char ***isotope, int niso);
static int read_zpands(char *filename, PREC_ZREC ***Z, PREC_ZREC **T,
		       int *nT, int nIso);



/*
  reversebytes: Invert the byte order. It can be used for any variable
  type, just specify the number of bytes of its class through the use of
  sizeof().

  @requires pointer to be a pointer to any type and byte to be an
  unsigned int.
*/
#if __BYTE_ORDER == __BIG_ENDIAN
#define reversebytes(pointer, nbytes)   do{                           \
             char one;                                                \
             int n,tot;                                               \
             tot=nbytes--/2;                                          \
             for(n=0;n<tot;n++){                                      \
               one=*((char *)(pointer)+n);                            \
               *((char *)(pointer)+n)=*((char *)(pointer)+nbytes-n);  \
               *((char *)(pointer)+nbytes-n)=one;                     \
             }                                                        \
             while(0)
#else
#define reversebytes(pointer, nbytes) do{}while(0)
#endif   /* __BYTE_ORDER */

/*
  dbreadBSf: Perform a binary search in file pointed by 'fp'(FILE *)
  between 'initial'(PREC_NREC) and 'final'(PREC_NREC) looking for
  'lookfor'(PREC_NREC) at the first item of a record, result is stored
  in 'resultp'(PREC_NREC *). Records are of length 'reclength'(int) and
  the first item of each of them is of type PREC_BREC.
*/
static inline void
dbreadBSf(FILE *fp,		/* File pointer */
	  PREC_NREC initial,	/* initial index */
	  PREC_NREC final,	/* last index */
	  double lookfor,	/* target value */
	  PREC_NREC *resultp,	/* result index */
	  int reclength)	/* Total length of record */
{
  int irec1,irec2;
  PREC_BREC temp;
  const int trglength=sizeof(PREC_BREC);

  irec1=initial;
  irec2=final;
  do{
    *(resultp)=(irec2+irec1)/2;
    fseek(fp,reclength*(*resultp),SEEK_SET);
    fread(&temp,trglength,1,fp);
    reversebytes(&temp,trglength);
    if(lookfor>temp)
      irec1=*(resultp);
    else
      irec2=*(resultp);
  }while (irec2-irec1>1);
  *resultp=irec1;
}


PREC_NREC
dbread_pands(char *filename,
	     struct linedb **lines, //2 pointers in order to be
				//able to allocate memory
	     float wlbeg,           //wavelengths in tli_fct
	     float wlend,           //units
	     /* Partition function data file */
	     char *Zfilename,
	     /* For the following 3 parameter, the memory is
		allocated in the dbread_* functions, and the
		size is returned in the last parameters. */
	     PREC_ZREC ***Z,        //Partition function:[iso][T]
	     PREC_ZREC **T,         //temps for Z
	     PREC_ZREC **isomass,   //Isotopes' mass in AMU
	     PREC_CS ***isocs,      //Isotope's cross-section: [iso][T]
	     int *nT,               //number of temperature
				//points 
	     int *nIso,             //number of isotopes
	     char ***isonames)      //Isotope's name
{
  PREC_NREC nrec;
  char *deffname  = "./oth/pands/h2ofast.bin";
  char *defzfname = "./oth/pands/h2opartfn.dat";

  struct linedb *line;
  struct recordstruct{
    PREC_BREC iwl;                    //Wavelength index
    short int ielow,igflog;           //Lower energy and log(gf) indices.
  }record;
  struct stat fs;
  double tablog[PANDS_NCODIDX+1];     //Array used to convert stored int
				      //to float
  double lnwl;                        //log of beginning wavelength to
				      //be used in BS
  PREC_NREC lnwl1, lnwl2, irec, frec; //upper and lower values of BS
  double ratiolog;
  PREC_NREC i;
  FILE *fp;

  wlbeg*=pands_fct/tli_fct;
  wlend*=pands_fct/tli_fct;

  *nIso=NUM_ISOT;
  if(!filename)
    filename=deffname;
  if(!Zfilename)
    Zfilename=defzfname;

  /*
    Following is to transform stored integer of gflog and ielow into their
    real floating point values
  */
  for(i=1;i<=PANDS_NCODIDX;i++)
    tablog[i]=pow(10,(i-16384)*0.001);

  ratiolog=log((double)1.0+(double)1.0/(2e+6));

  if(stat(filename,&fs)==-1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Data file '%s' cannot be accesed by stat() in "
		 "function dbread_pands().\nThis is important to obtain "
		 "its size and hence the number of lines to be\n "
		 "examinated\n"
		 ,filename);
    return -2;
  }
  nrec=fs.st_size;
  if(nrec/PANDS_RECLENGTH<nrec/(float)PANDS_RECLENGTH){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Data file '%s' does not contain an integer number of "
		 "%i-bytes records!.\nAre you sure it is the right %s "
		 "file?\n"
		 ,filename,PANDS_RECLENGTH,pands_name);
    return -3;
  }
  nrec/=PANDS_RECLENGTH;

  if((fp=fopen(filename,"r"))==NULL){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Data file '%s' cannot be opened\nas stream in function "
		 "dbread_pands(), stopping.\n"
		 ,filename);
    return -1;
  }
  fread(&lnwl1,4,1,fp);
  reversebytes(&lnwl1,4);

  fseek(fp,-PANDS_RECLENGTH,SEEK_END);
  fread(&lnwl2,4,1,fp);
  reversebytes(&lnwl2,4);

  lnwl=log(wlbeg)/ratiolog;
  if(lnwl>lnwl2){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "Last wavelength in the file (%g %s) is shorter "
		 "than\nrequested initial wavelength (%g %s)\n"
		 ,exp(lnwl2*ratiolog), pands_fct_ac, wlbeg
		 ,pands_fct_ac);
    return -4;
  }
  if(lnwl<lnwl1){
    wlbeg=exp(ratiolog*lnwl1);
    irec=0;
  }
  else
    dbreadBSf(fp,0,nrec,lnwl,&irec,8);
    //    pandsBS(0,nrec,lnwl,&irec);
  if(gabby_dbread>1)
    fprintf(stderr,
	    "Located beginning wavelength %g at position %li\n",wlbeg,irec);

  lnwl=log(wlend)/ratiolog;
  if(lnwl<lnwl1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
		 "First wavelength in the file (%g %s) is longer"
		 "than\nrequested ending wavelength (%g %s)\n"
		 ,exp(lnwl2*ratiolog), pands_fct_ac, wlbeg
		 ,pands_fct_ac);
    return -5;
  }
  if(lnwl>lnwl2){
    wlend=exp(ratiolog*lnwl1);
    frec=0;
  }
  else{
    dbreadBSf(fp,irec,nrec,lnwl,&frec,8);
    //    pandsBS(irec,nrec,lnwl,&frec);
    frec++;
  }
  if(gabby_dbread>1){
    fprintf(stderr,
	    "Located ending wavelength %g at position %li\n",wlend,frec);
    fprintf(stderr,
	    "About to initialize memory space to hold %li records.\n"
	    " I'll require %.2fMb of available memory.\n"
	    ,frec-irec
	    ,(frec-irec)*(float)sizeof(struct linedb)/1024.0/1024.0);
  }

  if(((*lines)=calloc(frec-irec,sizeof(struct linedb)))==NULL){
    transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
		 "Cannot allocate memory to hold all the data from %s\n"
		 "Required memory was %li bytes\n"
		 ,pands_name,(frec-irec)*sizeof(struct linedb));
    return -6;
  }
  if(gabby_dbread>1)
    fprintf(stderr,"Success in memory allocation\n");

  /* Main reading loop */
  if(gabby_dbread>0)
    fprintf(stderr,"reading... ");

  fseek(fp,irec*PANDS_RECLENGTH,SEEK_SET);

  i=0;
  do{
    line=(*lines)+i;

    /*TD: Is not safe to use structre to read data, because it may
      contain padding */

    fread(&record,1,8,fp);
    reversebytes(&(record.iwl),4);
    reversebytes(&(record.ielow),2);
    reversebytes(&(record.igflog),2);

#if 0
    int b;
    fprintf(stderr,"\n%08li  ",irec+i);
    for(b=0;b<8;b++)
      fprintf(stderr,"%02hhx ",(char)*(((char *)&record)+b));
#endif

    line->wl=exp(record.iwl*ratiolog)*pands_fct/tli_fct;

    //Isotopes (1h1h16o, 1h1h17o, 1h1h18o, 1h2h16o)
    if(record.ielow>0)
      line->isoid=record.igflog>0?0:1; //isotopes indices are Kurucz's - 1
    else
      line->isoid=record.igflog>0?2:3;
    line->recpos=i+irec;
    line->elow=(PREC_LNDATA)abs(record.ielow);
    line->gf=tablog[abs(record.igflog)];

    i++;
  }while(line->wl<wlend && i+irec<frec);
  if(gabby_dbread>0)
    fprintf(stderr,"done\n");

  fclose(fp);

  if(isonames!=NULL)
    isoname(isonames,*nIso);

  *isomass=(PREC_ZREC *)calloc(*nIso,sizeof(PREC_ZREC));
  (*isomass)[0]=18.01056468;
  (*isomass)[1]=19.01478156;
  (*isomass)[2]=20.01481046;
  (*isomass)[3]=19.01684143;

  if(Z!=NULL)
    if((nrec=read_zpands(Zfilename,Z,T,nT,*nIso))!=1)
      transiterror(TERR_SERIOUS,
		   "Function read_zpands() return error code '%i'\n",
		   nrec);

  return i;
}


#define MAX_LINE 200
/*\fcnfh
  read_zpands: Read from file 'filename'(char *), the partition function
  information into 'Z'(PREC_ZREC ***), it is a two dimensional array
  depending on isotope (1st dimension) and temperature(2nd
  dimension). Values of temperature are stored in 'T'(PREC_ZREC **), of
  size 'nT'(int *). 'nIso'(int) indicates the number of isotopes in
  consideration. It is assumed that the isotopes columns have the same
  order as the indices being used.
*/
static int read_zpands(char *filename, /* Doh! */
		       PREC_ZREC ***Z, /* Partition function */
		       PREC_ZREC **T,  /* Temperature points */
		       int *nT,        /* Number of temp. points */
		       int nIso)       /* Number of isotopes */
{
  FILE *fp;
  int i,cnt;
  char line[MAX_LINE], *sp,*sp2;
  int ignorelines;

  ignorelines=5;		/* Header */
  (*nT)=8; 			/* Initial value for allocation of
				   temperature info */

  if((fp=fopen(filename,"r"))==NULL){
    transiterror(TERR_SERIOUS,
		 "Data file '%s' cannot be opened as stream in function "
		 "read_zpands(), stopping.\n"
		 ,filename);
  }

  for(i=0;i<ignorelines;i++){
    fgets(line,MAX_LINE,fp);
  }

  *Z=(PREC_ZREC **)calloc(nIso,sizeof(PREC_ZREC *));
  for(i=0;i<nIso;i++)
    (*Z)[i]=(PREC_ZREC *)calloc((*nT),sizeof(PREC_ZREC));
  *T=(PREC_ZREC *)calloc((*nT),sizeof(PREC_ZREC));

  cnt=0;

  while(fgets(line,MAX_LINE,fp)!=NULL){
    sp=line;
    (*T)[cnt]=strtod(sp,&sp2);

    for(i=0;i<nIso;i++){
      if(sp2==sp){
	transiterror(TERR_SERIOUS,
		     "In function read_zpands(): line %i of file\n '%s'"
		     " has %i columns instead of %i.\n",
		     cnt+ignorelines,filename, i+1,nIso+1);
      }
      sp=sp2;
      (*Z)[i][cnt]=strtod(sp,&sp2);
    }

    if(++cnt==(*nT)){
      (*nT)<<=1;
      (*T)=(PREC_ZREC *)realloc((*T),(*nT)*sizeof(PREC_ZREC));
      for(i=0;i<nIso;i++)
	(*Z)[i]=(PREC_ZREC *)realloc((*Z)[i],(*nT)*sizeof(PREC_ZREC));
    }
  }
  (*nT)=cnt;
  (*T)=(PREC_ZREC *)realloc((*T),(*nT)*sizeof(PREC_ZREC));
  for(i=0;i<nIso;i++)
    (*Z)[i]=(PREC_ZREC *)realloc((*Z)[i],(*nT)*sizeof(PREC_ZREC));


  return 1;

}
#undef MAX_LINE

/*
  isoname: returns the name of the isotopes in the newly allocated array
  'isoname'(char ***)

  @todo    error check for calloc calls
  @returns 1 on success
*/
static int isoname(char ***isonames, int niso)
{
  int i;

  *isonames=(char **)calloc(niso,sizeof(char *));
  for(i=0;i<niso;i++){
    (*isonames)[i]=(char *)calloc(strlen(isotope[i])+1,sizeof(char));
    strcpy((*isonames)[i],isotope[i]);
  }

  return 1;
  
}

#ifdef TEST_RUN

/* \fcnfh
   
 */
static int
test_fill_defaults(char *file,
		   char *zfile,
		   char *dbread_wl_unit,
		   int *numisot,
		   char ***isonames,
		   PREC_NREC (*dbread_main)(char *,  struct linedb **,
					    float, float, char *, 
					    PREC_ZREC ***, PREC_ZREC **, 
					    PREC_ZREC **, PREC_CS ***,
					    int *, int *, char ***))
{
  strcpy(file,"../oth/pands/h2ofast.bin");
  strcpy(zfile,"../oth/pands/h2opartfn.dat");
  strcpy(dbread_wl_unit, pands_fct_ac);
  *numisot=NUM_ISOT;
  *isonames = (char **)calloc(NUM_ISOT, sizeof(char *));
  **isonames = (char *)calloc(sizeof(isotope)+NUM_ISOT, sizeof(char));

  fprintf(stderr, "Allocating %i bytes for isotope. Correct and check.\n"
	  , sizeof(isotope)+1);

  int i;
  for (i=0 ; i<NUM_ISOT ; i++){
    int len=strlen(strcpy(*isonames[i], isotope[i]))+1;
    if (i!=NUM_ISOT-1)
      *isonames[i+1] = *isonames[i] + len;
  }

  dbread_main = &dbread_pands;
  return 0;
}

#endif
