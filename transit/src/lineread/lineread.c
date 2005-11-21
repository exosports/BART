/*
 * lineread.c
 * linereaf.txc - output adequate line information for Transit.
 *                Part of Transit program.
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

/* TD: replace gabby_read in dbread_pands and lineread.h by verblevel */

#include <lineread.h>

/* Version history:
   0.5:  First light of procedure. Oct 7 2003. Patricio Rojo
   0.8:  Added 'dbid' and its effect on '.isoid' field. Added a starting
         value of 0 for 'dindex' marks in infofile. Added transitprint().
 	 100703. PMR
   0.9:  Changed 'dbid' function, now we store a correlative number in
         '.isoid' that starts at 0 in the first isotope of the first
         database and increase from there. 100903. PMR
   0.11: added option '-n' for no file output. 102603. PMR
   0.13: Verbose improved. 110903. PMR
   1.2:  Isotope information in the Infofile is now given in a
         sequential order and not per database as before. 012504. PMR
   1.3:  Version of TWII is now stored in the .inf file. 021204. PMR
   1.4:  Bug found and corrected that wrote a lot of masses instead of
         1 per isotope. 031203. PMR.
   1.5:  Change to directory structure and linking (several files
         instead of one). 030904. PMR.
   1.6:  Magic bytes are added as the first two in twii. 031504. PMR.
   1.7:  Two files merged into one, marks are no more. 031604. PMR.
   1.8:  Help fixed and default ouput changed. 032404. PMR.
   2.1:  Change from structure storage to array in TWII
         file. 032504. PMR. 
   3.1:  Change magic bits, which is 4 bytes now. Change file extension
         TWII => TLI. Big bug fix regarding reading of P&S, wasn't
         divided by freq before. 081204. PMR
*/
static int lineread_ver=3;	/* Different version implies
				   incompatible storage formats */
static int lineread_rev=1;

short gabby_dbread = 0;

/* Add sources and name for every new reader driver */
#define dbread_nfcn 2
PREC_NREC (*linefcn[dbread_nfcn])(char *,struct linedb **,float,float
				  ,char *, PREC_ZREC ***,PREC_ZREC **
				  ,PREC_ZREC **, PREC_CS ***
				  ,int *, int *, char ***)={
				    &dbread_pands, &dbread_text
				  };
char *dname[dbread_nfcn]={
  "Partridge & Schwenke (1997). Water", "TLI-ASCII"
};

//Wavelength units is always nanometers
double tli_fct=1e-7;
const char *tli_fct_name="nanometers";

/*
  synhelp: Help on syntax
*/
static void synhelp(float deltw,char *datafile, /*char *infofile,
		    float deltwmark,*/ int verblevel)
{
  fprintf(stderr,
	  "Syntax:\n\tlineread [..options..] <wavl_i> <wavl_f>\n\n"
	  " Where available options and their defaults are:\n"
	  "  -d<delt_wavl>   : range of wavelength in %s to be read\n"
	  "                    and write at once from each database (%.2f).\n"
	  "  -n              : Dummy run!, do not write anything to files\n"
	  "  -o              : output file (\"%s\"). A dash '-' indicates\n"
	  "                    that the output to be sent to the standard\n"
	  "                    output.\n"
	  "  -v and -q       : increase or decrease verbose level(%i).\n"
	  "  -h              : show this help.\n\n"
	  ,tli_fct_name,deltw,datafile/*,deltwmark,infofile*/,verblevel);
  exit(EXIT_FAILURE);
}


int main(int argc,char *argv[])
{
  /* Magic number, which in big endian would be abb3b6ff =
    {(char)0xff-'T',(char)0xff-'L',(char)0xff-'I',(char)0xff},
    or in little endian: ffb6b3ab */
  int sign=((0xff-'T')<<24)|((0xff-'L')<<16)|((0xff-'I')<<8)|(0xff);

  int i,j,rn;
  int left; //Number of databases with lines available for reorder.
  float deltw;
  double iniw,finw,parw;
  struct linedb **lineread,**crnt;
  PREC_ZREC ***Z,**T,**mass;
  PREC_CS ***cs;
  PREC_LNDATA tmin;
  int *nT,*nIso,*dbid,totaliso,adb;
  PREC_NREC *nlines,pmin;
  FILE *fpout;
  char *datafile;
  PREC_NREC dindex;
  void *tempp;
  int dummy;
  int verblevel;
  char ***isonames;


  if(dbread_nfcn<1)
    transiterror(TERR_CRITICAL,
		 "No drivers for reading database selected or found!!");


  crnt=    (struct linedb **)calloc(dbread_nfcn,sizeof(struct linedb *));
  lineread=(struct linedb **)calloc(dbread_nfcn,sizeof(struct linedb *));
  Z=       (PREC_ZREC ***)   calloc(dbread_nfcn,sizeof(PREC_ZREC **)   );
  cs=      (PREC_CS   ***)   calloc(dbread_nfcn,sizeof(PREC_CS   **)   );
  T=       (PREC_ZREC **)    calloc(dbread_nfcn,sizeof(PREC_ZREC *)    );
  mass=    (PREC_ZREC **)    calloc(dbread_nfcn,sizeof(PREC_ZREC *)    );
  nlines=  (PREC_NREC *)     calloc(dbread_nfcn,sizeof(PREC_NREC)      );
  isonames=(char ***)        calloc(dbread_nfcn,sizeof(char **)        );
  dbid=    (int *)           calloc(dbread_nfcn,sizeof(int)            );
  nIso=    (int *)           calloc(dbread_nfcn,sizeof(int)            );
  nT=      (int *)           calloc(dbread_nfcn,sizeof(int)            );

  const char *undefined_string="";
  deltw=40;
  verblevel=1;
  datafile=(char *)calloc(20,sizeof(char));
  strcpy(datafile,"-");
  fpout=NULL;
  dummy=0;

  if(argc<3)
    synhelp(deltw,datafile,verblevel);

  while(1){
    rn=getopt(argc,argv,"vhqd:o:m:n");
    if (rn==-1)
      break;

    switch(rn){
    case 'n':
      dummy=1;
      break;
    case 'h':
      synhelp(deltw,datafile,verblevel);
      break;
    case 'd':
      deltw=atof(optarg);
      break;
    case 'q':
      verblevel--;
      break;
    case 'v':
      verblevel++;
      break;
    case 'o':
      datafile=(char *)realloc(datafile,strlen(optarg)*sizeof(char));
      strcpy(datafile,optarg);
      break;
    default:
      synhelp(deltw,datafile,verblevel);
      break;
    }
  }
  iniw=atof(argv[optind++]);
  finw=atof(argv[optind]);
  gabby_dbread=verblevel;

  transitprint(1,verblevel,
	    "           LINEREAD v%i.%i. Part of Transit package\n"
	    "--------------------------------------------------------------\n",
	    lineread_ver,lineread_rev);
  transitprint(1,verblevel,
	       "Reading %i line database(s)\n\n",dbread_nfcn);

  rn=0;
  if(strcmp("-",datafile)==0){
    rn=1;
    fpout=stdout;
  }

  transitprint(2,verblevel,
	       "Extra blah blah enabled. Be prepared!\n\n"
	       "Wavelength during storage and verbose will always be in\n"
	       "%s.\n"
	       "Total wavelength range is %.2g to %.2g.\n"
	       "The databases are going to be read and written in the\n"
	       " standard TLI (Transit line information) format.\n"
	       "Processing chunks with a wavelength range of %.1f\n"
	       " %s each.\n\n"
	       ,tli_fct_name,iniw,finw,deltw,tli_fct_name);

  if(dummy)
    transitprint(1,verblevel,
		 "Dummy run: No output. However everything else actually runs\n"
		 "regularly\n");

  transitprint(1,verblevel,
	       "TLIf output file is: %s\n"
	       ,rn&1?"standard output":datafile);


  if(!dummy){
    if(fpout==NULL&&(fpout=fopen(datafile,"w"))==NULL){
      transiterror(TERR_SERIOUS,
		   "Data file '%s' cannot be opened for writing.\n"
		   ,datafile);
    }
  }

  if(finw<iniw)
    transiterror(TERR_SERIOUS,
		 "Final wavelength (%.2f) has to be greater than\n"
		 "initial wavelength (%.2f)\n",finw,iniw);

  if(!dummy){
    fwrite(&sign,sizeof(int),1,fpout);
    fwrite(&lineread_ver, sizeof(int),1,fpout);
    fwrite(&lineread_rev, sizeof(int),1,fpout);
    fwrite(&iniw,sizeof(double),1,fpout);
    fwrite(&finw,sizeof(double),1,fpout);
    rn=strlen(undefined_string);
    fwrite(&rn,sizeof(int),1,fpout);
    fwrite(undefined_string,sizeof(char),rn,fpout);
    rn=dbread_nfcn;
    fwrite(&rn,sizeof(int),1,fpout);
  }

  dindex=0;
  while(iniw<finw){
    parw=iniw+deltw;
    if(parw>finw)
      parw=finw;

    transitprint(1,verblevel,
		 "*******    Wavelength range %8.2f - %8.2f    *******\n"
		 ,iniw,parw);

    left=0;
    /* Reading of data in the range [iniw,parw] */
    for (i=0;i<dbread_nfcn;i++){
      transitprint(1,verblevel,"    Database %i (%s): \n",i+1,dname[i]);
      tempp=dindex?NULL:Z+left;	/* Only read Z if first time */
      if((nlines[left]=(linefcn[i])(NULL, lineread+left, iniw, 
				    parw, NULL, tempp, T+left,
				    mass+left, cs+left, nT+left, nIso+left,
				    isonames+left))>0){
	crnt[left]=lineread[left];

	if(left)
	  dbid[left]=dbid[left-1]+nIso[left-1];
	else
	  dbid[left]=0;
	left++;
      }
      else
	transiterror(TERR_WARNING,
		     "Database %i didn't have any line in the wavelength\n"
		     "range %f - %f, or there was an error.\n"
		     ,iniw,parw);

    }
    totaliso=dbid[left-1]+nIso[left-1];

    /* Sorting and output of tlifile is done here only while looking at
       the first range */
    if(!dindex&&!dummy){
      for(i=0;i<dbread_nfcn;i++){
	rn=strlen(dname[i]);
	fwrite(&rn,sizeof(int),1,fpout);
	fwrite(dname[i],sizeof(char),rn,fpout);
	fwrite(nT+i,sizeof(int),1,fpout);
	fwrite(nIso+i,sizeof(int),1,fpout);
	fwrite(T[i],sizeof(PREC_ZREC),nT[i],fpout);
      }
      fwrite(&totaliso,sizeof(int),1,fpout);
      for(adb=0;adb<dbread_nfcn;adb++){
	for(j=0;j<nIso[adb];j++){
	  fwrite(&adb,sizeof(int),1,fpout);
	  fwrite(mass[adb]+j,sizeof(PREC_ZREC),1,fpout);
	  transitDEBUG(22,verblevel,
		       "Just wrote mass %g at %li\n"
		       ,mass[adb][j],ftell(fpout));
	  rn=strlen(isonames[adb][j]);
	  fwrite(&rn,sizeof(int),1,fpout);
	  transitDEBUG(22,verblevel,
		       "Just wrote name's length %i at %li\n"
		       ,rn,ftell(fpout));
	  fwrite(isonames[adb][j],sizeof(char),rn,fpout);
	  fwrite(Z[adb][j],sizeof(PREC_ZREC),nT[adb],fpout);
	}
      }
    }

    transitprint(1,verblevel,"sorting... ");

    /* Merge and output of sorted data line list */
    while(left){
      tmin=crnt[0]->wl;
      pmin=0;
      for(i=1;i<left;i++){
	if(crnt[i]->wl<tmin){
	  tmin=crnt[i]->wl;
	  pmin=i;
	}
      }

      /* Instead we are leaving isoid as a correlative list and we are
	 storing its database information in the array 'dbinfo' */
      crnt[pmin]->isoid=dbid[pmin]+crnt[pmin]->isoid;

      if(!dummy){
	fwrite(&(crnt[pmin]->wl),sizeof(PREC_LNDATA),1,fpout);
	fwrite(&(crnt[pmin]->isoid),sizeof(short),1,fpout);
	fwrite(&(crnt[pmin]->elow),sizeof(PREC_LNDATA),1,fpout);
	fwrite(&(crnt[pmin]->gf),sizeof(PREC_LNDATA),1,fpout);
      }
      dindex++;

      if(++crnt[pmin]-lineread[pmin]>=nlines[pmin]){
	for(i=pmin;i<left-1;i++){
	  dbid[i]=dbid[i+1];
	  crnt[i]=crnt[i+1];
	}
	left--;
      }
    }
    //    Pprintf(2,"DD:read %f - %f\n",iniw,parw);
    iniw+=deltw;
    for (i=0;i<dbread_nfcn;i++){
      free(lineread[i]);
      free(Z[i]);
      free(T[i]);
      free(mass[i]);
    }
    transitprint(1,verblevel,"done\n");
  }
  if(!dummy)
    fclose(fpout);

  transitprint(1,verblevel,"\n");

  return EXIT_SUCCESS;
}


