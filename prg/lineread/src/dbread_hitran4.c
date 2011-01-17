#include <lineread.h>
#include <math.h>

#ifndef __BYTE_ORDER
#error __BYTE_ORDER not defined, include <endian.h>
#endif
#if __BYTE_ORDER == __BIG_ENDIAN
static inline void
reversebytes(void *pointer, int nbytes)
{
  char one;
  int n,tot;
  tot = nbytes--/2;
  for(n=0 ; n<tot ; n++){
    one = *((char *)(pointer)+n);
    *((char *)(pointer)+n) = *((char *)(pointer)+nbytes-n);
    *((char *)(pointer)+nbytes-n) = one;
  }
}
#else
#define reversebytes(pointer, nbytes) do{}while(0)
#endif   /* __BYTE_ORDER */

#define NUM_MOLEC 39 /*Numero de moleculas*/
#define PI 3.141592653589793
#define me 9.10938215E−31 /*Masa del electron en kg (? revisar unidades)*/
#define e -1.602176487E−19  /*Carga del electron en Coulombs*/
#define c 299792458 /*Velocidad de la luz en m/s*/
#define h 6.62606896E-34 /*J s*/

#define nT 1 /*Temperaturas*/
#define T 296 /*Temperatura en K*/
#define PREB_BREC int
#define HITRAN_RECLENGTH 160

/*Info de moleculas*/
typedef struct molec{
  int ind; /*Numero correlativo*/
  unsigned short numisot; /*Numero de isotopos*/
  double mass; /*Masa*/
  double xsect; /*Cross-section*/
  char *name; /*String con el nombre*/
  char *filename; /*String con el nombre de archivo*/

  struct isot[numisot]; /*Info de cada isotopo*/
}molec;

/*Info de isotopos*/
typedef struct iso{
  unsigned short namelength;
  char *name;
  PREC_MASS mass;
  PREC_Z z[nT];
  PREC_CS cs[nT];
}iso;

static const char *hitran_name="HITRAN (2008)";
static double hitran_fct=1e-4; /*cm a microns*/
static const char hitran_fct_ac="cm";
static const char t_u="K";
static const char *list="molecule_list.txt";
static FILE *molec;

static char *dbfilename=NULL;
static FILE *fp=NULL;
static _Bool partitionread=0;
static int verbose_db=10;

molec *molecs; /*Arreglo de moleculas*/

#define sfread(pointer, size, nmemb, f) sfread_fcn(pointer, size, nmemb, f, __LINE__)

static void sfread_fcn(void *pointer, size_t size, size_t nmemb, FILE *f, long int line){
	   size_t read=fread(pointer,size,nmemb,f);
	   if(read!=nmemb) //Error de fread
	   	    fprintf(stderr, "%s (%li): Number of elements read was only %li of %li\n", __FILE__, line, (long)read, (long)nmemb);
}

/*Busqueda binaria: igual que en pands, solo cambio*
como va leyendo los datos (formatos son diferentes)*/
static inline void
dbreadBSf(FILE *fpbs,		/* File pointer */
	  PREC_NREC initial,	/* initial index */
	  PREC_NREC final,	/* last index */
	  double lookfor,	/* target value */
	  PREC_NREC *resultp,	/* result index */
	  int reclength)	/* Total length of record */
{
  long int irec1,irec2;
  int temp1=2*sizeof(int);
  PREC_BREC temp2;
  int aux;
  const int trglength=sizeof(PREC_BREC);

  irec1=initial;
  irec2=final;
  do{
    *(resultp)=(irec2+irec1)/2;
    fseek(fpbs,reclength*(*resultp),SEEK_SET);
    sfread(&aux,2,1,fpbs);
    sfread(&aux,2,1,fpbs);
    sfread(&temp2,trglength,1,fpbs);
    reversebytes(&temp,trglength);
    if(lookfor>temp)
      irec1=*(resultp);
    else
      irec2=*(resultp);
  }while (irec2-irec1>1);
  *resultp=irec1;
}


static int db_open(char *name){

  if((fp=fopen(name,"r+"))==NULL)
    mperror(MSGP_USER,"Could not open file '%s' for reading\n",name);

  dbfilename=name;

  return LR_OK;
}

static int db_close(){
  if(fp)
    fclose(fp);

  return LR_OK;
}

/*Open file with list of molecules to read*/
/*moleclist: archivo con moleculas seleccionadas para leer.
Formato:
<nº molecula> <nombre molecula> <nombre de archivo> <nº isotopos> <masa molecula> <cross-section molecula>
Comento con # las que QUIERO usar*/
static int read_list(char *moleclist){
  FILE *list;
  char c;   
  int ignorelines=3;
  int maxchar=100;
  char line[maxchar];
  int i=0;
  unsigned short nIso=0;

  molecs=(molec *)calloc(NUM_MOLEC,sizeof(molec));

  list=fopen(moleclist,"r+");

  for(int j=0;j<=ignorelines;j++)
    fgets(line,maxchar,list);
	   
  while(!list.eof()){
    if(fgets(c,1,line)=="#"){ /*si esta comentado,*/ 
      sfread(&(molecs[i].ind),2,1,list);
      sfread(molecs[i].name,sizeof(char)*maxchar,1,list);
      sfread(molecs[i].filename,sizeof(char)*maxchar,1,list);
      sfread(&(molecs[i].numiso),2,1,list);
      sfread(&(molecs[i].mass),4,1,list);
      sfread(&(molecs[i].xsect),4,1,list);
      nIso+=molecs[i].numiso;
    }else{ /*si no esta comentado*/
      molecs[i].ind=-1; /*Con esto indico que esta molecula NO VA*/
      
    }
    i++;
  } //while!list.eof()
  
  fclose(list);	

  int k;
  for(k=0;k<NUM_MOLEC;k++)
    if(molecs[k].ind!=-1) db_open(molecs[k].filename);

} //read_list

/*Calcula func particion*/
/*Para cada isotopo! i=numero de molecula*/
/*static int db_part(int i, char **name, PREC_TEMP *T, unsigned short *niso, PREC_MASS **mass, PREC_Z ***Z, PREC_CS ***CS){*/
static _Bool db_part(int i, iso *iso){
  long maxline=300;
  char line[maxline], *lp, *lp2;
  int ignorelines=1, k=1;
  FILE *f;

  unsigned short nIso=*niso=numiso[i];

  for(int j=0;j<=ignorelines;j++)
    fgets(line,maxline,fp);

  /*Pido memoria para Z y CS, por isotopo*/
  PREC_Z **z=*Z=(PREC_Z **)calloc(nIso,sizeof(PREC_Z *));
  PREC_CS **cs=*CS=(PREC_CS **)calloc(nIso,sizeof(PREC_CS *));
  PREC_TEMP *T=(PREC_TEMP *)malloc(sizeof(PREC_TEMP));
  *z=(PREC_Z *)calloc(nIso,sizeof(PREC_Z));
  *cs=(PREC_CS *)calloc(nIso,sizeof(PREC_CS));
  for(int k=1;k<nIso;k++){
    z[k]=*z+k;
    cs[k]=*cs+k;
  }
	
  f=fopen(molecfile[i],"r+");

}

static long int db_transition(struct linedb **lineinfo, double wav1, double wav2, char *filename){

  FILE *file=fopen(filename,"r+");

  struct stat st;
  if(stat(dbfilename,&st) == -1){
    mperror(MSGP_USER|MSGP_ALLOWCONT,
	    "Data file '%s' cannot be accesed by stat()."
	    "\nThis is important to obtain "
	    "its size and hence the number of lines to be\n "
	    "examinated\n"
	    ,dbfilename);
    exit(EXIT_FAILURE);
  }
  off_t nrec = st.st_size, zrec = ftell(fp);
  nrec/=HITRAN_RECLENGTH;
  zrec/=HITRAN_RECLENGTH;
  
  PREC_BREC db1, db2;
  int aux;
  off_t irec, frec;
  
  sfread(&aux,2,1,file);
  sfread(&aux,2,1,file);
  sfread(&db1,4,1,file);
  reversebytes(&db1,4);

  fseek(file,-HITRAN_RECLENGTH,SEEK_END);
  
  sfread(&aux,2,1,file);
  sfread(&aux,2,1,file);
  sfread(&db2,4,1,file);
  reversebytes(&db2,4);

  if(wav1>db2 || lnwav2<db1) return 0;
  if(wav1<=db1) irec=zrec;
  else{
    dbreadBSf(fp,zrec,nrec,wav1,&irec,8);
    MESSAGEP(verbose_db,"%s driver: Found final wavelength (%f) at record %li, checking twins...",hitran_name,wav2,frec);

    int temp;

    do{
      if(irec<zrec) break;
      fseek(file,HITRAN_RECLENTGH*(irec--),SEEK_SET);
      sfread(&aux,2,1,file);
      sfread(&aux,2,1,file);
      sfread(&temp,4,1,file);
    }while(temp>=wav1);

    irec++;
    MESSAGEP(verbose_db,"done (%li)\n",irec);
  }/*else*/

  if(wav2>=lndb2){
    frec=nrec-1;
  }else{
    dbreadBSf(file,irec,nrec,wav2,&frec,8);
    MESSAGEP(verbose_db,"%s driver: Found final wavelength (%f) at record %li, checking twins...",wav2,frec);
    
    int temp;
    
    do{
      fseek(file,HITRAN_RECLENGTH*(frec++),SEEK_SET);
      sfread(&aux,2,1,file);
      sfread(&aux,2,1,file);
      sfread(&temp,4,1,file);
    }while(temp<=wav2);
    
    frec--;
    MESSAGEP(verbose_db, "done (%li)\n", frec);
  } /*else*/
  
  MESSAGEP(verbose_db, "%s driver: Target initial and final records found in relative positions %li and %li (of range %li-%li)\n",hitran_name,irec-zrec,frec-zrec,zrec,nrec);
  
  struct linedb *line;
  MESSAGEP(5,"\n%s driver: about to initialize memory space to hold %li records.\nI'll require %.2fMb of available memory.\n",hitran_name,frec-irec+1,(frec-irec+1)*(float)sizeof(struct linedb)/1024.0/1024.0);

  if(((*lineinfo)=(struct linedb *)calloc(frec-irec+1, sizeof(struct linedb)))==NULL){
    mperror(MSGP_SYSTEM|MSGP_ALLOWCONT,"%s driver: Cannot allocate memory to hold all the data from %s.\n Required memory was %li bytes\n",hitran_name,hitran_name,(frec-irec+1)*sizeof(struct linedb));

    exit(EXIT_FAILURE);
  }

  line=*lineinfo;
  MESSAGEP(5,"%s driver: Success in memory allocation\n",hitran_name);

  fseek(file,irec*HITRAN_RECLENGTH,SEEK_SET);

  int cnt=0;
  int aux;
  float aux2, S;
  do{
    sfread(&aux,2,1,file);
    sfread(&aux,2,1,file);
    sfread(line->wl,4,1,file);
    reversebytes(line->wl,4);
    sfread(S,4,1,file);
    line->gf=gf(S,line->wl);
    reversebytes(line->gf,4);
    sfread(&aux2,4,1,file);
    sfread(&aux2,4,1,file);
    sfread(&aux2,4,1,file);
    sfread(line->elow,4,1,file);
    reversebytes(line->elow,4);
    line++;
  }while(line-*lineinfo<frec-irec+1);

  return line-*lineinfo;
}

/*Calcula gf a partir del S leido (para cada transicion)*/
float gf(float s, float nu){
  float xp=-(h*nu)/(k*T);
  return (s*(me*c/PI*pow(e,2)))/(1-exp(xp));
}

static const driver_func pdriverf_hitran={"HITRAN (2008) driver",
					  &db_open,
					  &db_close,
					  &read_list,
					  &db_part,
					  &db_transition};

driver_func *initdb_hitran(){
  return (driver_func *)&pdriverf_hitran;
}
