#include <lineread.h>

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
#define PI (3.141592653589793)
#define me (9.10938215E-31) /*Masa del electron en kg*/
#define ec (-1.602176487E-19) /*Carga del electron en Coulombs*/
#define c (299792458) /*Velocidad de la luz en m/s*/
#define h (6.62606896E-34) /*J s*/
#define k (1.3806504E-23) /*cte. de Boltzman en J/K*/

#define nT 1 /*Temperaturas*/
#define T 296 /*Temperatura en K*/
#define PREB_BREC int
#define HITRAN_RECLENGTH 160

/*Info de moleculas*/
typedef struct molec{
  int ind; /*Numero correlativo*/
  unsigned short numiso; /*Numero de isotopos*/
  double mass; /*Masa*/
  double xsect; /*Cross-section*/
  char *name; /*String con el nombre*/
  char *filename; /*String con el nombre de archivo*/

  /*struct isot[numisot];*/ /*Info de cada isotopo*/
}molec;

static const char *hitran_name="HITRAN (2008)";
static double hitran_fct=1e-4; /*cm a microns*/
static const char *hitran_fct_ac="cm";
static const char *t_u="K";
static const char *listname="hitran-molecules.txt";
static const char *listpath; /*Path de la lista*/
static const char *hitranpath; /*Hitran path*/
static FILE *list, *currmolec;

static char *dbfilename=NULL;
/*static FILE *fp=NULL;*/
static _Bool partitionread=0;
static int verbose_db=10;

molec *molecs; /*Arreglo de moleculas*/
unsigned short **corriso;

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
  double temp;
  /*const int trglength=sizeof(PREC_NREC);*/

  irec1=initial;
  irec2=final;
  do{
    *(resultp)=(irec2+irec1)/2;
    fseek(fpbs,((*resultp)-1)*(HITRAN_RECLENGTH+2),SEEK_SET);
    fseek(fpbs,3,SEEK_CUR);
    fscanf(fpbs,"%lf",&temp);
    if(lookfor>temp)
      irec1=*(resultp);
    else
      irec2=*(resultp);
  }while (irec2-irec1>1);
  *resultp=irec1+1;
}

/*que es el deftarget?*/
static _Bool
db_find(const char *name)
{
  int len = strlen(name), lent = strlen(deftarget);
  
  if (len >= lent &&
      strncmp(deftarget, name+len-lent, lent) == 0)
    return 1;

  return 0;
}



/*Open file with list of molecules to read*/
static void open_list(){
  char full[100];
  sprintf(full,"%s/%s",listpath,listname);

  if((list=fopen(full,"r+"))==NULL)
    mperror(MSGP_USER,"Could not open molecule list '%s' for reading\n",listname);
}

/*moleclist: archivo con moleculas seleccionadas para leer.
Formato:
<nº molecula> <nombre molecula> <nombre de archivo> <xsect> <nº isotopos> <masa molecula>
Comento con # las que QUIERO usar*/
static void read_list(){

  char *c1;   
  int ignorelines=4;
  int maxchar=500;
  char *line, *auxline;
  char *comp;
  int i=0;
  /*unsigned short nIso=0;*/

  c1=(char *)malloc(sizeof(char));
  comp=(char *)malloc(maxchar*sizeof(char));
  line=(char *)malloc(maxchar*sizeof(char));
  auxline=(char *)malloc(maxchar*sizeof(char));
  molecs=(molec *)malloc(NUM_MOLEC*sizeof(molec));

  open_list(); /*Abro la lista*/

  sprintf(comp,"#");
  int j;
  for(j=0;j<=ignorelines;j++){
    fgets(line,maxchar,list);
    printf("%s",line);
  }

  int *charaux1;
  char  *charaux2, *charaux3;
  double *charaux6;
  int *charaux4, *aux2; 
  float *charaux5;

  charaux1=(int *)malloc(sizeof(int));
  charaux2=(char *)malloc(maxchar*sizeof(char));
  charaux3=(char *)malloc(maxchar*sizeof(char));
  charaux4=(int *)malloc(sizeof(int));
  charaux5=(float *)malloc(sizeof(float));
  charaux6=(double *)malloc(sizeof(double));
  aux2=(int *)malloc(sizeof(int));

  corriso=calloc(NUM_MOLEC,sizeof(unsigned short *));
  
  for(int y=0;y<NUM_MOLEC;y++)
    corriso[y]=calloc(NUM_MOLEC, sizeof(unsigned short));

  int ni;
  unsigned short nit=0;

  while(!feof(list)){
    fgets(c1,2,list);
    /*sleep(1);
      printf("c=%s, comp=%s\n",c,comp);*/
    if(strncmp(c1,comp,1)==0){ /*si esta comentado, */
      fscanf(list,"%d %s %s %lf %d %f",charaux1,charaux2,charaux3,charaux6,charaux4,charaux5);
      fseek(list,sizeof(char),SEEK_CUR);
      molecs[i].ind=*charaux1;
      molecs[i].name=charaux2;
      molecs[i].filename=charaux3;
      molecs[i].numiso=*charaux4;
      molecs[i].mass=*charaux5;
      molecs[i].xsect=*charaux6;

      corriso[i]=calloc(molecs[i].numiso,sizeof(unsigned short));

      for(ni=0;ni<molecs[i].numiso;ni++){
	corriso[i][ni]=nit;
	nit++;
      }

    }else{ /*si no esta comentado*/
      fseek(list,-1,SEEK_CUR);
      fscanf(list,"%d",aux2);
      fseek(list,sizeof(char)*23,SEEK_CUR);
      /*if(*aux2>=10) printf("num=%d\n",*aux2);*/
      fgets(auxline,maxchar,list);
      molecs[i].ind=-1; /*Con esto indico que esta molecula NO VA*/
      
    }
    i++;
  } /*while!list.eof()*/
  
  fclose(list);	

} /*read_list*/

static int db_open(char *name){

  read_list();

  if((currmolec=fopen(name,"r+"))==NULL)
    mperror(MSGP_USER,"Could not open file '%s' for reading\n",name);

  /*dbfilename=name;*/

  return LR_OK;
}

static int db_close(){
  if(currmolec)
    fclose(currmolec);

  return LR_OK;
}



/*Calcula func particion por isotopo*/
static _Bool db_part(){
  int y;
  double **z;

  z=calloc(NUM_MOLEC,sizeof(double *));
  
  for(y=0;y<NUM_MOLEC;y++)
    z[y]=calloc(NUM_MOLEC, sizeof(double));

  int i;
  int ison;
  int ind;
  double g, E;

  for(i=0;i<NUM_MOLEC;i++){
    if(molecs[i].ind!=0){
      /*currmolec=fopen(molecs[i].filename,"r+");*/
      db_open(molecs[i].filename);
      while(!feof(currmolec)){
	fscanf(currmolec,"%d %*f %*f %*f %*f %*f %lf",&ind,&E);
	fseek(currmolec,98,SEEK_CUR);
	fscanf(currmolec,"%lf",&g);
	if(ind>1000){
	  if(i<10) ind/=(1E5);
	  else ind/=(1E4);
	}
	ison=ind-10*i;
	if(isnan(g*exp(-E/(k*T)))==0) /*Si no es NaN*/
	  z[i][ison]+=g*exp(-E/(k*T));
      } /*while*/
      db_close();
    } /*if*/
  }/*for*/
  
  partitionread=1;
  return 0;

}

/*Calcula gf a partir del S leido (para cada transicion)*/
static double S_to_gf(double s, double n){
  double xp, a, b;
  xp=-(h*n)/(k*T);
  a=1-exp(xp);
  b=s*me*c;
  return (b/PI*ec*ec)/a;
}

static long int db_info(struct linedb **lineinfo, double wav1, double wav2){

  if(!partitionread)
    return -1;

  int i;
  for(i=0;i<NUM_MOLEC;i++){
    if(molecs[i].ind!=-1){

      /*currmolec=fopen(molecs[i].filename,"r+");*/
      db_open(molecs[i].filename);

      struct stat st;
      if(stat(dbfilename,&st)==-1){
	mperror(MSGP_USER|MSGP_ALLOWCONT,
		"Data file '%s' cannot be accesed by stat() in "
		"function dbread_hitran().\nThis is important to obtain "
		"its size and hence the number of lines to be\n "
		"examinated\n"
		,dbfilename);
	exit(EXIT_FAILURE);
      }
      
      off_t nrec=st.st_size, zrec=ftell(currmolec);
      
      if((zrec+nrec)/HITRAN_RECLENGTH < (zrec+nrec)/(float)HITRAN_RECLENGTH){
	mperror(MSGP_USER|MSGP_ALLOWCONT,
		"Data file '%s' does not contain an integer number of "
		"%i-bytes records!.\nAre you sure it is the right '%s' "
		"file?\n"
		,dbfilename, HITRAN_RECLENGTH, hitran_name);
	exit(EXIT_FAILURE);
      }
      
      nrec/=HITRAN_RECLENGTH;
      zrec/=HITRAN_RECLENGTH;
      
      /*Revisar si esto es necesario*/
      wav1/=hitran_fct/tli_fct;
      wav2/=hitran_fct/tli_fct;
      
      MESSAGEP(verbose_db, "HITRAN (2008) Driver: Going to look for wavelength range %g - %g (%s)\n",wav1,wav2,hitran_fct_ac);

      /****************************************************************/
      off_t irec, frec;
      double dbwav1, dbwav2;
      fseek(currmolec,3,SEEK_SET);
      fscanf(currmolec,"%lf",&dbwav1);
      fseek(currmolec,-159,SEEK_END);
      fscanf(currmolec,"%lf",&dbwav2);

      if(wav1<dbwav1 || wav2>dbwav2)
	return 0;
      if(wav1>=dbwav1)
	irec=zrec;
      else{
	dbreadBSf(currmolec,zrec,nrec,wav1,&irec,HITRAN_RECLENGTH);
	MESSAGEP(verbose_db,"HITRAN 2008 driver: found initial wavelength (%f) at record %li\n",wav1,irec);

	if(wav2<=dbwav2)
	  frec=nrec-1;
	else{
	  dbreadBSf(currmolec,irec,nrec,wav2,&frec,HITRAN_RECLENGTH);
	  MESSAGEP(verbose_db,"HITRAN 2008 driver: found final wavelength (%f) at record %li\n",wav2,frec);

	  MESSAGEP(verbose_db, "HITRAN 2008 driver: Target initial and final records found in relative positions %li and %li (of range %li-%li)\n",irec-zrec,frec-zrec,zrec,nrec);

      /***************************************************************/

      struct linedb *line;
      messagep(5,
	       "\nHITRAN driver: About to initialize memory space to hold %li records.\n"
	       "HITRAN driver: I'll require %.2fMb of available memory.\n"
	       ,frec-irec+1
	       ,(frec-irec+1)*(float)sizeof(struct linedb)/1024.0/1024.0);
      if(((*lineinfo)=(struct linedb *)calloc(frec-irec+1, sizeof(struct linedb))) == NULL){
	mperror(MSGP_SYSTEM|MSGP_ALLOWCONT,
		"HITRAN 2008 driver: Cannot allocate memory to hold all the data from %s\n"
		"HITRAN 2008 driver: Required memory was %li bytes\n"
		,hitran_name, (frec-irec+1)*sizeof(struct linedb));
	exit(EXIT_FAILURE);
      }
      line = *lineinfo;
      messagep(5, "HITRAN 2008 driver: Success in memory allocation\n");


      /*Voy al primer valor que necesito*/
      fseek(currmolec,(irec-1)*(HITRAN_RECLENGTH+2),SEEK_SET);

      /*struct recordstruct{
	PREC_NREC iwl;
	short int ielow,igflog;
	}record;*/

      int ind;
      char *auxline;
      double nu, sd, E;
      auxline=(char *)malloc(HITRAN_RECLENGTH*sizeof(char));

      do{
	/*fseek(currmolec,3,SEEK_CUR);*/
	fscanf(currmolec,"%d %lf %lf %*f %*f %*f %lf",&ind,&nu,&sd,&E);
	int isorest=ind%10;
	int rest=(ind-rest)/10;
	line->isoid=corriso[rest][isorest];
	line->recpos=irec+line-*lineinfo;
	line->wl=(PREC_LNDATA)(2.0*PI/nu);
	line->elow=(PREC_LNDATA)abs(E);
	line->gf=(PREC_LNDATA)S_to_gf(sd,nu);
	
	fgets(auxline,HITRAN_RECLENGTH,currmolec);
	
	line++;
      }while(line-*lineinfo<frec-irec+1);

      return line-*lineinfo;
	}
      }
    } /*if*/
    db_close();
  } /*for*/

}

static const driver_func pdriverf_hitran4={"HITRAN (2008) driver",
					   &db_find,
					   &db_open,
					   &db_close,
					   &db_info,
					   &db_part,
};

driver_func *initdb_hitran4(){
  return (driver_func *)&pdriverf_hitran4;
}
