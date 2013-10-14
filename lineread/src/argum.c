/*********************************************************************
 * 
 * argum.c - Line argument processing for lineread.
 *
 * Copyright (C) 2006 Patricio Rojo
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 * 
 **********************************************************************/

/* List of functions defined here:

static void *halfzero_realloc(void *ptr, size_t size)
    Reallocs an array and sets to zero its second half

int argum(int argc, char **argv, struct hints *hint)
    Initialize options-help documentation.
    Read input arguments and set hint variable values accordingly.

void hints_free(struct hints *hint)
    Free hint variables.
*/

#include <version_lr.h>
#include <pu/procopt.h>
#include <lineread.h>


#define CFGFILE ""

/* \fcnfh
   Reallocs an array and sets to zero its second half */
static void
*halfzero_realloc(void *ptr, size_t size){
  /* size must be an even number */
  if (size&1)
    mperror(MSGP_SYSTEM,
            "halfzero_realloc() was called with an odd\n"
            "(non-even) value on its second argument, 'size'.\n");
  void *tmp = realloc(ptr,size);
  memset((char *)tmp + size/2, 0, size/2);
  return tmp;
}


/* \fcnfh
   Initialize options-help documentation.
   Read arguments and set hint variable values accordingly */
int 
argum(int argc,             /* Number of arguments in call */
      char **argv,          /* Aray of arguments in call   */
      struct hints *hint){
  /* FINDME: Explain this block: */
  enum shortargs{
    LRA_DB,
  };

  /* Initialize options-help documentation: */
  struct optdocs var_docs[]={
    /* General Arguments   */
    {NULL,       0,  HELPTITLE,    NULL,
     NULL, "GENERAL ARGUMENTS"},
    {"dry-run", 'n', no_argument,  NULL,
     NULL, "Dry run. No output is written. But otherwise "
     "the program is executed fully"},
    {"quiet",   'q', no_argument,  NULL,
     NULL, "No output other than error messages are printed"},
    {"verbose", 'v', no_argument,  NULL,
     NULL, "Increase verbose by one level for each 'v'"},
    {"help",    'h', no_argument,  NULL,
     NULL, "Prints list of possible parameters"},
    {"version", 'V', no_argument,  NULL,
     NULL, "Print version number and exit"},
    {"paramf",  'p', ADDPARAMFILE, NULL,
     "filename", "Use filename to read parameters"},

    /* Wavelength arguments */
    {NULL,    0,  HELPTITLE,         NULL,
     NULL,    "WAVELENGTH ARGUMENTS"},
    {"wavi", 'i', required_argument, "1.9",
     "value", "Value of initial wavelength to consider (in microns)"},
    {"wavf", 'f', required_argument, "2.0",
     "value", "Value of final wavelength to consider (in microns)"},
    {"wavd", 'd', required_argument, "0.5",
     "value", "Range of wavelengths to read at a time (in microns)"},

    /* Data base arguments  */
    {NULL,        0,     HELPTITLE,         NULL,
     NULL, "DATABASE ARGUMENTS"},
    {"output",   'o',    required_argument, "-",
     "filename", "Output filename.  A dash (-) indicates standard "
     "output"},
    {"database", LRA_DB, required_argument, NULL,
     "filename", "Indicates another DB to process.  Specifying "
     "--database is optional as DBs can also be specified as non-option "
     "arguments.  Note, however, that using in the same call both "
     "options (with and without explicitly specifying --database) can "
     "yield undesired consequenses as all of those with '--database' "
     "could be considered before  those without it (matching "
     "--aux could therefore be messed up)."},
    {"aux",      'a',    required_argument, NULL,
     "[n:]filename", "Auxiliary file per database. If the optional "
     "'n' is specified, then this auxiliary file correspond to the "
     "nth database, otherwise it is paired consecutively to each "
     "database.  If more auxiliary files than databases are "
     "specified the extra names are discarded."},

    {NULL, 0, 0, NULL, NULL, NULL}
  };

  struct optcfg var_cfg;
  memset(&var_cfg, 0, sizeof(var_cfg));

  var_cfg.contact = "Patricio Rojo <pato@das.uchile.cl>";
  var_cfg.files   = CFGFILE;
  var_cfg.nonopt  = "<database1> [<database2> ...]";
  var_cfg.columns = 70;
  verblevel       = 2;

  /* Start populating hint variables: */
  int radb = 1;      /* Number of auxiliary data base read      */
  int allocaux = 2;  /* Temporary number of auxiliar data bases */
  hint->dbaux = (char **)calloc(allocaux, sizeof(char *));
  int allocdb = 2;   /* Temporary number of data bases          */
  hint->ndb = 0;
  hint->db = (char **)calloc(allocdb, sizeof(char *));
  hint->datafile = NULL;
  hint->dry = 0;

  /* FINDME: What is this block of code doing?
     I don't understand what *out is             */
  int extra = 40, maxver;
  char *out = (char *)malloc(1);
  do{
    maxver = extra;
    out = (char *)realloc(out, maxver*sizeof(char));
    extra = snprintf(out, maxver, "-rc%i", version_rc);
  }while(extra>maxver);

  procopt_debug=1; /* procopt debug level */
  /* Go over the input arguments and store the values in hint accordingly: */
  while(1){
    char *endptr, *ptr;
    int adb;  /* Index for auxiliar data base */
    int rn = procopt(argc, argv, var_docs, &var_cfg, NULL); /* See pu/     */
    if (rn==-1)
      break;

    switch(rn){ /* Switch to handle each input argument */
    case 'o':
      if (hint->datafile) free(hint->datafile); /* FINDME: Explain this    */
      hint->datafile = strdup(optarg);          /* FINDME: WHAT IS optarg? */
      break;

    case 'i':
      hint->iniw = strtod(optarg, &endptr);
      if (endptr==optarg)  /* Empty string case */
        mperror(MSGP_USER, "Invalid initial wavelength.  "
                           "Run 'lineread -h' for syntax help.\n");
      break;
    case 'f':
      hint->finw = strtod(optarg, &endptr);
      if (endptr==optarg)
        mperror(MSGP_USER, "Invalid final wavelength.  "
                           "Run 'lineread -h' for syntax help.\n");
      break;
    case 'd':
      hint->delw = strtod(optarg, &endptr);
      if (endptr==optarg) 
        mperror(MSGP_USER, "Invalid wavelength range.  "
                           "Run 'lineread -h' for syntax help.\n");
      break;

    case LRA_DB:  /* The database filename */
      if(allocdb==hint->ndb) /* Maximum size reached, duplicate size */
        hint->db = (char **)realloc(hint->db, (allocdb<<=1)*sizeof(char *));
      hint->db[hint->ndb++] = strdup(optarg);
      break;

    case 'a':
      /* Set auxiliar data base number */
      adb = strtol(optarg, &ptr, 10) + 1; /* Read pointer as int         */
      if (optarg!=ptr && *ptr==':')       /* Next char is a colon        */
        optarg = ptr+1;                   /* Set pointer after the colon */
      else
        adb = radb++;
      /* Realloc until number of auxiliar data bases is >= adb: */
      while(allocaux<adb)
        hint->dbaux = (char **)halfzero_realloc(hint->dbaux,
                                        (allocaux<<=1)*sizeof(char *));
      ptr = optarg; /* FINDME: Is this line necessary? */
      hint->dbaux[adb-1] = strdup(optarg);
      break;

    case 'n':
      hint->dry = 1;
      break;

    case 'q':
      verblevel = 0;
      break;

    case 'v':
      verblevel++;
      break;

    case 'V':
      fprintf(stderr, "This is 'lineread' version %i.%i%s "
                      "(produces TLI format version %i)\n",
              version, revision, version_rc>0?out:"", TLIversion);
      lineread_free();
      procopt_free();
      exit(EXIT_SUCCESS);
      
    case '?':
      rn = optopt;
      fprintf(stderr, 
              "Unknown, unsupported, or missing parameter to option "
              "of code %i(%c) passed as argument, use '-h' to see accepted "
              "options.\n", rn, (char)rn);
      lineread_free();
      exit(EXIT_FAILURE);
      break;
    default:    /* Ask for syntax help */
      fprintf(stderr,
              "Even though option of code %i(%c) had a valid structure "
              "element, it had no switch control statement. File %s\n"
              "need to be revised.\n",
              rn, (char)rn, __FILE__);
      lineread_free();
      exit(EXIT_FAILURE);
      break;
    case 'h':
      lineread_free();
      prochelp(EXIT_SUCCESS);
      break;
    }
  }

  messagep(4, "--------------------------\n"
              "  lineread v%i.%i%s\n"
              "--------------------------\n",
           version, revision, version_rc>0?out:"");
  free(out);
  procopt_free();

  /* Store the DBs given as non-option */
  argv += optind;
  argc -= optind;  /* Number of non-option arguments */
  messagep(4,"There are %i DBs specified with --database and %i as non-option.\n",
           hint->ndb, argc);
  /* No data bases case: */
  if (hint->ndb+argc<1){
    lineread_free();
    mperror(MSGP_USER,
            "No database specified. Run 'lineread -h' for syntax help.\n");
  }

  /* Add the non-option specified data bases to hint:    */
  hint->db = (char **)realloc(hint->db, (hint->ndb+argc)*sizeof(char *));
  while(argc--)
    hint->db[hint->ndb++] = strdup(*argv++);

  /* Discard extra auxiliary data-base files, reallocate 
     auxiliary file to match real data bases:            */
  while (allocaux>hint->ndb)
    if(hint->dbaux[--allocaux]) free(hint->dbaux[allocaux]);
  hint->dbaux = (char **)realloc(hint->dbaux, hint->ndb*sizeof(char *));
  /* Set to NULL if there are less aux DBs than ndb:     */ 
  for(int i=allocaux; i<hint->ndb; i++)  
    hint->dbaux[i] = NULL;

  /* Declare data-base driver index array:               */
  hint->dbd = (int *)calloc(hint->ndb, sizeof(int));
  return 0;
}


/* Free hint variables */
void
hints_free(struct hints *hint){
  for (int i=0 ; i<hint->ndb ; i++){
    /* Free auxiliary data-base files if exist: */
    if(hint->dbaux[i])
      free(hint->dbaux[i]);
    free(hint->db[i]);
  }
  free(hint->db);
  free(hint->dbaux);
  free(hint->dbd);
  free(hint->datafile);
}

#include <version_lr.h>


