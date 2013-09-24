/*
 * procopt.c - Practical command line option parsing
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


/*
 * To use procopt, it is needed to be familiar with the use of GNU's
 * getopt() (see section 3 of manpages).  The usage is identical as
 * both return the same, and thus procopt() should be called inside a
 * loop followed by a switch() statement.  See doc/procopt_sample.c
 * for a simple call.  The procopt.c routines extend the getopt
 * routines by allowing the specification of a filename from which the
 * parameters will be read.
 * 
 * The main function is 
 *
 *   int procopt(int argc, char **argv, struct optdocs *opts,
 *               struct optcfg *cfg, int *longidxp)
 *
 * 'cfg' receives any configuration parameter that wants to be
 * customized (it is not a mandatory argument and can be given a NULL
 * instead).  'longidxp' is analog to getopt()'s 'longindex' argument.
 * 'opts' is a read-only array of structures that will contain all the
 * information regarding each of the parameters that want to be
 * processed:
 *
 * struct optdocs {
 *         char *name; ***long name of the parameter***
 *         int val;    ***short option (if ASCII-alphanumeric) and/or
 *                        identifier for the switch statement ****
 *         int has_arg;***analog to getopt's has_arg field, except for
 *                        2 extra values it accept: HELPTITLE and
 *                        PARAMFILE.  The former indicates a that
 *                        the 'name' field indicates the name of the
 *                        subsection as it will appear in the help
 *                        output of prochelp().  PARAMFILE indicates
 *                        that that option takes the name of the
 *                        parameter file that will be processed for
 *                        other options.  Options with these 2 new
 *                        options will be processed internally and
 *                        transparently to users of procopt()***
 *         char *def;  ***Default value. If NULL, no default is set***
 *         char *prm;  ***Name of the option's parameter.  This is
 *                        shown among <> in the prochelp() output)***
 *         char *doc;  ***Description of the parameter***
 * };
 *
 *
 * 
 */



#include <pu/procopt.h>
#include <errno.h>

//variables that are going to have allocated a space
static char *prgname=NULL;
static struct option *getopts;
static char *shortopts;
static FILE **fp=NULL;
static char *paramfiles=NULL;
static char *line;

static struct optdocs *_opts;
static struct optcfg *_cfg;

unsigned short procopt_debug=1;
static _Bool freed=0;
static int givenparamf=-1;

//The following relate to settle the parameter list and processing defaults
static struct optdocs *_opts_def;
static int getalloc=8,shortalloc=8;
static struct option *co;
static char *cs;
static int getn,shortn;
static _Bool process_defaults=0;

static int getoptfrom(char *line, struct option *getopts, int *longindex);
static void cfg_free();

/* \fcnfh
   Fill short and long options and set defaults

   @returns value to be processed
            -1 if no default value or finished
*/
static int
fillanddef()
{
  if(getn==getalloc){
    getopts=(struct option *)realloc(getopts,
				     (getalloc<<=1)*sizeof(struct option));
    co=getopts+getn;
  }
  if(shortn>shortalloc){
    shortopts=(char *)realloc(shortopts,(shortalloc<<1)*sizeof(char));
    memset(shortopts+shortalloc, 0, shortalloc);
    shortalloc <<= 1;
  }
  if(_opts_def->name==NULL&&_opts_def->val==0){
    if(_opts_def++->has_arg==HELPTITLE)
      return fillanddef();
    if(procopt_debug>2){
      fprintf(stderr,
	      "----------------------------------------------\n"
	      "procopt_debug message:: Accepted short string '%s'\n"
	      "(set procopt_debug to less than 3 if you don't want\n"
	      "to see this again).\n"
	      "----------------------------------------------\n"
	      ,shortopts);
    }
    co->name=NULL;
    co->has_arg=0;
    co->flag=NULL;
    co->val=0;
    getn++;
    getopts=(struct option *)realloc(getopts,getn*sizeof(struct option));
    shortopts=(char *)realloc(shortopts,(shortn-1)*sizeof(char));
    process_defaults=0;
    return -1;
  }
  
  //save short options
  if(_opts_def->val>0x20 && _opts_def->val<0x80){
    cs=shortopts;
    while(*cs){
      if(*cs++==(char)_opts_def->val&&procopt_debug){
	fprintf(stderr,
		"procopt_debug error:: The short option '%c' appears more\n"
		" than once in the given 'struct optdocs'.\n"
		" Currently with %s\n"
		,_opts_def->val,_opts_def->name);
	exit(EXIT_FAILURE);
      }
    }
    *cs++=(char)_opts_def->val;
    shortn++;
    if(_opts_def->has_arg==required_argument ||
       _opts_def->has_arg==ADDPARAMFILE){
      *cs=':';
      shortn++;
    }
  }

  //save long option
  co->name=_opts_def->name;
  if(_opts_def->has_arg==ADDPARAMFILE){
    if(givenparamf!=-1&&procopt_debug){
      fprintf(stderr,
	      "procopt_debug error:: More than one option with the ADDPARAMFLAG\n"
	      "mode, only the last will be taken into account\n"
	      );
      exit(EXIT_FAILURE);
    }
    co->has_arg=required_argument;
    givenparamf=_opts_def->val;
  }
  else
    co->has_arg=_opts_def->has_arg;
  co->flag=NULL;
  co->val=_opts_def->val;

  //advance to next value
  getn++;
  _opts_def++;
  co++;

  //Run default if there is any.
  if(_opts_def[-1].def!=NULL){
    optarg=_opts_def[-1].def;
    return _opts_def[-1].val;
  }

  return -1;
}



/* \fcnfh
   Process comand line parameters. It process paramfile parameters and
   then the command line parameters.
*/
int
procopt(int argc,		/* Number of command line arguments */
	char **argv,		/* Command line parameters */
	struct optdocs *opts,	/* structure with information about
				   parameter options and syntax */
	struct optcfg *cfg,	/* various configuration parameters for
				   procopt */
	int *longidxp)		/* return longindexp. Seems to be not
				   working in libraries!, though */
{
  char *cfgfiles=NULL;
  int ret;

  if(freed){
    fprintf(stderr,
	    "procopt:: getprocopt() was called after a call to\n"
	    "getpropt_free(), which should be the very last {procopt}\n"
	    "function called\n");
    exit(EXIT_FAILURE);
  }

  //process next element of the array, returning to set default if there
  //is one.
  if(process_defaults){
    while((ret=fillanddef())==-1 && process_defaults);
    if(ret!=-1)
      return ret;
  }

  //If first time then initialize
  if(!prgname){
    //save structures to later call by prochelp
    _opts=opts;

#if 0  //Following commented as default is dealt in prochelp
    if (!cfg){
      _cfg = calloc(1, sizeof(struct optcfg));
      _cfg.freecfg=1;
    }
    else
#endif
    _cfg=cfg;

    if((cs=rindex(*argv,'/')))
      cs++;
    else
      cs=*argv;
    prgname=(char *)calloc(strlen(cs)+1,sizeof(char));
    strcpy(prgname,cs);

    //If there is no 'opts' structure then allocate an empty string as
    //'shortopts' and an empty structure for 'getopts'
    if(!opts){
      getopts=(struct option *)calloc(1,sizeof(struct option));
      shortopts=(char *)calloc(1,sizeof(char));
    }

    //otherwise look for the right values
    else{
      co=getopts=(struct option *)calloc(getalloc,sizeof(struct option));
      cs=shortopts=(char *)calloc(shortalloc,sizeof(char));
      //getn starts at zero, but 'shortn' starts at 2 because it has to
      //reserve a value for the final `\\0' and a possible colon `:'
      getn=0;
      shortn=2;
      if(cfg && (cfg->argmode=='+' || cfg->argmode=='-')){
	cs[0]=cfg->argmode;
	shortn++;
      }

      //Analize first element and set flag to process the rest
      _opts_def=opts;
      process_defaults=1;
      while((ret=fillanddef())==-1 && process_defaults);
      if(ret!=-1)
	return ret;
    }
  }

  //Set configuration if exists
  if(cfg)
    cfgfiles=cfg->files;

  //Return the corresponding option.
  return getopt_long_files(argc,argv,shortopts,getopts,longidxp,cfgfiles);
}


/* \fcnfh
 Works as getopt_long except by the extra parameter 'paramfile' which
 can process long versions of options from a file in the format
 \begin{verb} long_option [=] value\end{verb}.
 Reordering of arguments is also done in the same way as getopt_long
 does.

 @return shortoption value as it would have been returned by
                     getopt_long.
*/
int 
getopt_long_files(int argc,	/* number of arguments */
		  char **argv,	/* agruments list */
		  char *shortopts, /* short options accepted, format is
				      the same as getopt(). */
		  struct option *getopts, /* long options accepted,
					     format is the same as
					     getopt_long(). Although
					     optional arguments are not
					     supported */
		  int *longidxp, /* returns position in the array
				    getopts[] for the last selected
				    option */
		  char *paramfilelist) /* List of files from which to
					  process longoptions. Names are
					  separated by commas. Value of
					  this variable in the first call
					  to this function, is the only
					  one that matters. */
{
  static int fpn=-1, fpa=4;
  static _Bool needs_open=0;
  int ret;
  char *fn;

  if(freed){
    fprintf(stderr,
	    "procopt:: getopt_long_files() was called after a call to\n"
	    "getpropt_free(), which should be the very last {procopt}\n"
	    "function called\n");
    exit(EXIT_FAILURE);
  }
  
  //Store default files if it is first time that is called
  if(!paramfiles){
    if(paramfilelist){
      paramfiles=calloc(strlen(paramfilelist)+1,sizeof(char));
      strcpy(paramfiles,paramfilelist);
    }
    else
      paramfiles=calloc(1,sizeof(char));
    fp=(FILE **)calloc(fpa,sizeof(FILE *));
  }

  //if we have to work only in command line arguments
  if(fpn==-1 && !*paramfiles && !needs_open)
    ret= getopt_long(argc,argv,shortopts,getopts,longidxp);
  
  //otherwise if we are processing a file.
  else{
    //if we are opening a non default file check that we have allocated
    //enough file pointers. In this case a nonopenable file cause
    //error.
    if(needs_open){
      needs_open=0;
      if(fpn==fpa-1) fp=(FILE **)calloc(fpa<<=1,sizeof(FILE *));
      fp[++fpn]=fopen(optarg,"r");
      if(fp[fpn]==NULL){
	fprintf(stderr,
		"Unable to succesfully open parameter file '%s'\n"
		,optarg);
	exit(EXIT_FAILURE);
      }
    }
    //if we need to open a file from the defaults
    else if(fpn==-1){
      fn=paramfiles;
      while(*fn&&*fn!=',') fn++;
      char tmp=*fn;
      *fn='\0';
      fp[++fpn]=fopen(paramfiles,"r");
      if(tmp)
	strcpy(paramfiles,++fn);
      else
	paramfiles[0]='\0';
    }

    //if file was opened then process next line
    if(fp[fpn]){
      //skip empty lines and commented lines. If this is not the first
      //line read, free the array
      while(1){
	if(line)
	  free(line);
	//if we reach the end of file, close it and proceed with the
	//following field.
	if((line=fgets_alloc(fp[fpn],NULL))==NULL){
	  fclose(fp[fpn--]);
	  return getopt_long_files(argc,argv,shortopts,getopts,
				   longidxp,NULL);
	}
	if(*line&&*line!='#')
	  break;
      }
      //If line was read succesfully, find if option exist and set
      //optarg 
      ret=getoptfrom(line,getopts,longidxp);
    }
    //otherwise, if file was not opened then process next file
    else{
      fpn--;
      errno=0;
      return getopt_long_files(argc,argv,shortopts,getopts,
			       longidxp,NULL);
    }
  }

  //If the option was to process a new parameter file then do that. That
  //option is hence never returned out of this function.
  if(givenparamf>=0&&ret==givenparamf){
    needs_open=1;

    return getopt_long_files(argc,argv,shortopts,getopts,
			     longidxp,NULL);
  }

  return ret;
}


/* \fcnfh
   Gives command line summary
*/
void
prochelp(int status)		/* either 'EXIT_FAILURE' or
				   'EXIT_SUCCESS' */
{
  char *options   = "[options]";
  char *intro     = 
    " Where [options] are listed below\n";
  char *intro2 =
    "  (note that whenever there is a mandatory argument, it is\n"
    "mandatory for both short and long options)...\n";
  char *contintro = "----------------------------------\n"
    "Contact Information: %s\n";
  char *hprogname = xstrdup(prgname);
  char *noinfo = "(No info available)";
  int currind,tmp;
  char *doc;
  int cols;
  FILE *ostream;
  short helpmode =
    (_cfg && _cfg->helpmode)? _cfg->helpmode: 0;
  int procopt_columns =
    (_cfg && _cfg->columns)? _cfg->columns:__PADTO_COLUMNS;
  char *defword = 
    (_cfg && _cfg->defword)? (char *)_cfg->defword:" (default: ";
  char *postdefword = 
    (_cfg && _cfg->postdefword)? (char *)_cfg->postdefword:")";
  char *enddocchar = 
    (_cfg && _cfg->enddocchar)? (char *)_cfg->enddocchar:".";
  char *pretitle = 
    (_cfg && _cfg->pretitle)? (char *)_cfg->pretitle:"\n";
  char *posttitle = 
    (_cfg && _cfg->posttitle)? (char *)_cfg->posttitle:"\n";
  char *postoption =
    (_cfg && _cfg->postoption)? (char *)_cfg->postoption:"\n";
  int indentdoc   = (helpmode == 0)? 10: 21;
  if(_cfg){
    if(_cfg->endpadchar)
      printpad_endpadchar = _cfg->endpadchar;
    if(_cfg->indentdoc)
      indentdoc = _cfg->indentdoc;
    if(_cfg->intro)
      intro = xstrdup((char *)_cfg->intro);
    if(_cfg->intro2)
      intro2 = xstrdup((char *)_cfg->intro2);
    if(_cfg->options)
      options = xstrdup((char *)_cfg->options);
    if(_cfg->noinfo)
      noinfo = xstrdup((char *)_cfg->noinfo);
    if(_cfg->prg)
      hprogname = xstrdup((char *)_cfg->prg);
  }
  if(_cfg && _cfg->usestderr)
    ostream = stderr;
  else
    ostream = stdout;

  if(freed){
    fprintf(stderr,
	    "procopt:: prochelp() was called after a call to\n"
	    "getpropt_free(), which should be the very last {procopt}\n"
	    "function called\n");
    exit(EXIT_FAILURE);
  }
  if(!prgname && procopt_debug > 1){
    fprintf(stderr,
	    "----------------------------------------------\n"
	    "procopt_debug warning:: Prochelp was called before getprocopt and\n"
	    "hence there are no available help\n"
	    "(set procopt_debug to less than 2 if you don't want to see this again).\n"
	    "----------------------------------------------\n"
	    );
  }

  doc = getenv("COLUMNS");
  if(!doc || !(cols=atoi(doc)))
    cols = procopt_columns;

  char out[cols+1];

  if(!prgname){
    fprintf(stderr,
	    "procopt error:: prochelp() was called with status %i\n"
	    "before a call to getprocopt()\n\n"
	    ,status);
    procopt_free();
    exit(EXIT_FAILURE);
  }

  fprintf(ostream,"Usage:\n\t%s %s", prgname, options);

  if(_cfg&&_cfg->nonopt)
    fprintf(ostream, " %s", _cfg->nonopt);


  //process options now
  fprintf(ostream,"\n\n");
  if(_opts){

    //print intro
    fprintf(ostream, " %s", intro);
    fprintf(ostream, "%s\n", intro2);

    //process each option
    while(1){
      currind = 0;
      if(_opts->has_arg == HELPTITLE)
	fprintf(ostream, "%s%s%s", pretitle, _opts->doc, posttitle);
      else if(_opts->name == NULL && _opts->val == 0)
	break;
      else{
	//if short form is of a printable character
	if(_opts->val > 0x20 && _opts->val < 0x80){
	  fprintf(ostream, " -%c", _opts->val);
	  if(helpmode == 1){
	    if(_opts->has_arg == required_argument ||
	       _opts->has_arg == ADDPARAMFILE)
	      fprintf(ostream, " <%s>%n", _opts->prm, &currind);
	    else if(_opts->has_arg != no_argument){
	      fprintf(stderr,
		      "\n\nprocopt error:: a non-supported value (%i) was given\n"
		      "in .has_arg field of parameter '%c'\n"
		      ,_opts->has_arg, _opts->val);
	      exit(EXIT_FAILURE);
	    }
	  }
	  currind += 2;
	}
	//stop if it is not printable and there is no long alternative
	else if(!_opts->name){
	  fprintf(stderr,
		  "\n\nprocopt error:: in prochelp() only a non displayable\n"
		  "value was given (val: %i). Parameter name is '%s' and\n"
		  "document help is:\n%s\n"
		  ,_opts->val, _opts->prm, _opts->doc);
	  exit(EXIT_FAILURE);
	}
	//print long version if there is one
	if(_opts->name){
	  if (currind) {
	    if (helpmode == 0) fprintf(ostream, ",");
	    else if (helpmode == 1) fprintf(ostream, "\n");
	  }
	  fprintf(ostream, " --%s%n", _opts->name, &currind);
	  tmp = 0;
	  if(helpmode == 1){
	    if(_opts->has_arg == required_argument ||
	       _opts->has_arg == ADDPARAMFILE)
	      fprintf(ostream, " <%s>%n", _opts->prm, &tmp);
	    else if(_opts->has_arg != no_argument){
	      fprintf(stderr,
		      "\n\nprocopt error:: a non-supported value (%i) was "
		      "given\nin .has_arg field of parameter '%s'\n"
		      ,_opts->has_arg, _opts->name);
	      exit(EXIT_FAILURE);
	    }
	  }
	  currind += tmp;
	}
	//print argument in helpmode 0
	if (helpmode == 0){
	  if(_opts->has_arg == required_argument ||
	     _opts->has_arg == ADDPARAMFILE)
	    fprintf(ostream, "%c<%s>", currind>2? '=':' ', _opts->prm);
	  fprintf(ostream, "\n");
	  tmp = 0;
	  currind = 0;
	}
	//'tmp' is now spaces to indentation, and 'currind' is either the
	//'indentdoc' or wherever the pre-indent ends.
	tmp = indentdoc - currind;
	if(tmp > 0){
	  fprintf(ostream, "%*s", tmp, "");
	  currind += tmp;
	  tmp = 0;
	}
	else{
	  fprintf(ostream, " ");
	  tmp = 1;
	}

	//print info if available
	int len = ( _opts->doc?(strlen(_opts->doc)):0 ) +
	  strlen(enddocchar) +
	  ( _opts->def?(strlen(_opts->def) + strlen(defword) + 
			strlen(postdefword)):0 );
	if(len){
	  char *fdoc = (char *)calloc(len+1, sizeof(char));
	  if(_opts->doc)
	    strcpy(fdoc, _opts->doc);
	  if(_opts->def) {
	    strcat(fdoc, defword);
	    strcat(fdoc, _opts->def);
	    strcat(fdoc, postdefword);
	  }
	  if(fdoc[strlen(fdoc)-1] != enddocchar[0])
	    strcat(fdoc, enddocchar);
	  doc = fdoc;
	  while((doc = linepad(out, cols-currind-tmp, doc))){
	    tmp = 0;
	    fprintf(ostream, "%s\n", out);
	    fprintf(ostream, "%*s", indentdoc, "");
	    currind = indentdoc;
	  }
	  free(fdoc);
	}
	else{
	  strncpy(out, noinfo, cols);
	  out[cols] = '\0';
	}
	fprintf(ostream, "%s%s", out, postoption);
      }
      _opts++;

    }
  }

  if(_cfg && _cfg->contact)
    fprintf(ostream, _cfg->contintro?_cfg->contintro:contintro
	    ,_cfg->contact);

  if(_cfg){
    if(_cfg->intro)
      free(intro);
    if(_cfg->intro2)
      free(intro2);
    if(_cfg->options)
      free(options);
    if(_cfg->noinfo)
      free(noinfo);
  }
  free(hprogname);

  procopt_free();

  exit(status);
}


/* \fcnfh
   Mimics getopt_long (assuming always '.flags' equal to 0, but looking
   in line for a field like as 'name value'

   @returns 'val' from field in the struct option that contain the
                  right name
	    '?'   if name is an unknown option
	    ':'   if there is a missing parameter
*/
static int
getoptfrom(char *line,		/* line that should be of the form
				   'name value' */
	   struct option *getopt, /* array of the same form as the one
				       required for getopt_long */
	   int *longindex)	/* if not NULL, it points to a variable
				   which is set to the index of the long
				   option relative to longopts */
{
  char *opt;
  int index=0,flen,ret;

  opt=line;

  while(*opt&&*opt!=' '&&*opt!='\t')
    opt++;
  flen=opt-line;
  while(*opt==' '||*opt=='\t')
    opt++;

  ret='?';

  //search for each parameter
  while(getopt->name||getopt->has_arg||getopt->flag||getopt->val){
    if(getopt->name&&strncmp(getopt->name,line,flen)==0){
      if(longindex) 
	*longindex=index;
      switch(getopt->has_arg){
      case required_argument:
	if(*opt){
	  optarg=opt;
	  ret=getopt->val;
	}
	else
	  ret=':';
	break;
      case no_argument:
	ret=getopt->val;
	break;
      default:
	fprintf(stderr,
		"procopt:: Error in has_arg option at function\n"
		"getoptfrom(), which is called by getopt_long_files\n"
		"or getprocopt. Only required_argument or no_argument\n"
		"are accepted now\n");
	exit(EXIT_FAILURE);
      }
      break;
    }
    index++;
    getopt++;
  }

  return ret;
}


/* \fcnfh
   Frees all the memory allocated. This have to be called after any
   other call to getprocopt. Otherwise, because of the argv reordering
   everything will be mixed up.
 */
void
procopt_free()
{
  cfg_free();

  if(freed){
    fprintf(stderr,
	    "procopt:: getprocopt_free() was called twice\n");
    exit(EXIT_FAILURE);
  }

  //variables that are going to have allocated a space
  free(prgname);
  free(getopts);
  free(shortopts);
  free(fp);
  free(paramfiles);
  if(line) free(line);
}


/* \fcnfh
   Frees all the memory allocated in config.
 */
void
cfg_free()
{
  if(_cfg->freecfg)
    free(_cfg);
}
