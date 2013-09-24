/*********************************************************************
 * 
 * procopt_sample.c - Sample file for processing of arguments using
 *                    procopt routine.
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

/* Sample file for a simple call to procopt which will process -h
   or --help to display the help, or --paramf=<file> to search for
   more long options in <file>. */


#include <pu/procopt.h>

#define CFGFILE "argum-default"


int argum(int argc, char **argv)
{

  //General help-option structure
  struct optdocs var_docs[]={
    {NULL,0,HELPTITLE,NULL,
     NULL,"GENERAL ARGUMENTS"},
    {"help",'h',no_argument,NULL,
     NULL,"Prints list of possible parameters"},
    {"paramf",'p',ADDPARAMFILE,NULL,
     "filename","Use filename to read parameters in addition to"
     "default file(s): '" CFGFILE "'"},

    {NULL,0,0,NULL,NULL,NULL}
  };

  struct optcfg var_cfg;
  memset(&var_cfg,0,sizeof(var_cfg));
  var_cfg.contact="My Self <myself@my.web>";
  var_cfg.files=CFGFILE;
  var_cfg.columns=70;


  procopt_debug=1;
  while(1){
    int rn = procopt(argc,argv,var_docs,&var_cfg,NULL);
    if (rn==-1)
      break;

    switch(rn){
    case '?':
      rn=optopt;
      fprintf(stderr,
	      "Unknown, unsupported, or missing parameter to option "
	      "of code %i(%c) passed\n"
	      "as argument, use '-h' to see accepted options.\n"
	      ,rn,(char)rn);
      procopt_free();
      exit(EXIT_FAILURE);
      break;
    default:			//Ask for syntax help
      fprintf(stderr,
	      "Even though option of code %i(%c) had a valid structure\n"
	      "element, it had no switch control statement. File %s\n"
	      "need to be revised.\n"
	      ,rn,(char)rn, __FILE__);
      procopt_free();
      exit(EXIT_FAILURE);
      break;
    case 'h':
      prochelp(EXIT_SUCCESS);
      break;
    }
    
  }

  procopt_free();

  return 0;
}


