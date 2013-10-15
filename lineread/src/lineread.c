/*
 * lineread.c   - output adequate line information in TLI format.
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
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *
 */


#include <lineread.h>

static struct hints hint;

/* Wavelength units is always microns */
double tli_fct = 1e-4;
const char *tli_fct_name = "microns";


int 
main(int argc, char **argv){
  /* Set prgname FINDME: for what?                      */
  messagep_name("lineread");
  /* Read-in command parameters and populate hint       */
  argum(argc, argv, &hint);
  /* Read the data bases using the corresponding driver */
  db_drivers(&hint);
  /* Free drivers and hint structs                      */
  lineread_free();
  /* Free prgname variable                              */
  messagep_free();

  exit(EXIT_SUCCESS);
}

/* Free drivers, hint, and prgname memory: */
void 
lineread_free(){
  drivers_free(&hint);
  hints_free(&hint);
}
