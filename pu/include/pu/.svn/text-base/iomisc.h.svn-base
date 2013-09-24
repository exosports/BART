/*
 * iomisc.h - Miscellaneous input/output utilities header file
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

#ifndef _IOMISC_H
#define _IOMISC_H

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <termios.h>
#include <stdarg.h>
#include <string.h>

#define __PADTO_COLUMNS 78

extern char *linepad_break;
extern int printpad_columns;
extern char printpad_endpadchar;

#define readd(fp,cp) readds(fp,cp,NULL,0)


#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* iomisc.c */
/* src/iomisc.c */
extern int ncharchg(char *str, char car, char chg);
extern int nchar(char *str, char car);
extern char *readstr_sp_alloc(char *line, char **next, char fspace);
extern void freetoolongerr(void);
extern void settoolongerr(void (*errfcn)(int, char *, int), char *filename, long *currline);
extern int getad(int n, char sep, char *str, double **array);
extern int getnd(int n, char sep, char *str, ...);
extern int getnl(int n, char sep, char *str, ...);
extern void fprintpad(FILE *fp, int indent, char *fmt, ...);
extern double readds(FILE *fp, char *c, char *string, int maxstring);
extern double getds(char *in, char *c, char *string, int maxstring);
extern long readl(FILE *fp, char *c);
extern char *linepad(char *out, int nc, char *in);
extern double askforposd(char *fmt, ...);
extern long askforposl(char *fmt, ...);
extern char *fgets_alloc(FILE *fp, int *max);
extern void splitnzero_add(char ***array, char *string, char sep);
extern char **splitnzero_alloc(char *string, char sep);
extern void splitnzero_free(char **multi);
extern long countfields(char *l, char sep);
extern char fgetupto_err(char *line, int max, FILE *fp, void (*errfcn)(int,char *, int), char *name, long curr);
extern char fgetupto(char *line, int max, FILE *fp);

#undef P_

#endif /* _IOMISC_H */
