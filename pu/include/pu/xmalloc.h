/*
 * xmalloc.h - Functions that check correct allocation
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

#ifndef _XMALLOC_H
#define _XMALLOC_H

#define malloc(n) xmalloc(n)
#define realloc(p,n) xrealloc(p,n)
#define calloc(n,s) xcalloc(n,s)
#define strdup(p) xstrdup(p)

//functin definition (proto_xmalloc.h)
#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* xmalloc.c */
extern void *xmalloc P_((size_t n));
extern void *xcalloc P_((size_t n, size_t s));
extern void *xrealloc P_((void *p, size_t n));
extern char *xstrdup P_((char *str));

#undef P_


#endif /* _XMALLOC_H */
