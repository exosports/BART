#
# Makefile - prepared for transit.
# Copyright (C) 2004-2006 Patricio Rojo (pato@astro.cornell.edu)
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of version 2 of the GNU General 
# Public License as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


#uncomment for clean compilation without debugging info
#NODEBUG  = 1

SHELL = /bin/sh

srcdir = .
top_srcdir = .
tobasedir = .

SUBDIRS = src

#make specials
.PHONY=clean dep proto


#installation directories
ifndef prefix
prefix = ../../
endif
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${prefix}/share
sysconfdir = ${prefix}/etc
sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var
libdir = ${exec_prefix}/lib
infodir = ${prefix}/info
mandir = ${prefix}/man
includedir = ${prefix}/include

#programs
CPROTO  = /usr/bin/cproto
RANLIB  = ranlib
AR      = ar ru
MAKEDEP = $(CC) -MM
CLATEX  = $(HOME)/bin/clatex.pl

EXEEXT =
scriptdir = ./scripts/
mkinstalldirs = $(scriptdir)mkinstalldirs
CC = ${scriptdir}ccd-gcc.sh

#parameters
INCL        = ./include/
#CPUOPTION   = -mcpu=athlon-xp -march=athlon-xp
#CPUOPTION  = -mtune=athlon-tbird
#CPUOPTION  = -mtune=pentium3
EXTPRM      = -DTRANSIT #-O1 #-static
WARNFLAGS   = -Wall -Winline
#-W -Werror
INCLFLAGS   = -I/usr/X11R6/include -I/home/devel/include \
	      -I$(INCL) -I. $(INCLUDES)
LIBFLAGS    = -L/usr/X11R6/lib -L/home/devel/lib
CPROTOFLAGS = $(INCLFLAGS) -f 3 -m -e -i -X 0

NODBGFLAGS  = 
DBGFLAGS    = -g3 -gdwarf-2
ifdef NODEBUG
STDFLAGS    = $(NODBGFLAGS)
else
STDFLAGS    = $(DBGFLAGS)
endif
STDFLAGS   += $(INCLFLAGS) $(LIBFLAGS) $(WARNFLAGS) -pedantic \
	      -std=c99 $(EXTPRM) $(CPUOPTION)  $(CFLAGS)

#libraries to link
LIBS  = -lm  -lpu

INSTALL = /usr/bin/install -c
INSTALL_PROGRAM = ${INSTALL}
INSTALL_DATA = ${INSTALL} -m 644
install_sh_DATA = $(install_sh) -c -m 644
install_sh_PROGRAM = $(install_sh) -c
INSTALL_SCRIPT = ${INSTALL}
INSTALL_HEADER = $(INSTALL_DATA)

include Beauty.mk
include Rules.mk
