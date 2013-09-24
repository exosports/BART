
# Copyright (C) 1994, 1995-8, 1999 Free Software Foundation, Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

#  uncomment the following if testgtk requires the math library
#testgtk_LDADD=-lm


SHELL = /bin/sh

srcdir = .
top_srcdir = .

prefix = /usr/local
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
oldincludedir = /usr/include

DESTDIR =

pkgdatadir = $(datadir)/@PACKAGE@
pkglibdir = $(libdir)/@PACKAGE@
pkgincludedir = $(includedir)/@PACKAGE@

top_builddir = .

ACLOCAL = @ACLOCAL@
AUTOCONF = @AUTOCONF@
AUTOMAKE = @AUTOMAKE@
AUTOHEADER = @AUTOHEADER@

INSTALL = /usr/bin/install -c
INSTALL_PROGRAM = ${INSTALL} $(AM_INSTALL_PROGRAM_FLAGS)
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_SCRIPT = ${INSTALL}
transform = s,x,x,

NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
host_alias = 
host_triplet = @host@
CCR = gcc
HAVE_LIB = @HAVE_LIB@
ISODATE = @ISODATE@
LIB = @LIB@
LTLIB = @LTLIB@
MAKEINFO = @MAKEINFO@
PACKAGE = @PACKAGE@
VERSION = @VERSION@

NROFF = nroff
TAR = tar
GZIP_ENV = --best
DVIPS = dvips

distdir = $(PACKAGE)-$(VERSION)
top_distdir = $(distdir)

mkinstalldirs = ./mkinstalldirs
includedirpkg = $(includedir)/pu

CC = scripts/ccd-gcc.sh

ARFLAGS     = ruc
INCLFLAGS   = -I/usr/X11R6/include -I. $(INCLUDES)
LIBFLAGS    = -L/usr/X11R6/lib
EXTPRM      = -DTESTGTK
WARNFLAGS   = -Wall -Winline

NODBGFLAGS  = 
DBGFLAGS    = -g3 -gdwarf-2
ifdef NODEBUG
STDFLAGS    = $(NODBGFLAGS)
else
STDFLAGS    = $(DBGFLAGS)
endif
STDFLAGS   += $(INCLFLAGS) $(LIBFLAGS) $(WARNFLAGS) -pedantic \
	      -std=c99 $(EXTPRM) $(CPUOPTION) $(GSLPRM) $(CFLAGS)


include Beauty.mk
include Rules.mk

# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
