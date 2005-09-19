#installation directories
ifndef prefix
prefix = /usr/local
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
mkinstalldirs = $(BASEBDIR)/mkinstalldirs

#parameters
INCL        = $(BASEBDIR)/include/
CPUOPTION   = -mtune=athlon-xp
#CPUOPTION  = -mtune=athlon-tbird
#CPUOPTION  = -mtune=pentium3
EXTPRM      = -DTRANSIT -static
WARNFLAGS   = -Wall -Winline
GSLPRM      = -D_USE_GSL
GSLLIB      = -lgsl -lblas
#-W -Werror
INCLFLAGS   = -I/usr/X11R6/include -I/home/devel/include -I$(INCL) -I. $(INCLUDES)
LIBFLAGS    = -L/usr/X11R6/lib -L/home/devel/lib
CPROTOFLAGS = $(INCLFLAGS) $(GSLPRM) -f 3 -m -e -i -X 0

NODBGFLAGS  = 
DBGFLAGS    = -g3 -gdwarf-2
ifdef NODEBUG
STDFLAGS    = $(NODBGFLAGS)
else
STDFLAGS    = $(DBGFLAGS)
endif
STDFLAGS   += $(INCLFLAGS) $(LIBFLAGS) $(WARNFLAGS) -pedantic \
	      -std=c99 $(EXTPRM) $(CPUOPTION) $(GSLPRM) $(CFLAGS)

#libraries to link
LIBS  = -lX11 -lpu -lm $(GSLLIB)

include $(BASEBDIR)/beauty.mak

INSTALL = /usr/bin/install -c
INSTALL_PROGRAM = ${INSTALL}
INSTALL_DATA = ${INSTALL} -m 644
install_sh_DATA = $(install_sh) -c -m 644
install_sh_PROGRAM = $(install_sh) -c
INSTALL_SCRIPT = ${INSTALL}
INSTALL_HEADER = $(INSTALL_DATA)
