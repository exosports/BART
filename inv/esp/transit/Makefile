#uncomment for clean compilation without debugging info
#NODEBUG  = 1


#directories
SRC=.
BIN=../bin
INCL=./include
INCL_UT=$(INCL)/util

#paths for dependencies
vpath %.h include
vpath %.txc transit:lineread:util
vpath %.c transit:lineread:util

#programs
CPROTO  = /usr/bin/cproto
CLATEX  = /home/pato/bin/clatex.pl
RANLIB  = ranlib
AR      = ar ru
#objects and headers
OBJ_UT= sampling.o voigt.o
OBJ_TR= transitstd.o readlineinfo.o

#make specials
.PHONY=clean dep transit

#parameters
EXTPRM      =
WARNFLAGS   = -Wall -Winline
#-W -Werror
INCLFLAGS   = -I/usr/X11R6/include -I/home/devel/include -I$(INCL)
LIBFLAGS    = -L/usr/X11R6/lib -L/home/devel/lib
CPROTOFLAGS = $(INCLFLAGS) -f 3 -m -e

NODBGFLAGS  = -O2
DBGFLAGS    = -g3 -gdwarf-2
ifdef NODEBUG
STDFLAGS    = $(NODBGFLAGS)
else
STDFLAGS    = $(DBGFLAGS)
endif
STDFLAGS   += $(INCLFLAGS) $(LIBFLAGS) $(WARNFLAGS) -pedantic -std=c99 $(EXTPRM) $(CFLAGS)

#libraries to link
LIBS  = -lX11 -lP -lm



###############################
# $^ dependencies
# $< first dependencie
# $* stem of a pattern in target:pattern:depend
# $@ target
# $? modified dependencies
################################


$(BIN)/transit: $(OBJ_TR) $(OBJ_UT) transit.c
		$(CC) $(STDFLAGS) -o $@ $^ $(LIBS)

transit: $(BIN)/transit
	@touch $<

$(OBJ_UT): %.o: %.c
	$(CPROTO) $(CPROTOFLAGS) -o $(INCL_UT)/$*_proto.h $<
	$(CC) -c $(STDFLAGS) -o $@ $<

$(OBJ_TR): %.o: %.c
	$(CPROTO) $(CPROTOFLAGS) -o $(INCL)/$*_proto.h $<
	$(CC) -c $(STDFLAGS) -o $@ $<

%.c: %.txc
	$(CLATEX) $< -sl
	@sync


#/ast/plume/src/plume/libplume
#~/hacking/C   #templates
#-i -m -U __GNUC__
# DO NOT DELETE

lineread/lineread.o: lineread/dbread_pands.c
