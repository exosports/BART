#uncomment for clean compilation without debugging info
#NODEBUG  = 1


#directories
SRC=./src
BIN=../bin
INCL=./include
INCL_UT=$(INCL)/util

#paths for dependencies
vpath %.h include
vpath %.txc $(SRC)/transit:$(SRC)/lineread:$(SRC)/util

#programs
CPROTO  = /usr/bin/cproto
CLATEX  = /home/pato/bin/clatex.pl
RANLIB  = ranlib
AR      = ar ru

#objects and headers
DEP_FILES = sampling.c voigt.c lineread.c dbread_pands.c \
            readlineinfo.c transitstd.c transit.c 
PROTO_FILES = $(DEP_FILES)
OBJ_UT    = sampling.o voigt.o
OBJ_TR    = transitstd.o readlineinfo.o

#make specials
.PHONY=clean dep transit

#parameters
EXTPRM      =
WARNFLAGS   = -Wall -Winline
#-W -Werror
INCLFLAGS   = -I/usr/X11R6/include -I/home/devel/include -I$(INCL)
LIBFLAGS    = -L/usr/X11R6/lib -L/home/devel/lib
CPROTOFLAGS = $(INCLFLAGS) -f 3 -m -e -i

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


transit: $(OBJ_TR) $(OBJ_UT) transit.c
	$(CC) $(STDFLAGS) -o $@ $^ $(LIBS)
	cp -pf transit $(BIN)

#$(OBJ_UT): %.o: %.c
#	$(CC) -c $(STDFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(STDFLAGS) -o $@ $<

%.c: %.txc
	$(CLATEX) $< -sl

clean:
	rm -f $(OBJ_TR) $(OBJ_UT) $(DEP_FILES) texput.log transit

dep: $(DEP_FILES)
	makedepend -Y -f .depen -- $(STDFLAGS) -- $^

proto: $(PROTO_FILES:.c=.proto)

%.proto: %.c
	$(CPROTO) $(CPROTOFLAGS) -o $(INCL)/$*_proto.h $<


#/ast/plume/src/plume/libplume
#~/hacking/C   #templates
#-i -m -U __GNUC__

include .depen
