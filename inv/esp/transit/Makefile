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
            readlineinfo.c transitstd.c transit.c procopt.c \
	    makesample.c extinction.c readatminfo.c 
PROTO_FILES = $(DEP_FILES)
TR_OBJ    = transitstd.o readlineinfo.o makesample.o \
	    readatminfo.o extinction.o
LR_OBJ    = transitstd.o dbread_pands.o
UT_OBJ    = sampling.o voigt.o procopt.o

UT_HDR    = sampling.h profile.h procopt.h

UT_DEP    = $(UT_OBJ) $(HDR_UT)
TR_DEP    = $(TR_OBJ) transit.h transit.c
LR_DEP    = $(LR_OBJ) lineread.h lineread.c

#make specials
.PHONY=clean dep

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


transit: $(UT_DEP) $(TR_DEP)
	$(CC) $(STDFLAGS) -o $@ $(filter %.c %.o,$^) $(LIBS)
	cp -pf transit $(BIN)

lineread: $(UT_DEP) $(LR_DEP)
	$(CC) $(STDFLAGS) -o $@ $(filter %.c %.o,$^) $(LIBS)
	cp -pf lineread $(BIN)

#$(UT_OBJ): %.o: %.c
#	$(CC) -c $(STDFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(STDFLAGS) -o $@ $<

%.c: %.txc
	$(CLATEX) $< -sl

clean:
	rm -f $(TR_OBJ) $(UT_OBJ) $(DEP_FILES) texput.log transit

dep: $(DEP_FILES)
	makedepend -Y -f .depen -- $(STDFLAGS) -- $^

proto: $(PROTO_FILES:.c=.proto)

%.proto: %.c
	$(CPROTO) $(CPROTOFLAGS) -o $(INCL)/$*_proto.h $<


#/ast/plume/src/plume/libplume
#~/hacking/C   #templates
#-i -m -U __GNUC__

include .depen
