#uncomment for clean compilation without debugging info
#NODEBUG  = 1

#directories
ifndef PREFIX
PREFIX  = ..
endif
SRC     = ./src/
BIN     = $(PREFIX)/bin/
INCL    = ./include/
INCL_UT = $(INCL)/util/

#paths for dependencies
vpath %.h include
vpath %.txc $(SRC)/transit:$(SRC)/lineread:$(SRC)/util

#programs
CPROTO  = /usr/bin/cproto
RANLIB  = ranlib
AR      = ar ru
MAKEDEP = $(CC) -MM
CLATEX  = $(HOME)/bin/clatex.pl

#objects and headers
TR_OBJ    = transitstd.o readlineinfo.o makesample.o \
	    readatminfo.o extinction.o transit.o \
	    tau.o idxrefraction.o \
	    slantpath.o
LR_OBJ    = transitstd.o dbread_pands.o lineread.o
UT_OBJ    = sampling.o voigt.o procopt.o iomisc.o \
	    numerical.o xmalloc.o

DEP_FILES = $(TR_OBJ:.o=.c) $(LR_OBJ:.o=.c) $(UT_OBJ:.o=.c)
PROTO_FILES = $(DEP_FILES)

#make specials
.PHONY=clean dep proto

#parameters
CPUOPTION   = -mcpu=athlon-tbird
#CPUOPTION   = -mcpu=pentium3
EXTPRM      =
WARNFLAGS   = -Wall -Winline
GSLPRM      = -D_USE_GSL
GSLLIB      = -lgsl -lblas
#-W -Werror
INCLFLAGS   = -I/usr/X11R6/include -I/home/devel/include -I$(INCL)
LIBFLAGS    = -L/usr/X11R6/lib -L/home/devel/lib
CPROTOFLAGS = $(INCLFLAGS) $(GSLPRM) -f 3 -m -e -i -X 0

NODBGFLAGS  = -O2
DBGFLAGS    = -g3 -gdwarf-2
ifdef NODEBUG
STDFLAGS    = $(NODBGFLAGS)
else
STDFLAGS    = $(DBGFLAGS)
endif
STDFLAGS   += $(INCLFLAGS) $(LIBFLAGS) $(WARNFLAGS) -pedantic \
	      -std=c99 $(EXTPRM) $(CPUOPTION) $(GSLPRM) $(CFLAGS)

#libraries to link
LIBS  = -lX11 -lm $(GSLLIB)



###############################
# $^ dependencies
# $< first dependencie
# $* stem of a pattern in target:pattern:depend
# $@ target
# $? modified dependencies
################################

default: transit lineread

dbgsample: $(UT_OBJ) transitstd.c makesample.c
	$(CC) $(STDFLAGS) -DDBGSAMPLE -o $@ $(filter %.c %.o,$^) $(LIBS)
	cp -pf dbgsample $(BIN)

dbgreadline: $(UT_OBJ) transitstd.c readlineinfo.c
	$(CC) $(STDFLAGS) -DDBGREADLINEINFO -o $@ $(filter %.c %.o,$^) $(LIBS)
	cp -pf dbgreadline $(BIN)

transit: $(UT_OBJ) $(TR_OBJ)
	$(CC) $(STDFLAGS) -o $@ $(filter %.c %.o,$^) $(LIBS)
	cp -pf transit $(BIN)

lineread: $(UT_OBJ) $(LR_OBJ)
	$(CC) $(STDFLAGS) -o $@ $(filter %.c %.o,$^) $(LIBS)
	cp -pf lineread $(BIN)

#$(UT_OBJ): %.o: %.c
#	$(CC) -c $(STDFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(STDFLAGS) -o $@ $<

#if clatex exist then run it, otherwise just copy the file changing
#extensions
%.c: %.txc
ifeq ($(strip $(wildcard $(CLATEX))),)
	cp -f $< ./$(notdir $(basename $<)).c
else
	$(CLATEX) $< -sl
endif

clean:
	rm -f $(TR_OBJ) $(UT_OBJ) $(DEP_FILES) texput.log transit

dep: $(DEP_FILES)
	echo >.depen
	$(MAKEDEP) $(STDFLAGS) $^ >.depen

proto: $(PROTO_FILES:.c=.proto)

%.proto: %.c
	$(CPROTO) $(CPROTOFLAGS) -o $(INCL)/$*_proto.h $<


#/ast/plume/src/plume/libplume
#~/hacking/C   #templates
#-i -m -U __GNUC__

-include .depen
