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
#vpath %.c $(SRC)/transit:$(SRC)/lineread:$(SRC)/util

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

#$(OBJ_UT): %.o: %.c
#	$(CC) -c $(STDFLAGS) -o $@ $<

%.o: %.c
	if (!-e $(INC)/$*_proto.h)
		echo "You must run 'make proto' manually first"
		exit 1;
	endif
	$(CC) -c $(STDFLAGS) -o $@ $<

%.c: %.txc
	$(CLATEX) $< -sl

clean:
	rm -f $(OBJ_TR) $(OBJ_UT)

dep: $(DEP_FILES)
	makedepend -Y -- $(STDFLAGS) -- $^

proto: $(PROTO_FILES:.c=.proto)

%.proto: %.c
	$(CPROTO) $(CPROTOFLAGS) -o $(INCL)/$*_proto.h $<


#/ast/plume/src/plume/libplume
#~/hacking/C   #templates
#-i -m -U __GNUC__
# DO NOT DELETE

sampling.o: ./include/transit.h /home/devel/include/P/Putil.h
sampling.o: ./include/util/sampling.h ./include/sampling_proto.h
sampling.o: ./include/util/profile.h ./include/voigt_proto.h
sampling.o: ./include/transit_proto.h ./include/readlineinfo_proto.h
sampling.o: ./include/transitstd_proto.h
voigt.o: ./include/util/profile.h /home/devel/include/P/Putil.h
voigt.o: ./include/transit.h ./include/util/sampling.h
voigt.o: ./include/sampling_proto.h ./include/transit_proto.h
voigt.o: ./include/readlineinfo_proto.h ./include/transitstd_proto.h
voigt.o: ./include/voigt_proto.h
lineread.o: ./include/transit.h /home/devel/include/P/Putil.h
lineread.o: ./include/util/sampling.h ./include/sampling_proto.h
lineread.o: ./include/util/profile.h ./include/voigt_proto.h
lineread.o: ./include/transit_proto.h ./include/readlineinfo_proto.h
lineread.o: ./include/transitstd_proto.h ./include/lineread.h
lineread.o: ./include/lineread_proto.h dbread_pands.c
dbread_pands.o: ./include/transit.h /home/devel/include/P/Putil.h
dbread_pands.o: ./include/util/sampling.h ./include/sampling_proto.h
dbread_pands.o: ./include/util/profile.h ./include/voigt_proto.h
dbread_pands.o: ./include/transit_proto.h ./include/readlineinfo_proto.h
dbread_pands.o: ./include/transitstd_proto.h ./include/lineread.h
dbread_pands.o: ./include/lineread_proto.h
readlineinfo.o: ./include/transit.h /home/devel/include/P/Putil.h
readlineinfo.o: ./include/util/sampling.h ./include/sampling_proto.h
readlineinfo.o: ./include/util/profile.h ./include/voigt_proto.h
readlineinfo.o: ./include/transit_proto.h ./include/readlineinfo_proto.h
readlineinfo.o: ./include/transitstd_proto.h
transitstd.o: ./include/transit.h /home/devel/include/P/Putil.h
transitstd.o: ./include/util/sampling.h ./include/sampling_proto.h
transitstd.o: ./include/util/profile.h ./include/voigt_proto.h
transitstd.o: ./include/transit_proto.h ./include/readlineinfo_proto.h
transitstd.o: ./include/transitstd_proto.h
transit.o: ./include/transit.h /home/devel/include/P/Putil.h
transit.o: ./include/util/sampling.h ./include/sampling_proto.h
transit.o: ./include/util/profile.h ./include/voigt_proto.h
transit.o: ./include/transit_proto.h ./include/readlineinfo_proto.h
transit.o: ./include/transitstd_proto.h
