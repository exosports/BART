SRC=.
BIN=../bin
LINC=./include
LIBRARY=-L/home/devel/lib
INCLUDE=-I/home/devel/include -I./include
CDOCSTRIP=/home/pato/bin/clatex.pl
CPROTO=/usr/bin/cproto
EXTPRM=
STAT=-static
SRC_TR=$(SRC)/transit
SRC_UT=$(SRC)/util
SRC_LR=$(SRC)/lineread
LINC_TR=$(LINC)
LINC_UT=$(LINC)/util
LINC_LR=$(LINC)

LIB=-L/usr/X11R6/lib $(LIBRARY)
INC=-I/usr/X11R6/include $(INCLUDE)


WARN=-Wall -Winline
#-W -Werror
STDFLAG=$(LIB) $(INC) $(WARN) -pedantic -std=c99 -O2 $(EXTPRM) $(CFLAGS)
DBGFLAG=$(LIB) $(INC) $(WARN) -pedantic -std=c99 -g3 -gdwarf-2 $(EXTPRM) $(CFLAGS)
STDL=-lP -lm
STDLX=-LX11 -lP -lm
CPROTOFLAGS=-f 3 -m -e

readlineinfo: $(SRC_TR)/transitstd.doc.c
	$(CDOCSTRIP) $(SRC_TR)/$@.doc.c -sl
	$(CPROTO) $(CPROTOFLAGS) $(INC) -o $(LINC_TR)/$@_proto.h $(SRC_TR)/$@.c

transitstd: $(SRC_TR)/transitstd.doc.c
	$(CDOCSTRIP) $(SRC_TR)/$@.doc.c -sl
	$(CPROTO) $(CPROTOFLAGS) $(INC) -o $(LINC_TR)/$@_proto.h $(SRC_TR)/$@.c

sampling: $(SRC_UT)/sampling.doc.c transitstd
	$(CDOCSTRIP) $(SRC_UT)/$@.doc.c -sl
	$(CPROTO) $(CPROTOFLAGS) $(INC) -o $(LINC_UT)/$@_proto.h $(SRC_UT)/$@.c

profile: $(SRC_UT)/voigt.doc.c transitstd
	$(CDOCSTRIP) $(SRC_UT)/voigt.doc.c -sl
	$(CPROTO) $(CPROTOFLAGS) $(INC) -o $(LINC_UT)/$@_proto.h $(SRC_UT)/voigt.c

%: $(SRC_TR)/%.doc.c $(LINC_TR)/transit.h transitstd
	$(CDOCSTRIP) $(SRC_TR)/$@.doc.c -sl
	$(CPROTO) $(CPROTOFLAGS) $(INC) -o $(LINC_TR)/$@_proto.h $(SRC_TR)/$@.c
	$(CC) $(STDFLAG) -o $(BIN)/$@ $(SRC_TR)/$@.c $(STDL)

readline: $(SRC_RL)/readline.doc.c $(LINC_TR)/transit.h transitstd
	$(CDOCSTRIP) $(SRC_RL)/$@.doc.c -sl
	$(CPROTO) $(CPROTOFLAGS) $(INC) -o $(LINC_RL)/$@_proto.h $(SRC_RL)/$@.c
	$(CC) $(STDFLAG) -o $(BIN)/$@ $(SRC_TR)/$@.c $(STDL)

transit: $(SRC_TR)/transit.doc.c $(LINC_TR)/transit.h transitstd readlineinfo sampling profile
	$(CDOCSTRIP) $(SRC_TR)/$@.doc.c -sl
	$(CPROTO) $(CPROTOFLAGS) $(INC) -o $(LINC_TR)/$@_proto.h $(SRC_TR)/$@.c
	$(CC) $(DBGFLAG) -o $(BIN)/$@ $(SRC_TR)/$@.c $(STDL)

dbgspline: sampling
	$(CC) -DDBGSPLINE $(DBGFLAG) -o $(BIN)/dbgspline $(SRC_UT)/sampling.c $(STDL)

dbgsamp: sampling
	$(CC) -DDBGSAMPLING $(DBGFLAG) -o $(BIN)/dbgsamp $(SRC_UT)/sampling.c $(STDL)

dbgpands: $(SRC_RL)/dbread_pands.c
	$(CC) -DDBGPANDS $(DBGFLAG) -o $(BIN)/dbgpands $(SRC_RL)/dbread_pands.c $(STDL)

dbgreadline: readlineinfo
	$(CC) -DDBGREADLINEINFO $(DBGFLAG) -o $(BIN)/dbgreadline $(SRC_TR)/readlineinfo.c $(STDL)

dbgvoigt: $(SRC_UT)/voigt.doc.c
	$(CDOCSTRIP) $(SRC_UT)/voigt.doc.c -sl
	$(CPROTO) $(INC) -e -m $(SRC_UT)/voigt.c -o $(LINC)/transit/profile_proto.h
	$(CC) -DDBGVOIGT $(DBGFLAG) -o $(BIN)/dbgvoigt $(SRC_UT)/voigt.c $(STDL)

transitf: $(SRC_TR)/transit.doc.c $(LINC)/transit/transit.h transitstd readlineinfo sampling profile
	$(CC) -DNODEBUG_TRANSIT $(STDFLAG) -o $(BIN)/transit $(SRC_TR)/transit.c $(STDL)

#/ast/plume/src/plume/libplume
#~/hacking/C   #templates
#-i -m -U __GNUC__
