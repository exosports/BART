#
# Makefile - prepared for transit.
# Copyright (C) 2004 Patricio Rojo (pato@astro.cornell.edu)
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
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

#directories
ifndef PREFIX
PREFIX  = ..
endif
SRC     = ./src/
BIN     = $(PREFIX)/bin/
INCL    = ./include/

#paths for dependencies
vpath %.h include
vpath %.txc $(SRC)/transit:$(SRC)/transit/atmosphere:\
	    $(SRC)/lineread

#programs
CPROTO  = /usr/bin/cproto
RANLIB  = ranlib
AR      = ar ru
MAKEDEP = $(CC) -MM
CLATEX  = $(HOME)/bin/clatex.pl

#objects and headers
TR_OBJ    = transitstd.o readlineinfo.o makesample.o \
	    extinction.o transit.o \
	    tau.o idxrefraction.o argum.o \
	    slantpath.o geometry.o observable.o \
	    at_file.o readatm.o at_onept.o
LR_OBJ    = transitstd.o dbread_pands.o lineread.o
TEST_FCN  = slantpath dbread_pands

TEST_BIN    = $(patsubst %,test_%,$(TEST_FCN))
PROTO_FILES = $(TR_OBJ:.o=.c) $(LR_OBJ:.o=.c) \
	      $(patsubst %,%.c,$(TEST_BIN))
DEP_FILES   = $(PROTO_FILES) $(TEST_OBJ:.o=.c)

#make specials
.PHONY=clean dep proto

#parameters
CPUOPTION   = -mcpu=athlon-xp
#CPUOPTION   = -mcpu=athlon-tbird
#CPUOPTION   = -mcpu=pentium3
EXTPRM      = -DTRANSIT -static
WARNFLAGS   = -Wall -Winline
GSLPRM      = -D_USE_GSL
GSLLIB      = -lgsl -lblas
#-W -Werror
INCLFLAGS   = -I/usr/X11R6/include -I/home/devel/include -I$(INCL)
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


######################
## Beautiful output ##
######################
##Adapted from Kernel 2.6 sources##
#Different levels of output from more to less can be attained by
#specifying 'make OUTP=n <label>' where n can be
#  1(no output); or 2(regular output).
#If there is no especification then use reduced output
ifndef OUTP
 quiet = quiet_
else 
 ifeq ($(OUTP),1)
  quiet = silen_
 else
  ifeq ($(OUTP),2)
   quiet =
  endif
 endif
endif


# If quiet is set, only print short version of command
cmd = @$(if $($(quiet)cmd_$(1)),echo ' $($(quiet)cmd_$(1))' &&) $(cmd_$(1))

ifeq ($(strip $(wildcard include/proto_*.h)),)
	firstproto = 1
endif

#now for the different functions
quiet_cmd_c_o   = Compiling $@
      cmd_c_o   = $(CC) -c $(STDFLAGS) -o $@ $<

quiet_cmd_cp    = Copying to $@
      cmd_cp    = cp -f $< ./$(notdir $(basename $<)).c

quiet_cmd_txc_c = Clatexing to $@
      cmd_txc_c = $(CLATEX) $< -sl

quiet_cmd_exec  = Building executable "$@"
      cmd_exec  = $(CC) $(STDFLAGS) -o $@ $(filter %.c %.o,$^) $(LIBS)

quiet_cmd_test  = Building for testing executable "$@"
      cmd_test  = $(CC) $(STDFLAGS) -o $@ $(filter %.c %.o,$^) $(LIBS)

quiet_cmd_proto = Making prototypes for $<
      cmd_proto = $(CPROTO) $(CPROTOFLAGS) -o $(INCL)/proto_$*.h $<

quiet_cmd_dep   = Making dependencies ...
      cmd_dep   = $(MAKEDEP) $(STDFLAGS) $^ >.depen

quiet_cmd_clean = Deleting non-source files ...
      cmd_clean = rm -f $(TR_OBJ) $(LR_OBJ) $(DEP_FILES) \
                  texput.log transit lineread .depen

quiet_cmd_clean_proto = Deleting prototypes ...
      cmd_clean_proto = rm -f include/proto_*.h

###############################
# $^ dependencies
# $< first dependencie
# $* stem of a pattern in target:pattern:depend
# $@ target
# $? modified dependencies
################################


#make main executable
default: transit lineread


#explain makefile options regarding printing
help:
	@echo "Make option OUTP specifies output presented: if it is not"
	@echo "specified, then use nice printing, if 1 then no output (only "
	@echo "errors and warnings), and if 2 do a full call"


#install executables in target directory
install: default
	cp -pf lineread $(BIN)
	cp -pf transit $(BIN)


#compile test programs
test:	$(TEST_BIN)
$(TEST_BIN): % : %.o transitstd.o
	$(call cmd,test)


#compile debugging executables
dbgsample: transitstd.c makesample.c
	$(CC) $(STDFLAGS) -DDBGSAMPLE -o $@ $(filter %.c %.o,$^) $(LIBS)
dbgreadline: transitstd.c readlineinfo.c
	$(CC) $(STDFLAGS) -DDBGREADLINEINFO -o $@ $(filter %.c %.o,$^) $(LIBS)
	cp -pf dbgreadline $(BIN)


#make executables
transit : $(TR_OBJ)
lineread: $(LR_OBJ)
transit lineread:
	$(call cmd,exec)


#compilation of object files
%.o: %.c
	$(call cmd,c_o)


#if clatex exist then run it, otherwise just copy the file changing
#extensions
%.c: %.txc
ifeq ($(strip $(wildcard $(CLATEX))),)
	$(call cmd,cp)
else
	$(call cmd,txc_c)
endif


#delete unused files
clean:
	$(call cmd,clean)
cleanproto:
	$(call cmd,clean_proto)
cleanall: clean cleanproto


#find dependencies
dep: $(DEP_FILES)
	@echo >.depen
	$(call cmd,dep)


#making prototype files
proto: check_cproto $(PROTO_FILES:.c=.proto)
ifdef firstproto
	@echo
	@echo "** Probably reruning 'make proto' will clear any error message above **"
	@echo
endif
%.proto: %.c
	$(call cmd,proto)


#checking that apropiate cproto version is being used
check_cproto:
	@exec ./scripts/check_proto.sh $(CPROTO)


#include dependencies
-include .depen

