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

BASEBDIR=/home/pato/inv/esp/prg/transit/
include $(BASEBDIR)/defs.mak

SHELL = /bin/sh

srcdir = .
top_srcdir = .
tobasedir = .

SUBDIRS = src

#make specials
.PHONY=clean dep proto


###############################
# $^ dependencies
# $< first dependencie
# $* stem of a pattern in target:pattern:depend
# $@ target
# $? modified dependencies
################################

RECURSIVE_TARGETS = all-recursive dep-recursive clean-recursive \
	install-bin-recursive uninstall-bin-recursive proto-recursive \
	cleanproto-recursive cleanall-recursive

all:
	$(Q)$(MAKE) all-recursive

$(RECURSIVE_TARGETS):
	$(Q)target=`echo $@ | sed s/-recursive//`; \
	list='$(SUBDIRS)'; \
	cwd=`pwd|sed 's/$(subst /,\/,$(BASEBDIR:/=))/./'`; \
	for subdir in $$list; do \
	  $(call makingv,"$$target","$$cwd/$$subdir");\
	  (cd $$subdir && $(MAKE) $$target) || fail="yes"; \
	done; \
	test -z "$$fail"

#find recursiveness
dep: dep-recursive
install-bin: install-bin-recursive
uninstall-bin: uninstall-bin-recursive
proto: proto-recursive
cleanproto: cleanproto-recursive
clean: clean-recursive
cleanall: cleanall-recursive

install: install-bin
uninstall: uninstall-bin


#explain makefile options regarding printing
help:
	@echo "Make option OUTP specifies output presented: if it is not"
	@echo "specified, then use nice printing, if 1 then no output (only "
	@echo "errors and warnings), and if 2 do a full call"


#install executables in target directory
#	cp -pf lineread $(BIN)
#	cp -pf transit $(BIN)



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



#delete unused files

#cleanproto:
#	$(call cmd,clean_proto)
#cleanall: clean cleanproto


#if clatex exist then run it, otherwise just copy the file changing
#extensions
#%.c: %.txc
#ifeq ($(strip $(wildcard $(CLATEX))),)
#	$(call cmd,cp)
#else
#	$(call cmd,txc_c)
#endif

