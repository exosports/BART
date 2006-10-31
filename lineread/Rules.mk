#Adapted from recipes by Emile van Bergen, 2002

#standard things
.SUFFIXES:
.SUFFIXES: .c .o

.PHONY: clean all install

all: binaries libraries

#subdirectories, random order
d   := src
include $(d)/Rules.mk

d   := include
include $(d)/Rules.mk

#Commands

binaries: $(bin_PROGRAMS)
libraries: $(lib_STATIC) $(lib_DYNAMIC)
install: install-bin install-magic #install-header install-lib
uninstall: uninstall-bin uninstall-magic #uninstall-lib uninstall-header

#libraries and binaries
$(lib_STATIC):
	$(call cmd,slib)

$(lib_DYNAMIC):
	$(call cmd,dlib)

$(bin_PROGRAMS): 
	$(call cmd,exec)

$(test_PROGRAMS):
	$(call cmd,test)

#compilation of object files
$(PIC_OBJS): %_pic.o: %.c
	$(call cmd,c_o_pic)

%.o: %.c
	$(call cmd,c_o)

#deletion
clean:
	$(call cmd,clean)
cleanproto:  $(OBJS:.o=.delproto)
%.delproto:
	$(call cmd,clean_proto)
cleanall: clean cleanproto



install-header: install-headerPROGRAMS
uninstall-header: uninstall-headerPROGRAMS
install-headerPROGRAMS:
	$(call cmd,mkdirhead)
	@list='$(header_PROGRAM)'; for p in $$list; do \
	  if test -f $$p ; then \
	    f=`echo "$$p" | sed -e 's,^.*/,,'`; \
	    $(MAKE) INSTALLFROMTHISBIN=$$p \
	      INSTALLTOTHISDIR=$(DESTDIR)$(includedirpkg)/$$f \
		install-header-hook || fail="yes" ; \
	  else :; fi; \
	done;\
	test -z "$$fail"
install-header-hook:
	$(call cmd,installhead)
uninstall-headerPROGRAMS:
	@$(NORMAL_UNINSTALL)
	@list='$(header_PROGRAM)'; for p in $$list; do \
	  f=`echo "$$p" | sed -e 's,^.*/,,'`; \
	  $(MAKE) DELETETHISFILE=$$f uninstall-header-hook \
	    || fail="yes" ; \
	done; \
	test -z "$$fail"
uninstall-header-hook:
	$(call cmd,uninstallhead)
install-lib: install-libPROGRAMS
uninstall-lib: uninstall-libPROGRAMS
install-libPROGRAMS: $(lib_STATIC) $(lib_DYNAMIC)
	$(call cmd,mkdirlib)
	@list='$(lib_STATIC) $(lib_DYNAMIC)'; for p in $$list; do \
	  if test -f $$p \
	  ; then \
	    f=`echo "$$p" | sed -e 's,^.*/,,'`; \
	    $(MAKE) INSTALLFROMTHISBIN=$$p \
	      INSTALLTOTHISDIR=$(DESTDIR)$(libdir)/$$f install-lib-hook \
	      || fail="yes" ; \
	  else :; fi; \
	done;\
	test -z "$$fail"
install-lib-hook:
	$(call cmd,install)
uninstall-libPROGRAMS:
	@$(NORMAL_UNINSTALL)
	@list='$(lib_STATIC) $(lib_DYNAMIC)'; for p in $$list; do \
	  f=`echo "$$p" | sed -e 's,^.*/,,'`; \
	  $(MAKE) DELETETHISFILE=$$f uninstall-lib-hook \
	    || fail="yes" ; \
	done; \
	test -z "$$fail"
uninstall-lib-hook:
	$(call cmd,uninstalllib)
install-bin: install-binPROGRAMS
uninstall-bin: uninstall-binPROGRAMS
install-magic:
	@./scripts/magicadd.sh
uninstall-magic:
	@./scripts/magicadd.sh -u
install-binPROGRAMS: $(bin_PROGRAMS)
	$(call cmd,mkdirbin)
	@list='$(bin_PROGRAMS)'; for p in $$list; do \
	  p1=`echo $$p|sed 's/$(EXEEXT)$$//'`; \
	  if test -f $$p \
	  ; then \
	    p1=`echo "$$p1" | sed -e 's,^.*/,,'`; \
	    f=`echo $$p1|sed '$(transform);s/$$/$(EXEEXT)/'`; \
	    $(MAKE) INSTALLFROMTHISBIN=$$p \
	      INSTALLTOTHISDIR=$(DESTDIR)$(bindir)/$$f install-bin-hook \
	      || fail="yes" ; \
	  else :; fi; \
	done;\
	test -z "$$fail"
install-bin-hook:
	$(call cmd,install)
uninstall-binPROGRAMS:
	@$(NORMAL_UNINSTALL)
	@list='$(bin_PROGRAMS)'; for p in $$list; do \
	  f=`echo $$p|sed 's/$(EXEEXT)$$//;$(transform);s/$$/$(EXEEXT)/'`; \
	  f=`echo "$$f" | sed -e 's,^.*/,,'`; \
	  $(MAKE) DELETETHISFILE=$$f uninstall-bin-hook \
	    || fail="yes" ; \
	done; \
	test -z "$$fail"
uninstall-bin-hook:
	$(call cmd,uninstallbin)
clean-binPROGRAMS:
	@-test -z "$(bin_PROGRAMS)" || rm -f $(bin_PROGRAMS)


#making prototype files
proto: check_cproto $(OBJS:.o=.proto)
	@echo
	@echo "** Reruning 'make proto' should clear any error message above, if any **"
	@echo
%.proto: %.c
	$(call cmd,proto)

#checking that apropiate cproto version is being used
check_cproto:
	@exec ./scripts/check_proto.sh $(CPROTO)


#closure
.SECONDARY: $(CLEAN)
