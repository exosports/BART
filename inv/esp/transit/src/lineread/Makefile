#this has to be defined for every file in the tree
BASEBDIR=/home/pato/inv/esp/prg/transit/
include local.mak
include $(BASEBDIR)/defs.mak

ifeq ($(strip $(wildcard proto_*.h)),)
firstproto=1
endif

install-bin: install-binPROGRAMS
uninstall-bin: uninstall-binPROGRAMS
binPROGRAMS_INSTALL = $(INSTALL_PROGRAM)
install-binPROGRAMS: $(bin_PROGRAMS)
	$(call cmd,mkdir)
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
	$(call cmd,installbin)
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
	-test -z "$(bin_PROGRAMS)" || rm -f $(bin_PROGRAMS)


#Linking of binary
$(bin_PROGRAMS):
	$(call cmd,exec)

#Compilation of object files
%.o: %.c
	$(call cmd,c_o)


#Find dependencies
dep:
	@echo >.depen
	$(call cmd,dep)


#making prototype files
proto: check_cproto $(filter-out "/",$(SOURCES:.c=.proto))
ifdef firstproto
	@echo
	@echo "** Probably reruning 'make proto' will clear any error message above **"
	@echo
endif
%.proto: %.c
	$(call cmd,proto)

#checking that apropiate cproto version is being used
check_cproto:
	@exec $(BASEBDIR)/scripts/check_proto.sh $(CPROTO)


#Clean
cleanall: clean cleanproto
clean: clean_srcbuilddir

cleanproto:
	$(call cmd,clean_proto)

clean_srcbuilddir:
	$(call cmd,cleansrcb)

#include dependencies
-include .depen
