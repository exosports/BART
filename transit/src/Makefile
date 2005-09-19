#these 2 lines have to be defined for every file in the tree
BASEBDIR=/home/pato/inv/esp/prg/transit/
include $(BASEBDIR)/defs.mak

#Prepare recUrsive
SUBDIRS = transit lineread
RECURSIVE_TARGETS = all-recursive dep-recursive clean-recursive \
	install-bin-recursive uninstall-bin-recursive proto-recursive \
	cleanproto-recursive cleanall-recursive
all:	all-recursive
dep:	dep-recursive
clean:	clean-recursive
cleanall: cleanall-recursive
cleanproto: cleanproto-recursive
install-bin: install-bin-recursive
uninstall-bin: uninstall-bin-recursive
proto:	proto-recursive

SHELL=/bin/sh

$(RECURSIVE_TARGETS):
	$(Q)target=`echo $@ | sed s/-recursive//`; \
	list='$(SUBDIRS)'; \
	cwd=`pwd|sed 's/$(subst /,\/,$(BASEBDIR:/=))/./'`; \
	for subdir in $$list; do \
	  $(call makingv,"$$target","$$cwd/$$subdir");\
	  (cd $$subdir && $(MAKE) $$target) || fail="yes"; \
	done;\
	test -z "$$fail"
