#Check for verbosity
ifdef V
 ifeq ("$(origin V)", "command line")
    BUILD_VERBOSE = $(V)
  endif
endif
ifndef BUILD_VERBOSE
  BUILD_VERBOSE = 0
endif

###############################
# $^ dependencies
# $< first dependencie
# $* stem of a pattern in target:pattern:depend
# $@ target
# $? modified dependencies
################################

#####################
## Beautify output ##
#####################
##Adapted from Kernel 2.6 sources##
#Different levels of output from more to less can be attained by
#specifying V=n where n can be
#  0(reduced output), or 1(full output)
#If there is no especification then use reduced output
# ---------------------------------------------------------------------------
# (from Kernel's Makefile)
# Normally, we echo the whole command before executing it. By making
# that echo $($(quiet)$(cmd)), we now have the possibility to set
# $(quiet) to choose other forms of output instead, e.g.
#
#         quiet_cmd_cc_o_c = Compiling $(RELDIR)/$@
#         cmd_cc_o_c       = $(CC) $(c_flags) -c -o $@ $<
#
# If $(quiet) is empty, the whole command will be printed.
# If it is set to "quiet_", only the short version will be printed. 
# If it is set to "silent_", nothing wil be printed at all, since
# the variable $(silent_cmd_cc_o_c) doesn't exist.
#
# A simple variant is to prefix commands with $(Q) - that's usefull
# for commands that shall be hidden in non-verbose mode.
#
#	$(Q)ln $@ :<
#
# If BUILD_VERBOSE equals 0 then the above command will be hidden.
# If BUILD_VERBOSE equals 1 then the above command is displayed.
ifeq ($(BUILD_VERBOSE),1)
  quiet =
  Q =
else
  quiet=quiet_
  Q = @
  MAKEFLAGS += --no-print-directory
endif

# If the user is running make -s (silent mode), suppress echoing of
# commands
ifneq ($(findstring s,$(MAKEFLAGS)),)
  quiet=silent_
endif


# If quiet is set, only print short version of command
cmd = @$(if $($(quiet)cmd_$(1)),echo ' $($(quiet)cmd_$(1))' &&) $(cmd_$(1))

#now for the different functions
quiet_cmd_proto = Making prototypes for $<
      cmd_proto = $(CPROTO) $(CPROTOFLAGS) $(CP_LOCAL) -o $(dir $*)/proto_$(notdir $*).h $<

quiet_cmd_clean = Deleting non-source files ...
      cmd_clean = rm -f $(CLEAN)

quiet_cmd_cleansrcb = Deleting non-source files ...
      cmd_cleansrcb = rm -f $(OBJECTS) .depen $(PROGRAMS)

quiet_cmd_clean_proto = Deleting prototype $(notdir $*)
      cmd_clean_proto = rm -f $(dir $*)/proto_$(notdir $*).h

quiet_cmd_c_o   = Compiling $@
      cmd_c_o   = $(CC) -c $(CF_LOCAL) $(STDFLAGS) -o $@ $<

quiet_cmd_c_o_pic= Compiling PIC object $@
      cmd_c_o_pic= $(CC) -c $(CF_LOCAL) -fpic $(STDFLAGS) -o $@ $<

quiet_cmd_cp    = Copying to $@
      cmd_cp    = cp -f $< ./$(notdir $(basename $<)).c

quiet_cmd_txc_c = Clatexing to $@
      cmd_txc_c = $(CLATEX) $< -sl

quiet_cmd_exec  = Building executable "$@"
      cmd_exec  = $(CC) $(STDFLAGS) $(CF_LOCAL) -o $@ $(filter %.c %.o,$^) $(LIBS)  $(LL_LOCAL)

quiet_cmd_test  = Building executable "$@" for testing
      cmd_test  = $(CC) $(STDFLAGS) $(CF_LOCAL) $(CF_TEST) -o $@ $(filter %.c %.o,$^) $(LIBS) $(LL_LOCAL)

quiet_cmd_dlib  = Building the dynamic library $@
      cmd_dlib  = $(CC) -shared -Wl,--whole-archive $(filter %.o,$^) \
	-Wl,--no-whole-archive  -lm  -Wl,-soname -Wl,libpu.so.1 \
	-o libpu.so.1  -lm

quiet_cmd_slib  = Building the static library $@
      cmd_slib  = $(AR) $(ARFLAGS) $@ $(filter %.o,$^); touch $@

quiet_cmd_dep   = Making dependencies ...
      cmd_dep   = $(MAKEDEP) $(STDFLAGS) $^ >.depen

quiet_cmd_install = Installing $(INSTALLFROMTHISBIN) ...
      cmd_install = $(INSTALL_PROGRAM_ENV) $(INSTALL_PROGRAM) $(INSTALLFROMTHISBIN) $(INSTALLTOTHISDIR) || echo "=====> Make sure you have the right permission <===="

quiet_cmd_installhead = Installing header $(INSTALLFROMTHISBIN) ...
      cmd_installhead = $(INSTALL_PROGRAM_ENV) $(INSTALL_PROGRAM) -m 0644 $(INSTALLFROMTHISBIN) $(INSTALLTOTHISDIR) || echo "=====> Make sure you have the right permission <===="

quiet_cmd_mkdirbin = Creating "$(DESTDIR)$(bindir)" if it is not there...
      cmd_mkdirbin = $(mkinstalldirs) $(DESTDIR)$(bindir)

quiet_cmd_mkdirlib = Creating "$(DESTDIR)$(libdir)" if it is not there...
      cmd_mkdirlib = $(mkinstalldirs) $(DESTDIR)$(libdir)

quiet_cmd_mkdirhead = Creating "$(DESTDIR)$(includedirpkg)" if it is not there...
      cmd_mkdirhead = $(mkinstalldirs) $(DESTDIR)$(includedirpkg)

quiet_cmd_uninstallbin = Uninstalling $(DELETETHISFILE) ...
      cmd_uninstallbin = rm -f $(DESTDIR)$(bindir)/$(DELETETHISFILE)

quiet_cmd_uninstalllib = Uninstalling $(DELETETHISFILE) ...
      cmd_uninstalllib = rm -f $(DESTDIR)$(libdir)/$(DELETETHISFILE)

quiet_cmd_uninstallhead = Uninstalling $(DELETETHISFILE) ...
      cmd_uninstallhead = rm -f $(DESTDIR)$(includedirpkg)/$(DELETETHISFILE)

makingv         = echo "Processing '$(1)' in $(2)"
