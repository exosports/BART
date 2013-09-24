#Adapted from recipes by Emile van Bergen, 2002

#standard things
sp             := $(sp).x
dirstack_$(sp) := $(d)

#subdirectories in random order
#d       := $(d)/test
#include   $(d)/Rules.mk

#binaries in this directory
bin_PROGRAMS_$(d) := readdump$(EXEEXT)
lib_STT_$(d)      := 
lib_DYN_$(d)      := 
lib_$(d)          := 
local_$(d)        := $(lib_$(d)) $(bin_PROGRAMS_$(d))

#the following per binary
readdump_FILES_$(d) := readdump

readdump_OBJS_$(d)      := $(readdump_FILES_$(d):%=$(d)/%.o)
readdump_PIC_OBJS_$(d)  := $(readdump_FILES_$(d):%=$(d)/%_pic.o)
readdump_DEPS_$(d)      := $(readdump_FILES_$(d):%=$(d)/%.d)

$(lib_STT_$(d)): $(readdump_OBJS_$(d))
$(lib_DYN_$(d)): $(readdump_PIC_OBJS_$(d))
readdump$(EXEEXT): $(readdump_OBJS_$(d))

#gather all binaries together
PIC_OBJS_$(d)    := $(readdump_PIC_OBJS_$(d))
OBJS_$(d)        := $(readdump_OBJS_$(d))
DEPS_$(d)        := $(readdump_DEPS_$(d))

CLEAN := $(CLEAN) $(OBJS_$(d)) $(PIC_OBJS_$(d)) $(DEPS_$(d)) \
	 $(local_$(d))
bin_PROGRAMS :=  $(bin_PROGRAMS) $(bin_PROGRAMS_$(d))
lib_DYNAMIC := $(lib_DYNAMIC) $(lib_DYN_$(d))
lib_STATIC := $(lib_STATIC) $(lib_STT_$(d))
PIC_OBJS := $(PIC_OBJS) $(PIC_OBJS_$(d))


#extra dependencies and local variables
$(local_$(d)): $(d)/Rules.mk
$(local_$(d)): CF_LOCAL := -D_USE_GSL -Iinclude \
	`pkg-config --cflags gtk+-2.0`
$(local_$(d)): LL_LOCAL := -lm  `pkg-config --libs gtk+-2.0`

#Standard closure
-include $(DEPS_$(d))

sp   := $(basename $(sp))
d    := $(dirstack_$(sp))
