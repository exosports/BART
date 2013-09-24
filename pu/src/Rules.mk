#Adapted from recipes by Emile van Bergen, 2002

#standard things
sp             := $(sp).x
dirstack_$(sp) := $(d)

#subdirectories in random order
#d       := $(d)/test
#include   $(d)/Rules.mk

#binaries in this directory
bin_PROGRAMS_$(d) := 
lib_STT_$(d)      := libpu.a
lib_DYN_$(d)      := libpu.so.1
lib_$(d)          := $(lib_STT_$(d)) $(lib_DYN_$(d))
local_$(d)        := $(lib_$(d)) $(bin_PROGRAMS_$(d))

#the following per binary
libpu_FILES_$(d) := sampling voigt iomisc procopt numerical xmalloc messagep

libpu_OBJS_$(d)      := $(libpu_FILES_$(d):%=$(d)/%.o)
libpu_PIC_OBJS_$(d)  := $(libpu_FILES_$(d):%=$(d)/%_pic.o)
libpu_DEPS_$(d)      := $(libpu_FILES_$(d):%=$(d)/%.d)

$(lib_STT_$(d)): $(libpu_OBJS_$(d))
$(lib_DYN_$(d)): $(libpu_PIC_OBJS_$(d))
libpu$(EXEEXT):  $(libpu_OBJS_$(d))

#gather all binaries together
PIC_OBJS_$(d)    := $(libpu_PIC_OBJS_$(d))
OBJS_$(d)        := $(libpu_OBJS_$(d))
DEPS_$(d)        := $(libpu_DEPS_$(d))

CLEAN := $(CLEAN) $(OBJS_$(d)) $(PIC_OBJS_$(d)) $(DEPS_$(d)) \
	 $(local_$(d))
bin_PROGRAMS :=  $(bin_PROGRAMS) $(bin_PROGRAMS_$(d))
lib_DYNAMIC := $(lib_DYNAMIC) $(lib_DYN_$(d))
lib_STATIC := $(lib_STATIC) $(lib_STT_$(d))
PIC_OBJS := $(PIC_OBJS) $(PIC_OBJS_$(d))


#extra dependencies and local variables
$(local_$(d)): $(d)/Rules.mk
$(local_$(d)): CF_LOCAL := -D_USE_GSL -Iinclude
$(local_$(d)): LL_LOCAL := -lm

#Standard closure
-include $(DEPS_$(d))

sp   := $(basename $(sp))
d    := $(dirstack_$(sp))
