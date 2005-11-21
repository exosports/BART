#Adapted from recipes by Emile van Bergen, 2002

#standard things
sp             := $(sp).x
dirstack_$(sp) := $(d)

#subdirectories in random order
d       := $(d)/transit
include   $(d)/Rules.mk

d       := $(d)/lineread
include   $(d)/Rules.mk

#binaries in this directory
bin_PROGRAMS_$(d) := 
lib_STT_$(d)      := 
lib_DYN_$(d)      := 
lib_$(d)          := $(lib_STT_$(d)) $(lib_DYN_$(d))
local_$(d)        := $(lib_$(d)) $(bin_PROGRAMS_$(d))

#the following per binary
@P_BINARY@_FILES_$(d) := 

@P_BINARY@_OBJS_$(d)      := $(@P_BINARY@_FILES_$(d):%=$(d)/%.o)
@P_BINARY@_PIC_OBJS_$(d)  := $(@P_BINARY@_FILES_$(d):%=$(d)/%_pic.o)
@P_BINARY@_DEPS_$(d)      := $(@P_BINARY@_FILES_$(d):%=$(d)/%.d)

$(lib_STT_$(d)): $(@P_BINARY@_OBJS_$(d))
$(lib_DYN_$(d)): $(@P_BINARY@_PIC_OBJS_$(d))
@P_BINARY@$(EXEEXT):  $(@P_BINARY@_OBJS_$(d))

#gather all binaries together
PIC_OBJS_$(d)    := $(@P_BINARY@_PIC_OBJS_$(d))
OBJS_$(d)        := $(@P_BINARY@_OBJS_$(d))
DEPS_$(d)        := $(@P_BINARY@_DEPS_$(d))

CLEAN := $(CLEAN) $(OBJS_$(d)) $(PIC_OBJS_$(d)) $(DEPS_$(d)) \
	 $(local_$(d))
bin_PROGRAMS := $(bin_PROGRAMS) $(bin_PROGRAMS_$(d))
lib_DYNAMIC  := $(lib_DYNAMIC)  $(lib_DYN_$(d))
lib_STATIC   := $(lib_STATIC)   $(lib_STT_$(d))
PIC_OBJS     := $(PIC_OBJS)     $(PIC_OBJS_$(d))
OBJS         := $(OBJS)         $(OBJS_$(d))


#extra dependencies and local variables
$(OBJS_$(d):.o=.proto): CP_LOCAL :=  -I$(d)
$(local_$(d)): $(d)/Rules.mk
$(local_$(d)): CF_LOCAL := -D_USE_GSL -I$(d) -Iinclude \
	`pkg-config --cflags gtk+-2.0`
$(local_$(d)): LL_LOCAL := -lm `pkg-config --libs gtk+-2.0`

#Standard closure
#-include $(DEPS_$(d))

sp   := $(basename $(sp))
d    := $(dirstack_$(sp))
