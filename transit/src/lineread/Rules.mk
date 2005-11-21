#Adapted from recipes by Emile van Bergen, 2002

#standard things
sp             := $(sp).x
dirstack_$(sp) := $(d)

#subdirectories in random order
#d       := $(d)/test
#include   $(d)/Rules.mk

#binaries in this directory
bin_PROGRAMS_$(d) := lineread
lib_STT_$(d)      := 
lib_DYN_$(d)      := 
lib_$(d)          := $(lib_STT_$(d)) $(lib_DYN_$(d))
local_$(d)        := $(lib_$(d)) $(bin_PROGRAMS_$(d))

#the following per binary
lineread_FILES_$(d)    := ../transit/transitstd dbread_pands \
			  lineread dbread_text
lineread_OBJS_$(d)     := $(lineread_FILES_$(d):%=$(d)/%.o)
lineread_PIC_OBJS_$(d) := $(lineread_FILES_$(d):%=$(d)/%_pic.o)
lineread_DEPS_$(d)     := $(lineread_FILES_$(d):%=$(d)/%.o.d)    \
		          $(lineread_FILES_$(d):%=$(d)/%_pic.o.d)

#gather all binaries together
PIC_OBJS_$(d)    := $(lineread_PIC_OBJS_$(d))
OBJS_$(d)        := $(lineread_OBJS_$(d))
DEPS_$(d)        := $(lineread_DEPS_$(d))

CLEAN := $(CLEAN) $(OBJS_$(d)) $(PIC_OBJS_$(d)) $(DEPS_$(d)) $(local_$(d))
bin_PROGRAMS := $(bin_PROGRAMS) $(bin_PROGRAMS_$(d))
lib_DYNAMIC  := $(lib_DYNAMIC)  $(lib_DYN_$(d))
lib_STATIC   := $(lib_STATIC)   $(lib_STT_$(d))
PIC_OBJS     := $(PIC_OBJS)     $(PIC_OBJS_$(d))
OBJS         := $(OBJS)         $(OBJS_$(d))

#extra dependencies and local variables
$(lib_STT_$(d)): $(OBJS_$(d))
$(lib_DYN_$(d)): $(PIC_OBJS_$(d))
lineread$(EXEEXT):  $(lineread_OBJS_$(d))

LOCAL_INC_$(d) := -I$(d) -I$(d)/../transit
$(OBJS_$(d):.o=.proto): CP_LOCAL :=  $(LOCAL_INC_$(d))
$(local_$(d)): $(d)/Rules.mk
$(local_$(d)): CF_LOCAL := -D_USE_GSL $(LOCAL_INC_$(d)) \
	-DHAVE_INLINE -DGSL_RANGE_CHECK_OFF \
	`pkg-config --cflags gtk+-2.0`
$(local_$(d)): LL_LOCAL := -lm -lpu -lgsl `pkg-config --libs gtk+-2.0` \
	-lplplotd -lcfitsio -lblas

#Standard closure
-include $(DEPS_$(d))

sp   := $(basename $(sp))
d    := $(dirstack_$(sp))
