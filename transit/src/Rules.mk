#Adapted from recipes by Emile van Bergen, 2002

#standard things
sp             := $(sp).x
dirstack_$(sp) := $(d)

#subdirectories in random order
#d       := $(d)/atmosphere
#include   $(d)/Rules.mk

#binaries in this directory
bin_PROGRAMS_$(d) := transit
lib_STT_$(d)      := 
lib_DYN_$(d)      := 
lib_$(d)          := $(lib_STT_$(d)) $(lib_DYN_$(d))
local_$(d)        := $(lib_$(d)) $(bin_PROGRAMS_$(d))

#the following per binary
transit_FILES_$(d)    := transitstd readlineinfo makesample    \
                         extinction transit eclipse tau idxrefraction  \
                         argum slantpath geometry observable   \
                         atmosphere/at_file atmosphere/readatm \
			 atmosphere/at_onept cia
transit_OBJS_$(d)     := $(transit_FILES_$(d):%=$(d)/%.o)
transit_PIC_OBJS_$(d) := $(transit_FILES_$(d):%=$(d)/%_pic.o)
transit_DEPS_$(d)     := $(transit_FILES_$(d):%=$(d)/%.o.d)    \
		         $(transit_FILES_$(d):%=$(d)/%_pic.o.d)

#gather all binaries together
PIC_OBJS_$(d)    := $(transit_PIC_OBJS_$(d))
OBJS_$(d)        := $(transit_OBJS_$(d))
DEPS_$(d)        := $(transit_DEPS_$(d))

CLEAN := $(CLEAN) $(OBJS_$(d)) $(PIC_OBJS_$(d)) $(DEPS_$(d)) $(local_$(d))
bin_PROGRAMS := $(bin_PROGRAMS) $(bin_PROGRAMS_$(d))
lib_DYNAMIC  := $(lib_DYNAMIC)  $(lib_DYN_$(d))
lib_STATIC   := $(lib_STATIC)   $(lib_STT_$(d))
PIC_OBJS     := $(PIC_OBJS)     $(PIC_OBJS_$(d))
OBJS         := $(OBJS)         $(OBJS_$(d))

#extra dependencies and local variables
$(lib_STT_$(d)): $(OBJS_$(d))
$(lib_DYN_$(d)): $(PIC_OBJS_$(d))
transit$(EXEEXT):  $(transit_OBJS_$(d))

$(OBJS_$(d):.o=.proto): CP_LOCAL :=  -I$(d)
$(local_$(d)): $(d)/Rules.mk
$(local_$(d)): CF_LOCAL := -D_USE_GSL -I$(d) -Iinclude \
	-DHAVE_INLINE -DGSL_RANGE_CHECK_OFF #\
#	`pkg-config --cflags gtk+-2.0`
$(local_$(d)): LL_LOCAL := -lm -lgsl \
#`pkg-config --libs gtk+-2.0` \
	-lplplotd -lcfitsio -lblas

#Standard closure
-include $(DEPS_$(d))

sp   := $(basename $(sp))
d    := $(dirstack_$(sp))
