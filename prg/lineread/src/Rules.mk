#Adapted from recipes by Emile van Bergen, 2002

#standard things
sp             := $(sp).x
dirstack_$(sp) := $(d)

#subdirectories in random order
#d       := $(d)/test
#include   $(d)/Rules.mk

#binaries in this directory
bin_PROGRAMS_$(d) := lineread
test_PROGRAMS_$(d) := test_pands test_text
lib_STT_$(d)      := 
lib_DYN_$(d)      := 
lib_$(d)          := $(lib_STT_$(d)) $(lib_DYN_$(d))
local_$(d)        := $(lib_$(d)) $(bin_PROGRAMS_$(d)) $(test_PROGRAMS_$(d))

#the following per binary
lineread_FILES_$(d)    := lineread argum drivers messagep \
			  dbread_debug dbread_text dbread_pands
# While hitran4 is unstable, following is commented from previous line
# dbread_hitran4
lineread_OBJS_$(d)     := $(lineread_FILES_$(d):%=$(d)/%.o)
lineread_PIC_OBJS_$(d) := $(lineread_FILES_$(d):%=$(d)/%_pic.o)
lineread_DEPS_$(d)     := $(lineread_FILES_$(d):%=$(d)/%.o.d)    \
		          $(lineread_FILES_$(d):%=$(d)/%_pic.o.d)

test_FILES_$(d) := ../transit/transitstd.o test_dbread.c

test_pands_INCLUDE_$(d) := dbread_pands
test_pands_DEPS_$(d)   := $(test_pands_INCLUDE_$(d):%=$(d)/%.o.d)
test_pands_SRC_$(d)   := $(test_pands_INCLUDE_$(d):%=$(d)/%.c)

test_text_INCLUDE_$(d) := dbread_text
test_text_DEPS_$(d)    := $(test_text_INCLUDE_$(d):%=$(d)/%.o.d)
test_text_SRC_$(d)    := $(test_text_INCLUDE_$(d):%=$(d)/%.c)

#gather all binaries together
PIC_OBJS_$(d)    := $(lineread_PIC_OBJS_$(d))
OBJS_$(d)        := $(lineread_OBJS_$(d)) $(test_pands_OBJS_$(d)) $(test_text_OBJS_$(d))
DEPS_$(d)        := $(lineread_DEPS_$(d)) $(test_pands_DEPS_$(d)) $(test_text_DEPS_$(d))

CLEAN := $(CLEAN) $(OBJS_$(d)) $(PIC_OBJS_$(d)) $(DEPS_$(d)) $(local_$(d))
bin_PROGRAMS := $(bin_PROGRAMS) $(bin_PROGRAMS_$(d))
test_PROGRAMS := $(test_PROGRAMS) $(test_PROGRAMS_$(d))
lib_DYNAMIC  := $(lib_DYNAMIC)  $(lib_DYN_$(d))
lib_STATIC   := $(lib_STATIC)   $(lib_STT_$(d))
PIC_OBJS     := $(PIC_OBJS)     $(PIC_OBJS_$(d))
OBJS         := $(OBJS)         $(OBJS_$(d))

#extra dependencies and local variables
$(lib_STT_$(d)): $(OBJS_$(d))
$(lib_DYN_$(d)): $(PIC_OBJS_$(d))
$(test_PROGRAMS_$(d):%=%$(EXEEXT)): $(test_FILES_$(d):%=$(d)/%)
lineread$(EXEEXT):   $(lineread_OBJS_$(d))

LOCAL_INC_$(d) := -I$(d) 
$(OBJS_$(d):.o=.proto): CP_LOCAL :=  $(LOCAL_INC_$(d))
$(local_$(d)): $(d)/Rules.mk
$(local_$(d)): CF_LOCAL := $(LOCAL_INC_$(d)) \
	-DHAVE_INLINE 	-D_BSD_SOURCE #\
#	`pkg-config --cflags gtk+-2.0`
$(local_$(d)): LL_LOCAL := -lm -lpu \
#`pkg-config --libs gtk+-2.0` \
	-lplplotd -lcfitsio -lblas
test_pands: CF_TEST := -DTEST_RUN \
	$(foreach inc, $(test_pands_INCLUDE_$(d):%=$(d)/%.c), -include $(inc))
test_text:  CF_TEST := -DTEST_RUN \
	$(foreach inc, $(test_text_INCLUDE_$(d):%=$(d)/%.c), -include $(inc))

#Standard closure
-include $(DEPS_$(d))

sp   := $(basename $(sp))
d    := $(dirstack_$(sp))
