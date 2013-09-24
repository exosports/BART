#Adapted from recipes by Emile van Bergen, 2002

#standard things
sp             := $(sp).x
dirstack_$(sp) := $(d)

#subdirectories in random order
d       := $(d)/pu
include   $(d)/Rules.mk

FILES_PROGRAM_$(d)   := pu

header_PROGRAM_$(d) := $(FILES_PROGRAM_$(d):%=$(d)/%.h)

#all together
header_PROGRAM := $(header_PROGRAM) $(header_PROGRAM_$(d))

#Standard closure
#-include $(DEPS_$(d))

sp   := $(basename $(sp))
d    := $(dirstack_$(sp))
