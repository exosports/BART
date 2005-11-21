#Adapted from recipes by Emile van Bergen, 2002

#standard things
sp             := $(sp).x
dirstack_$(sp) := $(d)

#subdirectories in random order
#d       := $(d)/name
#include   $(d)/Rules.mk

#Standard closure
#-include $(DEPS_$(d))

sp   := $(basename $(sp))
d    := $(dirstack_$(sp))
