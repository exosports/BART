#sources
transit_SOURCES  = transitstd.c readlineinfo.c makesample.c \
		   extinction.c transit.c tau.c idxrefraction.c \
		   argum.c slantpath.c geometry.c observable.c \
		   at_file.c readatm.c at_onept.c cia.c
transit_OBJECTS  = $(transit_SOURCES:.c=.o)
bin_PROGRAMS     = transit$(EXEEXT)

OBJECTS          = $(transit_OBJECTS)
SOURCES          = $(transit_SOURCES)

vpath %.c atmosphere

all:	$(bin_PROGRAMS)

transit$(EXEEXT):  $(transit_OBJECTS)
dep:               $(transit_SOURCES)


