#sources
lineread_SOURCES = ../transit/transitstd.c dbread_pands.c lineread.c
lineread_OBJECTS = $(lineread_SOURCES:.c=.o)
bin_PROGRAMS     = lineread$(EXEEXT)

PROGRAMS         = $(bin_PROGRAMS)
OBJECTS          = $(lineread_OBJECTS)
SOURCES          = $(lineread_SOURCES)

INCLUDES = -I. -I$(tobasedir)/include -I../transit

vpath %.c atmosphere

all: $(bin_PROGRAMS)

lineread$(EXEEXT): $(lineread_OBJECTS)
dep:               $(SOURCES)
