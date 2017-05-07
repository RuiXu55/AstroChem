#===============================================================================
# FILE: Makefile
#
#-------------------  object files  --------------------------------------------
CORE_OBJ = ath_array.o \
	   ath_log.o \
           main.o \
           par.o \
           utils.o

CHEM_OBJ = chemistry/reprocess_list.o\
           chemistry/init_chemistry.o\
           chemistry/extract_reactions.o\
           chemistry/init_species.o\

ALL_OBJ = ${CORE_OBJ} ${CHEM_OBJ}

#-------------------  macro definitions  ---------------------------------------

BIN = ${DIRS}/extract
DIRS = bin
SRC = $(ALL_OBJ:.o=.c)

include Makeoptions

#--------------------  implicit rules  -----------------------------------------

.c.o:
	${CC} ${CFLAGS} -c $<

#---------------------  targets  -----------------------------------------------

all:	cmpchem compile

cmpchem:
	(cd chemistry; $(MAKE) compile)

.PHONY: compile
compile: ${BIN}

${BIN}: ${ALL_OBJ}
	${LDR} $(OPT) -o ${EXEDIR}extract ${ALL_OBJ} ${LIB}

help:
	@echo This is the Athena Makefile
	@echo Type 'make compile' to generate Athena executable
	@echo Type 'make clean'   to remove '*.o' files
	@echo OBJ=$(OBJ)

.PHONY: clean
clean:
	(cd chemistry; $(MAKE) clean)
	(rm -f *.o Makedepend $(BIN))

depend: Makedepend

Makedepend:
	${CC} -MM ${CFLAGS} ${SRC} > Makedepend

# this forces a rebuild of all objects, if code re-configure'd
#include Makedepend
