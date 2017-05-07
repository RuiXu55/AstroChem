#===============================================================================
# FILE: Makefile
#-------------------  macro definitions  ---------------------------------------
SHELL = /bin/sh
builddir = lib

LIBTOOL  = $(SHELL) $(builddir)/libtool
CPP      = gcc -E
CPPFLAGS =
CC       = gcc
CFLAGS   = 
LDFLAGS  = 

SUNDIALS_INCS = -I$(builddir)/include
SUNDIALS_LIBS = $(builddir)/src/cvode/libsundials_cvode.la   \
                $(builddir)/src/nvec_ser/libsundials_nvecserial.la

#####
EXE_DIR    := bin/
EXECUTABLE := $(EXE_DIR)react
SRC_FILES  := $(wildcard src/chemistry/*.c)\
              $(wildcard src/*.c)
OBJ_DIR    := obj/
OBJ_FILES  := $(addprefix $(OBJ_DIR), $(notdir $(SRC_FILES:.c=.o)))
SRC_DIR    := $(dir $(SRC_FILES) $(PROB_FILES))
VPATH      := $(SRC_DIR)


.PHONY : all dirs clean

all: dirs $(EXECUTABLE)

dirs : $(EXE_DIR) $(OBJ_DIR)

$(EXE_DIR):
	mkdir -p $(EXE_DIR)
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Create Objects from source files
$(OBJ_DIR)%.o : %.c
	$(CC) $(CPPFLAGS) $(SUNDIALS_INCS) \
	  $(CFLAGS) -c $< -o $@

# Link the objects to executable
$(EXECUTABLE) : $(OBJ_FILES)
	$(LIBTOOL) --mode=link $(CC) ${OPT} -g -o $@ \
          ${OBJ_FILES} $(CFLAGS) $(LDFLAGS) $(SUNDIALS_LIBS) 

# clean source file
.PHONY: clean
clean :
	rm -rf $(OBJ_DIR)*
	rm -rf $(EXECUTABLE)


