#===============================================================================
# FILE: Makefile
#-------------------  macro definitions  ---------------------------------------
SHELL = /bin/sh
builddir = lib

CPP      = gcc
CPPFLAGS =
CC       = gcc
CFLAGS   = -lm
LDFLAGS  = 


#####
EXE_DIR    := bin/
EXECUTABLE := $(EXE_DIR)react1
SRC_FILES  := $(wildcard src/chemistry/*.c)\
              $(wildcard src/*.c)
OBJ_DIR    := obj/
OBJ_FILES  := $(addprefix $(OBJ_DIR), $(notdir $(SRC_FILES:.c=.o)))
SRC_DIR    := $(dir $(SRC_FILES))
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
	$(CC) $(CPPFLAGS)\
	  $(CFLAGS) -c $< -o $@

# Link the objects to executable
$(EXECUTABLE) : $(OBJ_FILES)
	$(CC) -g -o $@ ${OBJ_FILES} $(CFLAGS) $(LDFLAGS)

# clean source file
.PHONY: clean
clean :
	rm -rf $(OBJ_DIR)*
	rm -rf $(EXECUTABLE)

