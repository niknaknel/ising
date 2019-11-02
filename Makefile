# sample Makefile.
# It compiles every .c files in the src/ directory to object files in obj/ directory, and build the ./my_executable my_executable.

# using gcc :
COMPILER ?= $(GCC_PATH)gcc

DEBUG ?= -DDEBUG
FLAGS ?= $(DEBUG) -fopenmp -O2 -Wall -Wno-variadic-macros -pedantic -g $(GCC_SUPPFLAGS)

LDFLAGS ?= -g
LDLIBS = -lm

EXECUTABLE = ising

SRCS=$(wildcard src/*.c)
OBJS=$(SRCS:src/%.c=obj/%.o)

all: release

release: $(OBJS)
	$(COMPILER) $(LDFLAGS) -fopenmp -o $(EXECUTABLE) $(OBJS) $(LDLIBS) 

obj/%.o: src/%.c
	@mkdir -p $(@D)
	$(COMPILER) $(FLAGS) -fopenmp -o $@ -c $<

clean:
	rm -f obj/*
	rm ${EXECUTABLE} 

cleandata:
	rm cachegrind.out.*

dist-clean: clean
	rm -f $(EXECUTABLE) *~ .depend *.zip
	
#automatically handle include dependencies
#depend: .depend
#
#.depend: $(SRCS)
#	rm -f ./.depend
#	@$(foreach SRC, $(SRCS), $(COMPILER) $(FLAGS) -MT $(SRC:src/%.cpp=obj/%.o) -MM $(SRC) >> .depend;)
#
#include .depend	
