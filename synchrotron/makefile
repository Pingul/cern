CC = clang++
CPLUS_FLAGS +=-Wall -pedantic
INCLUDE_FLAGS = -I/opt/local/include -I../common -Isource
LIBRARY_PATHS = -L/opt/local/lib 
LIBRARIES = -ltbb
H_DEPS = $(wildcard source/*.hh) ../common/timer.hh
CPP = source/main.cc ../common/timer.cc
#H_DEPS = source/2d-synchrotron.hh source/settings.hh source/accelerator.hh source/hamiltonian.hh ../common/timer.hh

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
# This is meant for simulations on lxplus
	LIBRARY_PATHS += -L$(TBB_LIBRARY_RELEASE)
	INCLUDE_FLAGS += -I$(TBB_INCLUDE)
	CC = g++
# Not linking statically to libstdc++ crashes for some reason
	CPLUS_FLAGS += -static-libstdc++
# IO_TO_DIR makes the 2d-synchrotron look for input and put output in the same directly as the binary
# was called
	CPLUS_FLAGS += -DIO_TO_SAME_DIR
endif

ifeq ($(UNAME), Darwin)
	CPLUS_FLAGS += -std=c++11 -stdlib=libc++ 
endif


all: 2dsynch

2dsynch: ${H_DEPS} source/main.cc ../common/timer.cc
	$(CC) $(CPLUS_FLAGS) -o 2dsynch $(INCLUDE_FLAGS) $(LIBRARY_PATHS) $(LIBRARIES) $(CPP)

misc: source/misc.cc $(H_DEPS)
	$(CC) $(CPLUS_FLAGS) -o misc $(INCLUDE_FLAGS) $(LIBRARY_PATHS) $(LIBRARIES) source/misc.cc

tool: source/exportTool.cc $(H_DEPS)
	$(CC) $(CPLUS_FLAGS) -o tool $(INCLUDE_FLAGS) $(LIBRARY_PATHS) $(LIBRARIES) source/exportTool.cc

clean :
	rm -f 2dsynch misc

