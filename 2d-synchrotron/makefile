CC = clang++
CPLUS_FLAGS +=-Wall -pedantic
INCLUDE_FLAGS = -I/opt/local/include -I../common
LIBRARY_PATHS = -L/opt/local/lib 
LIBRARIES = -ltbb

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
	CPLUS_FLAGS += -DIO_TO_DIR
endif

ifeq ($(UNAME), Darwin)
	CPLUS_FLAGS += -std=c++11 -stdlib=libc++ 
endif

all: main

main: main.cc 2d-synchrotron.hh ../common/timer.cc
	$(CC) $(CPLUS_FLAGS) -o main $(INCLUDE_FLAGS) $(LIBRARY_PATHS) $(LIBRARIES) main.cc ../common/timer.cc

clean :
	rm -f main

