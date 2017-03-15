#tbbinc=/afs/cern.ch/user/s/swretbor/test/tbb/include
#tbblib=/afs/cern.ch/user/s/swretbor/test/tbb/build/linux_intel64_gcc_cc4.8.4_libc2.12_kernel2.6.32_release
UNAME := $(shell uname)

$(warning NAME: $(UNAME))

CC = clang++
CPLUS_FLAGS +=-Wall -pedantic
INCLUDE_FLAGS = -I/opt/local/include -I../common
LIBRARY_PATHS = -L/opt/local/lib 
LIBRARIES = -ltbb

ifeq ($(UNAME), Linux)
	LIBRARY_PATHS += -L$(TBB_LIBRARY_RELEASE)
	INCLUDE_FLAGS += -I$(TBB_INCLUDE)
	CC = g++
	CPLUS_FLAGS += -static-libstdc++
endif

all: main

main: main.cc 2d-synchrotron.hh ../common/timer.cc
	$(CC) $(CPLUS_FLAGS) -o main $(INCLUDE_FLAGS) $(LIBRARY_PATHS) $(LIBRARIES) main.cc ../common/timer.cc

clean :
	rm -f main

