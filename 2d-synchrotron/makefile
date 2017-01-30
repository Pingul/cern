all: main

main: main.cc 2d-synchrotron.hh ../common/timer.cc
	clang++ -std=c++11 -stdlib=libc++ -Wall -pedantic -o main -I/opt/local/include -L/opt/local/lib -ltbb -I../common main.cc ../common/timer.cc


clean :
	rm -f main

