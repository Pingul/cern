all: main

main: main.cc 
	clang++ -std=c++11 -stdlib=libc++ -Wall -pedantic -o main -I/opt/local/include -L/opt/local/lib -ltbb -I../common main.cc


clean :
	rm -f main

