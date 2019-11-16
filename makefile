CC = gcc
CXX = g++
CFLAGS  = 	-I. -O3 -Wall -Wextra -pedantic -ansi -c -Wno-write-strings

all:
	g++ -march=native -I./include -I ./include/types/ src/* -o xavier

demo: demo.cpp
	$(CXX) -march=native -std=c++14 -O3 -fpermissive -o demo demo.cpp -DDEBUG

clean:
	rm -f *.o
	rm -f demo

