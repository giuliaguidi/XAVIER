CC= gcc
CXX= g++
INCLUDE= -Itools/libgaba -I../include -Itools -Itools/seqan/include -IComplete-Striped-Smith-Waterman-Library/src
LIBS= -Ltools/libgaba -lgaba -lparasail
CFLAGS = -I. -O3 -std=c++14 -march=native -Wall -Wextra -pedantic -Wno-write-strings -fpermissive

build_benchmark: tools/ksw2
	$(CXX) -march=native -std=c++14 $(INCLUDE) benchmark.cpp ../lib/xavierlib.a  $(CFLAGS) -o test $(LIBS)