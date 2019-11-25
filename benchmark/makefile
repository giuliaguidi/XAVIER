CC= gcc
CXX= g++
INCLUDE= -Itools/libgaba -I../include -Itools -I/global/cscratch1/sd/gguidi/xavier/benchmark/tools/parasail/include -Itools/seqan/include
LIBS= -Ltools/libgaba -L/global/cscratch1/sd/gguidi/xavier/benchmark/tools/parasail/lib -lgaba -lparasail #-ledlib
CFLAGS = -I. -O3 -std=c++14 -march=native -Wall -Wextra -pedantic -Wno-write-strings -fpermissive

build_benchmark: tools/ksw2
	$(CXX) -march=native -std=c++14 $(INCLUDE) benchmark.cpp ../lib/xavierlib.a  $(CFLAGS) -o test $(LIBS)