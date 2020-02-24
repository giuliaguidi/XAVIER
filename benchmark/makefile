CC = 	gcc
CXX= 	g++
INCLUDE= 	-Itools/libgaba -I../include -Itools -Itools/seqan/include
LIBS= 		-Ltools/libgaba -lgaba
CFLAGS = 	-I. -O3 -mavx2 -Wall -Wextra -pedantic -Wno-write-strings -fopenmp

# the build target executable:
TARGET = test-score test-runtime gen-seq
all: $(TARGET)

gen-seq: generateSeq.cpp 
	$(CXX) -mavx2 -std=c++14 -fpermissive $(INCLUDE) generateSeq.cpp ../lib/xavierlib.a $(CFLAGS) -o gen-seq-exe

test-score: testScore.cpp tools/ksw2  
	$(CXX) -mavx2 -std=c++14 -fpermissive $(INCLUDE) testScore.cpp ../lib/xavierlib.a $(CFLAGS) -o test-score-exe $(LIBS)

test-runtime: testRuntime.cpp tools/ksw2  
	$(CXX) -mavx2 -std=c++14 -fpermissive $(INCLUDE) testRuntime.cpp ../lib/xavierlib.a $(CFLAGS) -o test-runtime-exe $(LIBS)

clean:
	rm -rf *.o *-exe
