CC= 	gcc
CXX= 	g++
INCLUDE= 	-Itools/libgaba -I../include -Itools -Itools/seqan/include
LIBS= 		-Ltools/libgaba -lgaba
CFLAGS = 	-I. -O3 -std=c++14 -march=native -Wall -Wextra -pedantic -Wno-write-strings -fpermissive -fopenmp

# the build target executable:
TARGET = test-random test-overlap test-runtime
all: $(TARGET)

test-random: test-score-random.cpp tools/ksw2
	$(CXX) -march=native -std=c++14 $(INCLUDE) test-score-random.cpp ../lib/xavierlib.a $(CFLAGS) -o test-random-exe $(LIBS)

test-overlap: test-score-overlap.cpp tools/ksw2
	$(CXX) -march=native -std=c++14 $(INCLUDE) test-score-overlap.cpp ../lib/xavierlib.a $(CFLAGS) -o test-overlap-exe $(LIBS)

test-runtime: test-runtime-overlap.cpp tools/ksw2
	$(CXX) -march=native -std=c++14 $(INCLUDE) test-runtime-overlap.cpp ../lib/xavierlib.a $(CFLAGS) -o test-runtime-exe $(LIBS)

# test-tn: test-tn.cpp tools/ksw2
# 	$(CXX) -march=native -std=c++14 $(INCLUDE) test-tn.cpp ../lib/xavierlib.a $(CFLAGS) -o test-tn-exe $(LIBS)

clean:
	rm -rf *.o test*exe
