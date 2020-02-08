CC= 	gcc
CXX= 	g++
INCLUDE= 	-Itools/libgaba -I../include -Itools -Itools/seqan/include
LIBS= 		-Ltools/libgaba -lgaba
CFLAGS = 	-I. -O3 -std=c++14 -march=native -Wall -Wextra -pedantic -Wno-write-strings -fpermissive

test-pc: test-pc.cpp tools/ksw2
	$(CXX) -march=native -std=c++14 $(INCLUDE) test-pc.cpp ../lib/xavierlib.a $(CFLAGS) -o test-pc-exe $(LIBS)

test-tp: test-tp.cpp tools/ksw2
	$(CXX) -march=native -std=c++14 $(INCLUDE) test-tp.cpp ../lib/xavierlib.a $(CFLAGS) -o test-tp-exe $(LIBS)

test-tn: test-tn.cpp tools/ksw2
	$(CXX) -march=native -std=c++14 $(INCLUDE) test-tn.cpp ../lib/xavierlib.a $(CFLAGS) -o test-tn-exe $(LIBS)

clean:
	rm -rf *.o test*exe
