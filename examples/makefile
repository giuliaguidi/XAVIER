CFLAGS = -std=c++14 -march=native -Wall -Wextra -pedantic -Wno-write-strings -fsanitize=address -O1 -fno-omit-frame-pointer
build_demo:
	g++ demo.cpp ../lib/xavierlib.a -I../include $(CFLAGS) -o demo