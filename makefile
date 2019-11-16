CC = gcc
CXX = g++
CFLAGS = -march=native
#CFLAGS = -O3 -std=c++14  -march=native -Wall -Wextra -pedantic -ansi -Wno-write-strings


ODIR = obj
LDIR = lib
OUT = $(LDIR)/xavierlib.a
SDIR = src
INC = -Iinclude -Iinclude/types

_OBJS = score.o \
		seed.o \
		state.o \
		trace.o \
		vectors.o \
		xavier.o

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))


$(ODIR)/%.o: $(SDIR)/%.cpp
	mkdir -p $(ODIR)
	$(CXX) $< -c -o $@ $(INC) $(CFLAGS)

$(OUT): $(OBJS)
	mkdir -p $(LDIR)
	ar rvs $(OUT) $^

.PHONY: clean

clean:
	rm -rf $(ODIR) $(LDIR) $(OUT)