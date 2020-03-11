CC = gcc
CXX = g++
CFLAGS = -std=c++17 -mavx2 -Wall -Wextra -pedantic -Wno-write-strings -O3

ODIR = obj
LDIR = lib
OUT = $(LDIR)/libxavier.a
SDIR = src
INC = -Iinclude

_OBJS = score.o \
		seed.o \
		aligner.o \
		trace.o \
		vectors.o \
		xavier.o

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

$(ODIR)/%.o: $(SDIR)/%.cpp
	mkdir -p $(ODIR)
	$(CXX) $(CFLAGS) $(INC) -c -o $@ $<

$(OUT): $(OBJS)
	mkdir -p $(LDIR)
	ar rvs $(OUT) $^
	ranlib $(OUT)

.PHONY: clean

clean:
	rm -rf $(ODIR) $(LDIR) $(OUT)