LZ4INC = -I../lz4/lib
ROOTLIBS = $(shell root-config --libs)
CXXFLAGS = $(shell root-config --cflags) -Ilibcpp $(LZ4INC)
LIBFLAG = -c -O3 -fPIC -m64 -fmessage-length=0 -g $(LZ4INC) -D__LZ4__ -D__LIBZ__
LIB = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
CXX = g++

LZ4 = lz4/lib/lz4.o
.PHONY: clean

all: $(LZ4) $(LIB)

$(LIB): %.o: %.cpp
	$(CXX) -std=c++11 $(LIBFLAG) -o $@ $<

$(LZ4):
	make lib -C ../lz4

clean:
	-rm -f $(LIB)
	
purge: clean
	make clean -C ../lz4
