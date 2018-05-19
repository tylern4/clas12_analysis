LZ4INC = -I../lz4/lib
ROOTLIBS = $(shell root-config --libs)
CXXFLAGS = $(shell root-config --cflags) -I../libcpp $(LZ4INC)
LIBFLAG = -c -O3 -fPIC -m64 -fmessage-length=0 -g $(LZ4INC) -D__LZ4__ -D__LIBZ__
LIB = $(patsubst %.cpp,%.o,$(wildcard ../libcpp/*.cpp))
CXX = g++
MAINS = deltat.cpp histogram.cpp physics.cpp main.cpp
TARGET = clas12_ana
MAIN = $(MAINS:%.cpp=%.o)

LZ4 = ../lz4/lib/lz4.o
.PHONY: clean

all: $(TARGET)

$(LIB): %.o: %.cpp
	$(CXX) -std=c++11 $(LIBFLAG) -o $@ $<

$(LZ4):
	make lib -C ../lz4

$(MAIN): %.o : %.cpp
		$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(MAIN) $(LZ4) $(LIB)
	$(CXX) $(MAIN) $(LIB) $(LZ4) $(CXXFLAGS) $(ROOTLIBS) -o $(TARGET)

clean:
	-rm -f $(TARGET) $(MAIN)
purge: clean
	-rm -f $(LIB)
	make clean -C ../lz4
