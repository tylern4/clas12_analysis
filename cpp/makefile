UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
    FOPENMP = -fopenmp -lgfortran
endif

ROOTLIBS	= $(shell root-config --libs)
CXX = g++
CXXFLAGS = -O3 -march=native -fopenmp -fPIC -w -g $(FOPENMP) $(shell root-config --cflags)
TARGET =	    test

MAINS = main.cpp
MAIN = $(MAINS:%.cpp=%.o)

.PHONY: all clean

all:	$(TARGET)

$(MAIN): %.o : %.cpp
		$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(MAIN)
	$(CXX) $(MAIN) $(CXXFLAGS) $(ROOTLIBS) -o $(TARGET)


clean:
	-rm -f $(TARGET) $(MAIN)
