ROOTLIBS = $(shell root-config --libs) $(shell pkg-config --libs hipocpp)
CXXFLAGS = $(shell root-config --cflags) $(shell pkg-config --cflags hipocpp)
LIBFLAG = -c -O3 -fPIC -m64 -fmessage-length=0 -g
LIB = $(patsubst %.cpp,%.o,$(wildcard ../libcpp/*.cpp))
CXX = g++
MAINS = main.cpp
TARGET = clas12_ana
MAIN = $(MAINS:%.cpp=%.o)

UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	CXXFLAGS += -fopenmp
endif

.PHONY: clean

all: $(TARGET)

$(MAIN): %.o : %.cpp
		$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(MAIN)
	$(CXX) $(MAIN) $(CXXFLAGS) $(ROOTLIBS) -o $(TARGET)

clean:
	-rm -f $(TARGET) $(MAIN)
purge: clean
	-rm -f $(LIB)
	make clean -C ../lz4
