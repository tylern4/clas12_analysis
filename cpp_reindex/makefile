ROOTLIBS	= $(shell root-config --libs)
CXX = g++
CXXFLAGS = -O3 -march=native -fPIC -w -g $(shell root-config --cflags)
TARGET =	    clas12_ana

MAINS = clas12branches.cpp main.cpp histogram.cpp reaction.cpp deltat.cpp physics.cpp
MAIN = $(MAINS:%.cpp=%.o)

.PHONY: all clean

all:	$(TARGET)

$(MAIN): %.o : %.cpp
		$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(MAIN)
	$(CXX) $(MAIN) $(CXXFLAGS) $(ROOTLIBS) -o $(TARGET)


clean:
	-rm -f $(TARGET) $(MAIN)
