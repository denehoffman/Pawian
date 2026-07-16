CXX := c++
TARGET := printAmpToolsIntensity
SRC := main.cc ROOTDataReader.cc Zlm.cc wignerD.cc
OBJ := $(SRC:.cc=.o)

AMPTOOLS_INC := $(AMPTOOLS_HOME)/AmpTools
AMPTOOLS_LIB := $(AMPTOOLS_HOME)/AmpTools/lib

CXXFLAGS := -O2 -std=c++17 -I. -I$(AMPTOOLS_INC) $(shell root-config --cflags)
LDFLAGS := -L$(AMPTOOLS_LIB) -lAmpTools $(shell root-config --libs)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all clean
