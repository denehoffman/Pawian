CXX ?= c++
TARGETS := fitZlm amptoolsReference
COMMON_OBJ := ROOTDataReader.o Zlm.o wignerD.o
FIT_OBJ := main.o $(COMMON_OBJ)
REFERENCE_OBJ := reference_main.o $(COMMON_OBJ)

AMPTOOLS_INC := $(AMPTOOLS_HOME)/AmpTools
AMPTOOLS_LIB := $(AMPTOOLS_HOME)/AmpTools/lib
AMPTOOLS_REFERENCE_REVISION := $(shell python3 -c 'import json; print(json.load(open("reference/compatibility-artifacts.json"))["reference"]["revision"])')

CXXFLAGS := -O2 -std=c++17 -I. -I$(AMPTOOLS_INC) $(shell root-config --cflags)
CXXFLAGS += -DAMPTOOLS_REFERENCE_REVISION=\"$(AMPTOOLS_REFERENCE_REVISION)\"
CXXFLAGS += -DPAWIAN_REVISION=\"$(shell git -C ../../.. describe --always --dirty)\"
LDFLAGS := -L$(AMPTOOLS_LIB) -lAmpTools $(shell root-config --libs)

all: $(TARGETS)

fitZlm: $(FIT_OBJ)
	$(CXX) $(FIT_OBJ) -o $@ $(LDFLAGS)

amptoolsReference: $(REFERENCE_OBJ)
	$(CXX) $(REFERENCE_OBJ) -o $@ $(LDFLAGS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f main.o reference_main.o $(COMMON_OBJ) $(TARGETS)

.PHONY: all clean
