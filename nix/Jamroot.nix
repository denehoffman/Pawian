import os ;
import testing ;

path-constant TOP : . ;

local MINUIT2_STANDALONE = [ os.environ MINUIT2_STANDALONE ] ;
local MINUIT2_LIBRARY = [ os.environ MINUIT2_LIBRARY ] ;
local ROOTSYS = [ os.environ ROOTSYS ] ;
local BOOSTINCLUDE = [ os.environ BOOSTINCLUDE ] ;
local BOOSTLIBPATH = [ os.environ BOOSTLIBPATH ] ;
local rlibs = [ SHELL "$(ROOTSYS)/bin/root-config --libs" ] ;
ROOTLIBS = [ MATCH "(.*)[\n]" : $(rlibs) ] ;

BOOSTLIBS =
  -lboost_chrono
  -lboost_date_time
  -lboost_filesystem
  -lboost_program_options
  -lboost_regex
  -lboost_serialization
  -lboost_system
  -lboost_timer
  -lboost_unit_test_framework
  ;

lib boost_test : : <name>boost_unit_test_framework ;

project :
  requirements <include>./
  <include>$(TOP)
  <include>$(BOOSTINCLUDE)
  <include>$(MINUIT2_STANDALONE)/include/Minuit2
  <include>$(ROOTSYS)/include
  <cxxflags>--std=c++17
  <cxxflags>-ftemplate-depth=256
  <cxxflags>-DBOOST_BIND_GLOBAL_PLACEHOLDERS
  <link>shared
  <linkflags>-L$(TOP)/lib
  <linkflags>$(ROOTLIBS)
  <linkflags>-L$(BOOSTLIBPATH)
  <linkflags>$(BOOSTLIBS)
  <linkflags>-lgsl
  <linkflags>-lgslcblas
  <linkflags>-llog4cpp
  <cxxflags>-fPIC
  <cxxflags>-Wall
  ;

actions rootlibs
{
  $(ROOTSYS)/bin/root-config --libs
}

lib Minuit2 : : <file>$(MINUIT2_LIBRARY) : : ;

build-project ErrLogger ;
build-project qft++ ;
build-project qft++Extension ;
build-project ConfigParser ;
build-project Particle ;
build-project Utils ;
build-project Event ;
build-project FitParams ;
build-project PwaDynamics ;
build-project PwaUtils ;
build-project MinFunctions ;
build-project AppUtils ;
build-project pbarpUtils ;
build-project epemUtils ;
build-project resUtils ;
build-project ggUtils ;
build-project gammapUtils ;
build-project pipiScatteringUtils ;
build-project gslUtils ;
build-project qaErrorExtract ;
build-project PspGen ;
build-project Examples ;
build-project KMatrixExtract ;
build-project PwaApps ;

alias gammap-install : Examples/gammap//install ;
