#pragma once

#include "Event/EventReader.hh"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <utility>


class EventList;

class D0ToKsPipPimReader : public EventReader
{
public:
  D0ToKsPipPimReader();
  D0ToKsPipPimReader(const std::vector<std::string>& files, int particles, int skip, bool useWeight=false);

  virtual ~D0ToKsPipPimReader();

  virtual bool fillAll(EventList& evtList);

  
private:
  std::vector<std::string> fileNames;
  std::vector<std::string>::const_iterator currentFile;
  std::ifstream currentStream;
  int numParticles;
  int linesToSkip;
  
};


