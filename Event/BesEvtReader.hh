#ifndef BesEvtReader_HH
#define BesEvtReader_HH

#include "Event/EventReader.hh"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>


class EventList;

class BesEvtReader : public EventReader
{
public:
  BesEvtReader();
  BesEvtReader(const std::vector<std::string>& files, int particles, int skip);

  virtual ~BesEvtReader();

  virtual bool fillAll(EventList& evtList);
 
private:
  std::vector<std::string> fileNames;
  std::vector<std::string>::const_iterator currentFile;
  std::ifstream currentStream;
  int numParticles;
  int linesToSkip;
};

#endif
