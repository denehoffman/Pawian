#ifndef Psi2SToKpKmPiGamReader_HH
#define Psi2SToKpKmPiGamReader_HH

#include "Event/EventReader.hh"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>


class EventList;

class Psi2SToKpKmPiGamReader : public EventReader
{
public:
  Psi2SToKpKmPiGamReader();
  Psi2SToKpKmPiGamReader(const std::vector<std::string>& files, int particles, int skip);

  virtual ~Psi2SToKpKmPiGamReader();

  virtual bool fillAll(EventList& evtList);
 
private:
  std::vector<std::string> fileNames;
  std::vector<std::string>::const_iterator currentFile;
  std::ifstream currentStream;
  int numParticles;
  int linesToSkip;
};

#endif
