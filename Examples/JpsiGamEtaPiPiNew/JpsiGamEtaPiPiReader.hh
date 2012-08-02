#ifndef JpsiGamEtaPiPiReader_HH
#define JpsiGamEtaPiPiReader_HH

#include "Event/EventReader.hh"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <utility>


class EventList;

class JpsiGamEtaPiPiReader : public EventReader
{
public:
  JpsiGamEtaPiPiReader();
  JpsiGamEtaPiPiReader(const std::vector<std::string>& files, int particles, int skip);

  virtual ~JpsiGamEtaPiPiReader();

  virtual bool fillAll(EventList& evtList);
  virtual bool fillMassRange(EventList& evtList, std::pair<double,double> massRange){
    _useMassRange=true;
    _massRange = massRange;
    return fillAll(evtList);
  };
  
private:
  std::vector<std::string> fileNames;
  std::vector<std::string>::const_iterator currentFile;
  std::ifstream currentStream;
  int numParticles;
  int linesToSkip;
  bool _useMassRange;
  std::pair<double,double> _massRange;
  
};

#endif
