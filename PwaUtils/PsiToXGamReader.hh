#ifndef PsiToXGamReader_HH
#define PsiToXGamReader_HH

#include "Event/EventReader.hh"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <utility>


class EventList;

class PsiToXGamReader : public EventReader
{
public:
  PsiToXGamReader();
  PsiToXGamReader(const std::vector<std::string>& files, int particles, int skip, bool useWeight=false);

  virtual ~PsiToXGamReader();

  virtual bool fillAll(EventList& evtList);
  bool fillMassRange(EventList& evtList, std::pair<double,double> massRange){
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
