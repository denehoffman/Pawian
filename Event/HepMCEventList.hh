#pragma once

#include <vector>
#include <string>
#include "HepMC/GenEvent.h"

class HepMCEventList
{
public:
  HepMCEventList();  
  HepMCEventList(const std::vector<std::string>& files);
  HepMCEventList(const std::string& file);
  ~HepMCEventList();

  void add(HepMC::GenEvent*);
  void removeEvents(unsigned int nBegin, unsigned int nEnd);

  HepMC::GenEvent* nextEvent();
  void rewind();
  int size();

private:
  bool fillFromFiles(const std::vector<std::string>& files);
  bool fillFromFile(const std::string& file);

  std::vector<HepMC::GenEvent*> eventList;
  std::vector<HepMC::GenEvent*>::const_iterator currentEvent;
  
};
