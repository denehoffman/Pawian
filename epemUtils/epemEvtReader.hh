// epemEvtReader class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once

#include "Event/EventReader.hh"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <utility>


class EventList;

class epemEvtReader : public EventReader
{
public:
  epemEvtReader();
  epemEvtReader(const std::vector<std::string>& files, int particles, int skip, bool useWeight=false);

  virtual ~epemEvtReader();

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
};


