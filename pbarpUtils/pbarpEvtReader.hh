// pbarpEvtReader class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include "Event/EventReader.hh"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <utility>


class EventList;

class pbarpEvtReader : public EventReader
{
public:
  pbarpEvtReader();
  pbarpEvtReader(const std::vector<std::string>& files, int particles, int skip, bool useWeight=false);

  virtual ~pbarpEvtReader();

  virtual bool fillAll(EventList& evtList);
  
private:
  std::vector<std::string> fileNames;
  std::vector<std::string>::const_iterator currentFile;
  std::ifstream currentStream;
  int numParticles;
  int linesToSkip;
};


