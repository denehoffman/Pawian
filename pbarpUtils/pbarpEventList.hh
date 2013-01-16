// pbarpEventList class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/EvtDataBaseList.hh"


class EventList;

class pbarpEventList : public EvtDataBaseList {

public:

  // create/copy/destroy:

  ///Constructor 
 pbarpEventList();



  /** Destructor */
  virtual ~pbarpEventList();

  // Getters:
  

protected:


  virtual void read4Vecs(EventList& evtList, std::vector<EvtData*>& theEvtList, double& evtWeightSum, int maxEvts);
  
private:

  std::string getName(std::vector<Particle*>& theVec);
};

