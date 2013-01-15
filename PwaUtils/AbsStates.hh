// AbsStates class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>

class AbsStates {

public:
  AbsStates();
  virtual ~AbsStates();
  virtual void print(std::ostream& os) const;

 protected:
  virtual bool calcStates() = 0;
};
