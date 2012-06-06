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
