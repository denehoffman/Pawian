#pragma once

#include <string>
#include <vector>

class EventList;

class EventReader
{
public:
  EventReader();
  EventReader(bool useWeight);
  virtual ~EventReader();

  virtual bool fillAll(EventList& evtList) = 0;

protected:
 bool _useWeight; 
};
