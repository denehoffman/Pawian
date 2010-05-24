#ifndef EVENTREADER_HH
#define EVENTREADER_HH

#include <string>
#include <vector>

class EventList;

class EventReader
{
public:
  EventReader();
  virtual ~EventReader();

  virtual bool fillAll(EventList& evtList) = 0;
 
};

#endif
