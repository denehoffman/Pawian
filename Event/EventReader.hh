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

  virtual bool fillMassRange(EventList& evtList, std::pair<double,double> massRange){
    _useMassRange=true;
    _massRange = massRange;
    return fillAll(evtList);
  };

protected:
  bool _useWeight;
  bool _useMassRange;
  std::pair<double,double> _massRange; 
};
