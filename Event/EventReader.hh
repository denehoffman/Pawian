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

  virtual bool fillMassRange(EventList& evtList, double massMin, double massMax, std::vector<unsigned int>& particleIndices){
    _useMassRange=true;
    _massMin = massMin;
    _massMax = massMax;
    _particleIdxesMassRange=particleIndices;
    return fillAll(evtList);
  };

protected:
  bool _useWeight;
  bool _useMassRange;
  double _massMin;
  double _massMax;
  std::vector<unsigned int> _particleIdxesMassRange;

  bool isMassrangeParticle(int index); 
};
