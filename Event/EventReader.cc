#include <algorithm>   

#include "Event/EventReader.hh"

EventReader::EventReader():
_useWeight(false)
{
}
EventReader::EventReader(bool useWeight):
_useWeight(useWeight)
,_useMassRange(false)
,_massMin(0.)
,_massMax(50.)
{
}
EventReader::~EventReader()
{
}

bool EventReader::isMassrangeParticle(int index){
  bool result=false;

  if (std::find(_particleIdxesMassRange.begin(), _particleIdxesMassRange.end(), index) != _particleIdxesMassRange.end())
    result=true;
    
  return result;
}
