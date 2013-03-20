#include "Event/EventReader.hh"

EventReader::EventReader():
_useWeight(false)
{
}
EventReader::EventReader(bool useWeight):
_useWeight(useWeight)
,_useMassRange(false)
,_massRange(std::make_pair(0.,10.) )
{
}
EventReader::~EventReader()
{
}
