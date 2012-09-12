#include "Event/EventReader.hh"

EventReader::EventReader():
_useWeight(false)
{
}
EventReader::EventReader(bool useWeight):
_useWeight(useWeight)
{
}
EventReader::~EventReader()
{
}
