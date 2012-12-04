#include "Examples/D0ToKsPipPim/D0ToKsPipPimReader.hh"

#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"

D0ToKsPipPimReader::D0ToKsPipPimReader()
{}

D0ToKsPipPimReader::D0ToKsPipPimReader(const std::vector<std::string>& files, int particles, int skip, bool useWeight):
  EventReader(useWeight),
  numParticles(particles),
  linesToSkip(skip)
{
  if (0 == files.size()) {
    Alert << "empty list of event files" ;  // << endmsg;
    exit(1);
  }
  std::vector<std::string>::const_iterator iter = files.begin();
  for (; iter != files.end(); ++iter)
    fileNames.push_back(*iter);
  currentFile = fileNames.begin();
}

D0ToKsPipPimReader::~D0ToKsPipPimReader()
{}

bool D0ToKsPipPimReader::fillAll(EventList& evtList)
{

  while (currentFile != fileNames.end()) {
    currentStream.open(currentFile->c_str());
    if (!currentStream) {
      Alert << "can not open " << *currentFile ;  // << endmsg;
      exit(1);
    }

    while (!currentStream.eof()) {
      double e,px,py,pz;
      Event* newEvent = new Event();
      int parts;

      if(_useWeight)
      {
         double weight;
         currentStream >> weight;
         newEvent->addWeight(weight);
      }      

      Vector4<double> fvX(0,0,0,0); //X four-vector
      for (parts = 0; parts < numParticles; parts++) {
	currentStream >> e >> px >> py >> pz;
        newEvent->addParticle(e,px,py,pz);
	Vector4<double> tmp = newEvent->p4(parts);
	if(parts>0) fvX= fvX+tmp;
      }
      
      if (!currentStream.fail()) {
       	evtList.add(newEvent);
       	for (parts = 0; parts < linesToSkip; parts++)
       	  currentStream >> e >> px >> py >> pz;
      }
    }
    //remove last event
    evtList.removeEvents(evtList.size(), evtList.size());
    currentStream.close();
    ++currentFile;
  }
  evtList.rewind();

  return true; // success
}

