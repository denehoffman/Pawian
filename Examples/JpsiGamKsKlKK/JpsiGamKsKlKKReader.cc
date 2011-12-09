#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKReader.hh"

#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"

JpsiGamKsKlKKReader::JpsiGamKsKlKKReader()
{}

JpsiGamKsKlKKReader::JpsiGamKsKlKKReader(const std::vector<std::string>& files, int particles, int skip):
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

JpsiGamKsKlKKReader::~JpsiGamKsKlKKReader()
{}

bool JpsiGamKsKlKKReader::fillAll(EventList& evtList)
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
      for (parts = 0; parts < numParticles; parts++) {
	currentStream >> px >> py >> pz >> e;
        newEvent->addParticle(e,px,py,pz);
      }
      if (!currentStream.fail()) {
	evtList.add(newEvent);
	for (parts = 0; parts < linesToSkip; parts++)
	  currentStream >> px >> py >> pz >> e;
      }
    }
    currentStream.close();
    ++currentFile;
  }
  evtList.rewind();

  return true; // success
}

