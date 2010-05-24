#include "Event/CBElsaReader.hh"

#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "ErrLogger/ErrLineLog.hh"

CBElsaReader::CBElsaReader()
{}

CBElsaReader::CBElsaReader(std::vector<std::string>& files, int particles, int skip):
  numParticles(particles),
  linesToSkip(skip)
{
  if (0 == files.size())
    ErrMsg(fatal) << "empty list of event files" << endmsg;
  std::vector<std::string>::iterator iter = files.begin();
  for (; iter != files.end(); ++iter)
    fileNames.push_back(*iter);
  currentFile = fileNames.begin();
}

CBElsaReader::~CBElsaReader()
{}

bool CBElsaReader::fillAll(EventList& evtList)
{
  
  while (currentFile != fileNames.end()) {
    currentStream.open(currentFile->c_str());
    if (!currentStream)
      ErrMsg(fatal) << "can not open " << *currentFile << endmsg;

    while (!currentStream.eof()) {
      double e,px,py,pz;
      Event* newEvent = new Event();
      int parts;
      for (parts = 0; parts < numParticles; parts++) {
	currentStream >> e >> px >> py >> pz;
	newEvent->addParticle(e*1.e-3,px*1.e-3,py*1.e-3,pz*1.e-3);
      }
      if (!currentStream.fail()) {
	evtList.add(newEvent);
	for (parts = 0; parts < linesToSkip; parts++)
	  currentStream >> e >> px >> py >> pz;
      }
    }
    currentStream.close();
    ++currentFile;
  }
  evtList.rewind();

  return true; // success
}

