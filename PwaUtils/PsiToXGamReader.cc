#include "PwaUtils/PsiToXGamReader.hh"

#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"

PsiToXGamReader::PsiToXGamReader()
{}

PsiToXGamReader::PsiToXGamReader(const std::vector<std::string>& files, int particles, int skip, bool useWeight):
  EventReader(useWeight),
  numParticles(particles),
  linesToSkip(skip),
  _useMassRange(false),
  _massRange(std::make_pair(0.,10.) )
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

PsiToXGamReader::~PsiToXGamReader()
{}

bool PsiToXGamReader::fillAll(EventList& evtList)
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
	currentStream >> px >> py >> pz >> e;
        newEvent->addParticle(e,px,py,pz);
	Vector4<double> tmp = newEvent->p4(parts);
	if(parts>0) fvX= fvX+tmp;
      }
      
      if(_useMassRange){
	if(fvX.Mass()<_massRange.first || fvX.Mass()>_massRange.second  ) continue;
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

