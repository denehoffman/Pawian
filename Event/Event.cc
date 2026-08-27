//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

#include "Event/Event.hh"
#include <vector>
#include <fstream>

Event::Event() : weight(1.),
  polMag(0.),
  polAngle(0.),
  evtNo(0)
{}
Event::Event(int evtNumber) :
  weight(1.),
  polMag(0.),
  polAngle(0.),
  dataPoint(0.),
  dataPointErr(0.),
  // phase(0.),
  // phaseErr(0.),
  // eta(0.),
  // etaErr(0.),
  evtNo(evtNumber)
{}
Event::~Event()
{
  std::vector<EvtPartData*>::iterator iter;
  for (iter = particles.begin(); iter != particles.end(); ++iter)
    delete *iter;
}

Vector4<float>* Event::p4(unsigned int i)
{
  if (particles.size() > i)
    return particles[i]->vector4;
  else {
    Alert << "accessing " << i << "th of "
	  << particles.size() << " particles" << endmsg;
    exit(1);
  }
  return 0;
}


float* Event::pid(unsigned int i)
{
  if (particles.size() > i)
    return particles[i]->pidVector;
  else {
    Alert << "accessing pid of " << i << "th of "
	  << particles.size() << " particles" << endmsg;
    exit(1);
  }
}


void Event::addParticle(double e, double px, double py, double pz)
{
  EvtPartData* evtData = new EvtPartData(e,px,py,pz);
  particles.push_back(evtData);
  return;
}

// void Event::addScatterInfo(double thePhase, double thePhaseErr, double theEta, double theEtaErr){
//   phase=thePhase;
//   phaseErr=thePhaseErr;
//   eta=theEta;
//   etaErr=theEtaErr;
// }

void Event::addPolarizationInfo(double thePolMag, double thePolAngle) {
  polMag = thePolMag;
  polAngle = thePolAngle;
}

void Event::addScatterInfo(double theDataPoint, double theDataPointErr){
  dataPoint=theDataPoint;
  dataPointErr=theDataPointErr;
}

void Event::addWeight(double theWeight)
{
  weight = theWeight;
}

int Event::size()
{
  return particles.size();
}

void Event::dumpEvt(std::ostream& os, bool useMeV, std::string orderInFile){
  os << weight << std::endl;
  os << "Polarization: "<< polMag * 100. << "%" << "\t" << polAngle << " rad" << std::endl;
  for(unsigned int i=0; i<particles.size(); ++i){
    Vector4<float>* p4Vec=p4(i);
    if(useMeV){
      if (orderInFile=="Px Py Pz E") os << p4Vec->Px()*1000. << "\t" <<  p4Vec->Py()*1000. << "\t" <<  p4Vec->Pz()*1000. << "\t" << p4Vec->E()*1000.  << std::endl;
      else os << p4Vec->E()*1000. << "\t" << p4Vec->Px()*1000. << "\t" <<  p4Vec->Py()*1000. << "\t" <<  p4Vec->Pz()*1000. << std::endl;
    }
    else{
    if (orderInFile=="Px Py Pz E") os << p4Vec->Px() << "\t" <<  p4Vec->Py() << "\t" <<\
  p4Vec->Pz() << "\t" << p4Vec->E()  << std::endl;
    else os << p4Vec->E() << "\t" << p4Vec->Px() << "\t" <<	p4Vec->Py() << "\t" <<	p4Vec->Pz() << std::endl;
    }
  }

}

bool Event::operator<(const Event& compare) const{
  bool result=false;
  if ( evtNo < compare.eventNo()) result=true;
  return result;
}
