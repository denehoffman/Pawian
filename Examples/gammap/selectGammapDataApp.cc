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

#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <map>
#include <iterator>
#include <memory>

#include "ErrLogger/ErrLogger.hh"

#include "Utils/PawianCollectionUtils.hh"
#include "Utils/ErrLogUtils.hh"

#include "Event/Event.hh"

int main(int __argc,char *__argv[]){

  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
                    || strcmp( __argv[1], "--help" ) == 0 ) ){
    
    InfoMsg << "switch\n"
         << "-inputPath   path and name of the input file\n"
         << "-outputPath   path and name of the output file"
         << endmsg;
    return 1;
  }

    std::string inputFileName="/data/sleipnir1/bertram/pwaData/gammap/omegap/g12_4v.txt";
    std::string outpuFileName="./out.txt";

    int optind = 1;
  
    while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
      bool found=false;
      std::string sw = __argv[optind];
      if (sw=="-inputPath"){
	optind++;
	inputFileName = __argv[optind];
	found=true;
      }
      if (sw=="-outputPath"){
	optind++;
	outpuFileName = __argv[optind];
	found=true;
      }
      if (!found){
	WarningMsg << "Unknown switch: " 
		<< __argv[optind] << endmsg;
	optind++;
      }
      while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
    }

  //  std::string inputFileName="/data/sleipnir1/bertram/pwaData/gammap/omegap/procsim_4v.txt";
  std::ifstream iStream;
  iStream.open(inputFileName.c_str());
  if (!iStream) {
      Alert << "can not open " << inputFileName ;  // << endmsg;
      exit(1);
    }
  int evtCounter=0;
  std::ofstream oStream;
  oStream.open(outpuFileName.c_str());
  while (!iStream.eof() && evtCounter<50000) {
    double e,px,py,pz;
    double weight;
    iStream >> weight;
    Event* newEvent = new Event();
    for(int parts=0; parts <5; ++parts){
      iStream >> px >> py >> pz >>e;
      newEvent->addParticle(e, px, py, pz);
    }
    Vector4<double> photon4Vec=newEvent->p4(0);
    Vector4<double> proton4Vec=newEvent->p4(1);
    Vector4<double> pip4Vec=newEvent->p4(2);
    Vector4<double> pim4Vec=newEvent->p4(3);
    Vector4<double> pi04Vec=newEvent->p4(4);

    // std::cout << "\npx(gammap): " << photon4Vec.Px() << "\tpx(finalps): " << (proton4Vec+pip4Vec+pim4Vec+pi04Vec).Px() << std::endl;
    // std::cout << "py(gammap): " << photon4Vec.Py() << "\tpy(finalps): " << (proton4Vec+pip4Vec+pim4Vec+pi04Vec).Py() << std::endl;
    // std::cout << "pz(gammap): " << photon4Vec.Pz() << "\tpz(finalps): " << (proton4Vec+pip4Vec+pim4Vec+pi04Vec).Pz() << std::endl;
    // std::cout << "E(gammap): " << photon4Vec.E()+0.93827 << "\tE(finalps): " << (proton4Vec+pip4Vec+pim4Vec+pi04Vec).E() << std::endl;
    // std::cout << "pi0: px=" << pi04Vec.Px() << " py=" << pi04Vec.Py() << " pz=" << pi04Vec.Pz() << " E=" << pi04Vec.E() << " M=" <<  pi04Vec.Mass() << std::endl;

    if( pi04Vec.Mass()>0.100 && pi04Vec.Mass()< 0.170){
      oStream << weight
	      <<"\n" << proton4Vec.Px() << "\t" << proton4Vec.Py() << "\t" << proton4Vec.Pz() << "\t" << proton4Vec.E()
	      <<"\n" << pip4Vec.Px() << "\t" << pip4Vec.Py() << "\t" << pip4Vec.Pz() << "\t" << pip4Vec.E()
	      <<"\n" << pim4Vec.Px() << "\t" << pim4Vec.Py() << "\t" << pim4Vec.Pz() << "\t" << pim4Vec.E()
	      <<"\n" << pi04Vec.Px() << "\t" << pi04Vec.Py() << "\t" << pi04Vec.Pz() << "\t" << pi04Vec.E()
	      << std::endl;

      if(evtCounter%1000==0){
	std::cout << "write event " << evtCounter << " into output file" << std::endl;
      }
      evtCounter++;
    }
    delete newEvent;
  }

  iStream.close();
  oStream.close();  
  return 1;
}


