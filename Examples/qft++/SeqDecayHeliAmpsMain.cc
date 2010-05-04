#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include "Examples/qft++/SeqDecHeliAmps.hh"

int findParityAndRemove(std::string& theString){
  int result=0;
  std::string plusStr("+");
  std::string minusStr("-");
  size_t strfound=0;
  if ( (strfound=theString.find(plusStr)) != std::string::npos){
   theString.erase (strfound, strfound+1);
   result++;
  }
  if ( (strfound=theString.find(minusStr)) != std::string::npos){
   theString.erase (strfound, strfound+1);
   result--;
  }
  return result;
}

int main(int __argc,char *__argv[]){

  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
		    || strcmp( __argv[1], "--help" ) == 0 ) ){

    std::cout << "This application calculates all helicity amplitudes for the decay chain:" << std::endl;
    std::cout << " mother -> intermediate + recoil  ; intermediate -> daughter1 + daughter2" << std::endl << std::endl;
    std::cout << "The final state particles are: " << "\n" 
              << "recoil= pi0 --> JP=0-" << "\n" 
              << "daughter1 = eta --> JP=0-" << "\n"
	      << "daughter2 = pi0 --> JP=0-" << "\n"
	      << "the mass of the intermediate resonance will be taken into acccount with a simple Breit-Wigner function (m=1.320 GeV, sigma(m)=0.100 GeV)" << std::endl << std::endl;  
    std::cout << "The input parameters have to be passed with the following arguments: " << "\n"
	      << "-m \"JP\" : spin-parity of the mother particle" << "\n"
	      << "-i \"JP\" : spin-parity of the intermediate resonance" << "\n"
	      << "-motherJzs \"M1 M2 M3 ...\": contributed Jz's of the mother particle" << "\n"
	      << "to start the application one has to type e.g.:" << "\n"
	      << "./SeqDecayHeliAmpsApp -m \"1-\" -i \"2+\" -motherJzs \"-1 1\" "  << "\n"  << "\n"
	      << "with \"-v\" one can switch into the verbose mode" << std::endl;

    return 0;
  }

  bool verboseMode=false;
  std::string  jpMother="";
  std::string  jpIntermediate="";
  std::string  motherJzs="";

  int optind=1;
  // decode arguments
  while ((optind < __argc) && (__argv[optind][0]=='-')) {
    
    std::string sw = __argv[optind];
    if (sw=="-m") {
      optind++;
      jpMother = __argv[optind];
    }
    else if (sw=="-i") {
      optind++;
      jpIntermediate = __argv[optind];
    }
    else if (sw=="-motherJzs") {
      optind++;
      motherJzs = __argv[optind];
    }
    else if (sw=="-v") {
      optind++;
      verboseMode=true;
    }
    else
      cout << "Unknown switch: " 
	   << __argv[optind] << endl;
    optind++;
  }

  std::stringstream motherJzsStream(motherJzs);
 
  std::vector<double> motherJzVec;

  if ( motherJzsStream.good() ){
    while ( ! motherJzsStream.eof() ){
      double tmpjz=-100.;
      motherJzsStream >> tmpjz;
      motherJzVec.push_back(tmpjz);    
    }
  }
  

  int pMother=findParityAndRemove(jpMother);
  if ( fabs(pMother)!=1){
    std::cout << "unvalid parity for mother particle!!!";
    return 0;    
  }

  std::stringstream jpMotherStream(jpMother);
  int jMother=-1;
  jpMotherStream >> jMother;  

  if ( jMother<0 || jMother>9){
    std::cout << "Spin of the mother particle negative or >9 ";
    return 0;    
  }

  std::cout << "Mother Resonance: J=" << jMother  << " P=" << pMother << std::endl;  


  int pInter=findParityAndRemove(jpIntermediate);
  if ( fabs(pInter)!=1){
    std::cout << "unvalid parity for the intermediate resonance!!!";
    return 0;    
  }

  std::stringstream jpInterStream(jpIntermediate);
  int jInter=-1;
  jpInterStream >> jInter;  

  if ( jInter<0 || jInter>9){
    std::cout << "Spin of the intermediate resonance negative or >9 ";
    return 0;    
  }

  std::cout << "Intermediate Resonance: J=" << jInter  << " P=" << pInter << std::endl;


  jpChain motherChain;
  motherChain.mother.J=jMother;
  motherChain.mother.P=pMother;
  motherChain.daughter1.J=jInter;
  motherChain.daughter1.P=pInter;
  motherChain.daughter2.J=0;
  motherChain.daughter2.P=-1;

  jpChain daughterChain;
  daughterChain.mother.J=motherChain.daughter1.J;
  daughterChain.mother.P=motherChain.daughter1.P;
  daughterChain.daughter1.J=0.;
  daughterChain.daughter1.P=-1;
  daughterChain.daughter2.J=0.;
  daughterChain.daughter2.P=-1;  

  SeqDecHeliAmps(motherChain, daughterChain, motherJzVec, verboseMode);
  return 0;
}

