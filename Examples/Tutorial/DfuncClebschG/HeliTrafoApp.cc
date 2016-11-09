//************************************************************************//
//                                                                        //
//  Copyright 2016 Bertram Kopf (bertram@ep1.rub.de)                      //
//                 - Ruhr-Universit??t Bochum                              //
//                                                                        //
//  This file is part of Pawian.                                          //
//                                                                        //
//  Pawian is free software: you can redistribute it and/or modify        //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or     //
//  (at your option) any later version.                                   //
//                                                                        //
//  Pawian is distributed in the hope that it will be useful,             //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of        //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
//  GNU General Public License for more details.                          //
//                                                                        //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.       //
//                                                                        //
//************************************************************************//

#include <getopt.h>
#include <math.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/KinUtils.hh"
#include "Utils/PawianConstants.hh"
#include "qft++/topincludes/c++-template-utils.hh"
#include "PspGen/EvtKine.hh"
#include "PspGen/EvtVector4R.hh"

using namespace std;
//_____________________________________________________________________________

/// Prints program usage to the screen
void PrintUsage();

void Print4Vec(Vector4<double>& the4Vec, std::string name="current4Vec");

/// Prints a line accross the screen
void PrintLine(char __c){ 
  for(int i = 0; i < 80; i++) cout << __c; 
  cout << endl;
}
//_____________________________________________________________________________

int main(int __argc,char *__argv[]){

  /*__________________________Parse the Command Line_________________________*/
  int c;
  //extern char* optarg;
  //extern int optind;
  
  while((c = getopt(__argc,__argv,"h")) != -1){
    switch(c){
    case 'h': // help option
      PrintUsage();
      return EXIT_SUCCESS;
      break;
    default:
      break;
    }
  }  

  //creating 4 vectors
  Vector4<double> ref1Vec(3., 0., 0., 1.3);
  Vector4<double> mother1Vec(2., 0.3, 0.4, 0.5);
  Vector4<double> recoilMother1Vec = ref1Vec-mother1Vec;

  double mDecA1=1.;
  double pDecA1=0.3;
  Vector4<double> decA1Vec(sqrt(mDecA1*mDecA1+pDecA1*pDecA1), 0., 0., pDecA1);
  KinUtils::SetTheta(decA1Vec, 40./PawianConstants::radToDeg);
  KinUtils::SetPhi(decA1Vec, -115./PawianConstants::radToDeg);
  Vector4<double> decA2Vec=mother1Vec-decA1Vec;
  

  double mDecB1=1.;
  double pDecB1=0.3;
  Vector4<double> decB1Vec(sqrt(mDecB1*mDecB1+pDecB1*pDecB1), 0., 0., pDecB1);
  Vector4<double> decB2Vec = mother1Vec-decB1Vec;

  double eDecC1=sqrt(1.+0.5*0.5);
  double pDecC1=0.5;

  
  Vector4<double> decC1Vec(eDecC1, 0., 0., pDecC1); //parallel to mother1Vec 
  KinUtils::SetTheta(decC1Vec, mother1Vec.Theta());
  KinUtils::SetPhi(decC1Vec, mother1Vec.Phi());
  Vector4<double> decC2Vec = mother1Vec-decC1Vec;

  cout << "************************************" << endl;
  Print4Vec(ref1Vec, "ref1Vec");
  Print4Vec(mother1Vec, "mother1Vec");
  Print4Vec(recoilMother1Vec, "recoilMother1Vec");
  cout << endl;
  Print4Vec(decA1Vec, "decA1Vec");
  Print4Vec(decA2Vec, "decA2Vec");
  cout << endl;
  Print4Vec(decB1Vec, "decB1Vec || to z ref axis");
  Print4Vec(decB2Vec, "decB2Vec");
  cout << endl;
  Print4Vec(decC1Vec, "decC1Vec || to mother1");
  Print4Vec(decC2Vec, "decC2Vec");

  Vector4<double> ref1HeliVec = helicityVec(ref1Vec, mother1Vec, ref1Vec);
  Vector4<double> mother1HeliVec = helicityVec(ref1Vec, mother1Vec, mother1Vec);
  Vector4<double> recoilMother1HeliVec = helicityVec(ref1Vec, mother1Vec, recoilMother1Vec);
  Vector4<double> decA1HeliVec = helicityVec(ref1Vec, mother1Vec, decA1Vec);
  Vector4<double> decA2HeliVec = helicityVec(ref1Vec, mother1Vec, decA2Vec);
  Vector4<double> decB1HeliVec = helicityVec(ref1Vec, mother1Vec, decB1Vec);
  Vector4<double> decB2HeliVec = helicityVec(ref1Vec, mother1Vec, decB2Vec);
  Vector4<double> decC1HeliVec = helicityVec(ref1Vec, mother1Vec, decC1Vec);
  Vector4<double> decC2HeliVec = helicityVec(ref1Vec, mother1Vec, decC2Vec);

  Print4Vec(ref1HeliVec, "ref1HeliVec");
  Print4Vec(mother1HeliVec, "mother1HeliVec");
  Print4Vec(recoilMother1HeliVec, "recoilMother1HeliVec: anti || to z-axis");
  cout << endl;
  Print4Vec(decA1HeliVec, "decA1HeliVec");
  Print4Vec(decA2HeliVec, "decA2HeliVec");

  cout << "\ncheck theta and Phi with EvtGen:" << endl;
  EvtVector4R ref1EvtGen(ref1Vec.E(), ref1Vec.Px(), ref1Vec.Py(), ref1Vec.Pz());
  EvtVector4R mother1VecEvtGen(mother1Vec.E(), mother1Vec.Px(), mother1Vec.Py(), mother1Vec.Pz());
  EvtVector4R decA1VecEvtGen(decA1Vec.E(), decA1Vec.Px(), decA1Vec.Py(), decA1Vec.Pz());

  double cosThetadecA1Heli = costDecHeli(ref1Vec, mother1Vec, decA1Vec);
  cout << "thetaA1Heli lorinv: " << acos(cosThetadecA1Heli)*PawianConstants::radToDeg << endl;

  double costhetaA1Heli =  EvtDecayAngle(ref1EvtGen, mother1VecEvtGen, decA1VecEvtGen);
  cout << "thetaA1Heli EvtGen: " << acos(costhetaA1Heli)*PawianConstants::radToDeg << endl;

  EvtVector4R refAxis(0., 0., 0., 1.);
  double phiA1Heli =  EvtDecayAnglePhi(refAxis, ref1EvtGen, mother1VecEvtGen, decA1VecEvtGen);
  cout << "phiA1Heli EvtGen: " << phiA1Heli*PawianConstants::radToDeg << endl;  

  EvtVector4R ref1aEvtGen(ref1EvtGen.get(0)/2., ref1EvtGen.get(1)/2., ref1EvtGen.get(2)/2, ref1EvtGen.get(3)/2.+0.5);
  double phiA1Heli1a =  EvtDecayAnglePhi(ref1aEvtGen, ref1EvtGen, mother1VecEvtGen, decA1VecEvtGen);
  cout << "phiA1Heli1a EvtGen: " << phiA1Heli1a*PawianConstants::radToDeg << endl;  

  Vector4<double> refAxisPawian(0., 0., 0., 1.);
  Vector4<double> ref1RecoilVec=refAxisPawian-ref1Vec;
  Print4Vec(ref1RecoilVec, "ref1RecoilVec");
  Vector4<double> decA1HeliVecNew = KinUtils::heliVec(refAxisPawian, ref1Vec, mother1Vec, decA1Vec);
  Vector4<double> recoilMother1HeliVecNew = KinUtils::heliVec(refAxisPawian, ref1Vec, mother1Vec, recoilMother1Vec);
  Vector4<double> ref1HeliVecNew = KinUtils::heliVec(refAxisPawian,ref1Vec, mother1Vec, ref1Vec);
  Vector4<double> ref1RecoilHeliVecNew = KinUtils::heliVec(refAxisPawian,ref1Vec, mother1Vec, ref1RecoilVec);
  Vector4<double> motherref1HeliVecNew = ref1HeliVecNew+ref1RecoilHeliVecNew;

  Print4Vec(decA1HeliVecNew, "decA1HeliVecNew");
  Print4Vec(recoilMother1HeliVecNew, "recoilMother1HeliVecNew");
  Print4Vec(ref1HeliVecNew, "ref1HeliVecNew");
  Print4Vec(ref1RecoilHeliVecNew, "ref1RecoilHeliVecNew");

  Vector4<double> crossMotherRefMotherHeli=KinUtils::perpTo(motherref1HeliVecNew, refAxisPawian); 
  Vector4<double> crossNormMotherRefMother_MotherHeli=KinUtils::perpTo(crossMotherRefMotherHeli, refAxisPawian);
  cout << "crossMotherRefMotherHeli (y-direction): " << crossMotherRefMotherHeli << "\n" << endl;
  cout << "crossNormMotherRefMother_MotherHeli (x-direction): " << crossNormMotherRefMother_MotherHeli << "\n" << endl; 
    
  Vector4<double> decA1HeliVecEvtGenNew = KinUtils::heliEvtGenVec(refAxisPawian, ref1Vec, mother1Vec, decA1Vec);
  Vector4<double> recoilMother1HeliVecEvtGenNew = KinUtils::heliEvtGenVec(refAxisPawian, ref1Vec, mother1Vec, recoilMother1Vec);
  Vector4<double> ref1HeliVecEvtGenNew = KinUtils::heliEvtGenVec(refAxisPawian,ref1Vec, mother1Vec, ref1Vec);
  Vector4<double> ref1RecoilHeliVecEvtGenNew = KinUtils::heliEvtGenVec(refAxisPawian,ref1Vec, mother1Vec, ref1RecoilVec);
  Vector4<double> motherref1HeliVecEvtGenNew = ref1HeliVecEvtGenNew+ref1RecoilHeliVecEvtGenNew;
  Print4Vec(decA1HeliVecEvtGenNew, "decA1HeliVecEvtGenNew");
  Print4Vec(recoilMother1HeliVecEvtGenNew, "recoilMother1HeliVecEvtGenNew");
  Print4Vec(ref1HeliVecEvtGenNew, "ref1HeliVecEvtGenNew");
  Print4Vec(ref1RecoilHeliVecEvtGenNew, "ref1RecoilHeliVecEvtGenNew");
  Print4Vec(motherref1HeliVecEvtGenNew, "motherref1HeliVecEvtGenNew");

  Vector4<double> crossMotherRefMotherHeli1=KinUtils::perpTo(motherref1HeliVecEvtGenNew, refAxisPawian); 
  Vector4<double> cross_crossMotherRefMotherHeli_MotherHeli1=KinUtils::perpTo(crossMotherRefMotherHeli1, refAxisPawian);
  cout << "crossMotherRefMotherHeli1 (y-direction): " << crossMotherRefMotherHeli1 << "\n" << endl;
  cout << "cross_crossMotherRefMotherHeli_MotherHeli1 (x-direction): " << cross_crossMotherRefMotherHeli_MotherHeli1 << "\n" << endl; 


  cout << "************************************" << endl;
  cout << "************************************" << endl;

  Vector4<double> motherRef2Vec(sqrt(3.1*3.1+1.), 0., 0., 1.);
  KinUtils::SetTheta(motherRef2Vec, 55./PawianConstants::radToDeg);
  KinUtils::SetPhi(motherRef2Vec, -100./PawianConstants::radToDeg); 
  Vector4<double> ref2Vec(2., -0.1, 0.2, .3);
  Vector4<double> ref2VecRecoil=motherRef2Vec-ref2Vec;
  Vector4<double> mother2Vec(1., 0.1, 0.3, -0.2);
  Vector4<double> daughter2Vec(1., 0.1, 0.1, 0.1);


  Vector4<double> ref2VecHeli = KinUtils::heliVec(motherRef2Vec, ref2Vec, mother2Vec, ref2Vec);
  Vector4<double> ref2VecRecoilHeli = KinUtils::heliVec(motherRef2Vec, ref2Vec, mother2Vec, ref2VecRecoil);
  Vector4<double> mother2VecHeli = KinUtils::heliVec(motherRef2Vec, ref2Vec, mother2Vec, mother2Vec);
  Print4Vec(ref2VecHeli, "ref2VecHeli");
  Print4Vec(ref2VecRecoilHeli, "ref2VecRecoilHeli");
  Print4Vec(mother2VecHeli, "mother2VecHeli");


  Vector4<double> ref2VecHeliNew = KinUtils::heliEvtGenVec(motherRef2Vec, ref2Vec, mother2Vec, ref2Vec);
  Vector4<double> ref2VecRecoilHeliNew = KinUtils::heliEvtGenVec(motherRef2Vec, ref2Vec, mother2Vec, ref2VecRecoil);
  Vector4<double> mother2VecHeliNew = KinUtils::heliEvtGenVec(motherRef2Vec, ref2Vec, mother2Vec, mother2Vec);
  Print4Vec(ref2VecHeliNew, "ref2VecHeliNew");
  Print4Vec(ref2VecRecoilHeliNew, "ref2VecRecoilHeliNew");
  Print4Vec(mother2VecHeliNew, "mother2VecHeliNew");
  return EXIT_SUCCESS;
}

//_____________________________________________________________________________

void PrintUsage(){ 
  cout << "Usage: tensor " << endl;
  cout << "This executable provides a number of example usages of the tensor "
       << "package. Run\nthe executable to see what's being done, then look at"
       << " the source file to see \nhow it's done in the code."
       << endl;
}


void Print4Vec(Vector4<double>& the4Vec, std::string name){
  cout << name << ": " << the4Vec << "\tMass: " << the4Vec.M() 
       << "\t Theta: " << the4Vec.Theta()*PawianConstants::radToDeg << "\t Phi: " << the4Vec.Phi()*PawianConstants::radToDeg << endl;
}

//_____________________________________________________________________________
