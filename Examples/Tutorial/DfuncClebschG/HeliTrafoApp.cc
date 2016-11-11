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
  

  cout << "************************************" << endl;
  Print4Vec(ref1Vec, "ref1Vec");
  Print4Vec(mother1Vec, "mother1Vec");
  Print4Vec(recoilMother1Vec, "recoilMother1Vec");
  cout << endl;
  Print4Vec(decA1Vec, "decA1Vec");
  Print4Vec(decA2Vec, "decA2Vec");

  Vector4<double> ref1HeliVec = helicityVec(ref1Vec, mother1Vec, ref1Vec);
  Vector4<double> mother1HeliVec = helicityVec(ref1Vec, mother1Vec, mother1Vec);
  Vector4<double> recoilMother1HeliVec = helicityVec(ref1Vec, mother1Vec, recoilMother1Vec);
  Vector4<double> decA1HeliVec = helicityVec(ref1Vec, mother1Vec, decA1Vec);
  Vector4<double> decA2HeliVec = helicityVec(ref1Vec, mother1Vec, decA2Vec);

  Vector4<double> refAxisPawian(0., 0., 0., 1.);
  Vector4<double> ref1RecoilVec=refAxisPawian-ref1Vec;
  Print4Vec(ref1RecoilVec, "ref1RecoilVec");
  
  Vector4<double> ref1HeliVecNew = KinUtils::heliVec(refAxisPawian, ref1Vec, mother1Vec, ref1Vec);
  Vector4<double> mother1HeliVecNew = KinUtils::heliVec(refAxisPawian, ref1Vec, mother1Vec, mother1Vec);
  Vector4<double> recoilMother1HeliVecNew = KinUtils::heliVec(refAxisPawian, ref1Vec, mother1Vec, recoilMother1Vec);
  Vector4<double> decA1HeliVecNew = KinUtils::heliVec(refAxisPawian, ref1Vec, mother1Vec, decA1Vec);
  Vector4<double> decA2HeliVecNew = KinUtils::heliVec(refAxisPawian, ref1Vec, mother1Vec, decA2Vec);

  Vector4<double> ref1RecoilHeliVecNew = KinUtils::heliVec(refAxisPawian,ref1Vec, mother1Vec, ref1RecoilVec);
  Vector4<double> motherref1HeliVecNew = ref1HeliVecNew+ref1RecoilHeliVecNew;

  Vector4<double> ref1HeliVecEvtGen = KinUtils::heliEvtGenVec(refAxisPawian, ref1Vec, mother1Vec, ref1Vec);
  Vector4<double> mother1HeliVecEvtGen = KinUtils::heliEvtGenVec(refAxisPawian, ref1Vec, mother1Vec, mother1Vec);
  Vector4<double> recoilMother1HeliVecEvtGen = KinUtils::heliEvtGenVec(refAxisPawian, ref1Vec, mother1Vec, recoilMother1Vec);
  Vector4<double> decA1HeliVecEvtGen = KinUtils::heliEvtGenVec(refAxisPawian, ref1Vec, mother1Vec, decA1Vec);
  Vector4<double> decA2HeliVecEvtGen = KinUtils::heliEvtGenVec(refAxisPawian, ref1Vec, mother1Vec, decA2Vec);

  Vector4<double> ref1RecoilHeliVecEvtGen = KinUtils::heliEvtGenVec(refAxisPawian,ref1Vec, mother1Vec, ref1RecoilVec);
  Vector4<double> motherref1HeliVecEvtGen = ref1HeliVecEvtGen+ref1RecoilHeliVecEvtGen;
  
  cout << endl;
  Print4Vec(ref1HeliVec, "ref1HeliVec");
  Print4Vec(ref1HeliVecEvtGen, "ref1HeliVecEvtGen");
  Print4Vec(ref1HeliVecNew, "ref1HeliVecNew");
  cout << endl;
  Print4Vec(mother1HeliVec, "mother1HeliVec");
  Print4Vec(mother1HeliVecEvtGen, "mother1HeliVecEvtGen");
  Print4Vec(mother1HeliVecNew, "mother1HeliVecNew");
  cout << endl;
  Print4Vec(recoilMother1HeliVec, "recoilMother1HeliVec: anti || to z-axis");
  Print4Vec(recoilMother1HeliVecEvtGen, "recoilMother1HeliVecEvtGen: anti || to z-axis");
  Print4Vec(recoilMother1HeliVecNew, "recoilMother1HeliVec: anti || to z-axis");
  cout << endl;
  Print4Vec(decA1HeliVec, "decA1HeliVec");
  Print4Vec(decA1HeliVecEvtGen, "decA1HeliVecEvtGen");
  Print4Vec(decA1HeliVecNew, "decA1HeliVecNew");
  cout << endl;
  Print4Vec(decA2HeliVec, "decA2HeliVec");
  Print4Vec(decA2HeliVecEvtGen, "decA2HeliVecEvtGen");
  Print4Vec(decA2HeliVecNew, "decA2HeliVecNew");

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

  Vector4<double> crossMotherRefMotherHeliNew=KinUtils::perpTo(motherref1HeliVecNew, refAxisPawian); 
  Vector4<double> crossNormMotherRefMother_MotherHeliNew=KinUtils::perpTo(crossMotherRefMotherHeliNew, refAxisPawian);
  cout << "\ncrossMotherRefMotherHeliNew (y-direction): " << crossMotherRefMotherHeliNew << endl;
  cout << "crossNormMotherRefMother_MotherHeliNew (x-direction): " << crossNormMotherRefMother_MotherHeliNew << "\n" << endl; 
    

  cout << "\n************************************" << endl;
  cout << "rotate reference vector!!!!" << endl; 

  Vector4<double> refRotRefVec=ref1Vec;
  KinUtils::SetTheta(refRotRefVec, 25./PawianConstants::radToDeg);
  KinUtils::SetPhi(refRotRefVec, -23./PawianConstants::radToDeg);

  Vector4<double> motherRotRefVec(2., 0.3, 0.4, 0.5);
  Vector4<double> recoilMotherRotRefVec = refRotRefVec-motherRotRefVec;

  Vector4<double> decA1RotRefVec=decA1Vec;
  Vector4<double> decA2RotRefVec=motherRotRefVec-decA1RotRefVec;
  cout << "************************************" << endl;
  Print4Vec(refRotRefVec, "refRotRefVec");
  Print4Vec(motherRotRefVec, "motherRotRefVec");
  Print4Vec(recoilMotherRotRefVec, "recoilMotherRotRefVec");
  cout << "************************************" << endl;
  cout << endl;
  Print4Vec(decA1RotRefVec, "decA1RotRefVec");
  Print4Vec(decA2RotRefVec, "decA2RotRefVec");


  Vector4<double> refRotRefHeliVec = helicityVec(refRotRefVec, motherRotRefVec, refRotRefVec);
  Vector4<double> motherRotRefHeliVec = helicityVec(refRotRefVec, motherRotRefVec, motherRotRefVec);
  Vector4<double> recoilMotherRotRefHeliVec = helicityVec(refRotRefVec, motherRotRefVec, recoilMotherRotRefVec);
  Vector4<double> decA1RotRefHeliVec = helicityVec(refRotRefVec, motherRotRefVec, decA1RotRefVec);
  Vector4<double> decA2RotRefHeliVec = helicityVec(refRotRefVec, motherRotRefVec, decA2RotRefVec);

  Vector4<double> refRotRefRecoilVec=refAxisPawian-refRotRefVec;
  Print4Vec(refRotRefRecoilVec, "refRotRefRecoilVec");
  
  Vector4<double> refRotRefHeliVecNew = KinUtils::heliVec(refAxisPawian, refRotRefVec, motherRotRefVec, refRotRefVec);
  Vector4<double> motherRotRefHeliVecNew = KinUtils::heliVec(refAxisPawian, refRotRefVec, motherRotRefVec, motherRotRefVec);
  Vector4<double> recoilMotherRotRefHeliVecNew = KinUtils::heliVec(refAxisPawian, refRotRefVec, motherRotRefVec, recoilMotherRotRefVec);
  Vector4<double> decA1RotRefHeliVecNew = KinUtils::heliVec(refAxisPawian, refRotRefVec, motherRotRefVec, decA1RotRefVec);
  Vector4<double> decA2RotRefHeliVecNew = KinUtils::heliVec(refAxisPawian, refRotRefVec, motherRotRefVec, decA2RotRefVec);

  Vector4<double> refRotRefRecoilHeliVecNew = KinUtils::heliVec(refAxisPawian,refRotRefVec, motherRotRefVec, refRotRefRecoilVec);
  Vector4<double> motherrefRotRefHeliVecNew = refRotRefHeliVecNew+refRotRefRecoilHeliVecNew;

  Vector4<double> refRotRefHeliVecEvtGen = KinUtils::heliEvtGenVec(refAxisPawian, refRotRefVec, motherRotRefVec, refRotRefVec);
  Vector4<double> motherRotRefHeliVecEvtGen = KinUtils::heliEvtGenVec(refAxisPawian, refRotRefVec, motherRotRefVec, motherRotRefVec);
  Vector4<double> recoilMotherRotRefHeliVecEvtGen = KinUtils::heliEvtGenVec(refAxisPawian, refRotRefVec, motherRotRefVec, recoilMotherRotRefVec);
  Vector4<double> decA1RotRefHeliVecEvtGen = KinUtils::heliEvtGenVec(refAxisPawian, refRotRefVec, motherRotRefVec, decA1RotRefVec);
  Vector4<double> decA2RotRefHeliVecEvtGen = KinUtils::heliEvtGenVec(refAxisPawian, refRotRefVec, motherRotRefVec, decA2RotRefVec);

  Vector4<double> refRotRefRecoilHeliVecEvtGen = KinUtils::heliEvtGenVec(refAxisPawian,refRotRefVec, motherRotRefVec, refRotRefRecoilVec);
  Vector4<double> motherrefRotRefHeliVecEvtGen = refRotRefHeliVecEvtGen+refRotRefRecoilHeliVecEvtGen;
  
  cout << endl;
  Print4Vec(refRotRefHeliVec, "refRotRefHeliVec");
  Print4Vec(refRotRefHeliVecEvtGen, "refRotRefHeliVecEvtGen");
  Print4Vec(refRotRefHeliVecNew, "refRotRefHeliVecNew");
  cout << endl;
  Print4Vec(motherRotRefHeliVec, "motherRotRefHeliVec");
  Print4Vec(motherRotRefHeliVecEvtGen, "motherRotRefHeliVecEvtGen");
  Print4Vec(motherRotRefHeliVecNew, "motherRotRefHeliVecNew");
  cout << endl;
  Print4Vec(recoilMotherRotRefHeliVec, "recoilMotherRotRefHeliVec: anti || to z-axis");
  Print4Vec(recoilMotherRotRefHeliVecEvtGen, "recoilMotherRotRefHeliVecEvtGen: anti || to z-axis");
  Print4Vec(recoilMotherRotRefHeliVecNew, "recoilMotherRotRefHeliVec: anti || to z-axis");
  cout << endl;
  Print4Vec(decA1RotRefHeliVec, "decA1RotRefHeliVec");
  Print4Vec(decA1RotRefHeliVecEvtGen, "decA1RotRefHeliVecEvtGen");
  Print4Vec(decA1RotRefHeliVecNew, "decA1RotRefHeliVecNew");
  cout << endl;
  Print4Vec(decA2RotRefHeliVec, "decA2RotRefHeliVec");
  Print4Vec(decA2RotRefHeliVecEvtGen, "decA2RotRefHeliVecEvtGen");
  Print4Vec(decA2RotRefHeliVecNew, "decA2RotRefHeliVecNew");

  Vector4<double> crossMotherRefMotherRotRefHeliNew=KinUtils::perpTo(motherrefRotRefHeliVecNew, refAxisPawian); 
  Vector4<double> crossNormMotherRefMother_MotherRotRefHeliNew=KinUtils::perpTo(crossMotherRefMotherRotRefHeliNew, refAxisPawian);

  cout << "\ncrossMotherRefMotherRotRefHeliNew (y-direction): " << crossMotherRefMotherRotRefHeliNew << endl;
  cout << "crossNormMotherRefMother_MotherRotRefHeliNew (x-direction): " << crossNormMotherRefMother_MotherRotRefHeliNew << "\n" << endl;


   cout << "\n************************************" << endl;
  cout << "rotate mother reference vector!!!!" << endl;
  Vector4<double> mother_MotherRefRot_Axis=ref1Vec+Vector4<double>(1., 0.2, 0.2, 0.);

  cout << "************************************" << endl;
  Print4Vec(ref1Vec, "ref1Vec");
  Print4Vec(mother1Vec, "mother1Vec");
  Print4Vec(recoilMother1Vec, "recoilMother1Vec");
  cout << endl;
  Print4Vec(decA1Vec, "decA1Vec");
  Print4Vec(decA2Vec, "decA2Vec");

  Vector4<double> ref1Recoil_MotherRefRot_Vec=mother_MotherRefRot_Axis-ref1Vec;
  Print4Vec(mother_MotherRefRot_Axis, "mother_MotherRefRot_Axis");
  Print4Vec(ref1Recoil_MotherRefRot_Vec, "ref1RecoilRotVec");

  Vector4<double> ref1_MotherRefRot_HeliVecNew = KinUtils::heliVec(mother_MotherRefRot_Axis, ref1Vec, mother1Vec, ref1Vec);
  Vector4<double> mother1_MotherRefRot_HeliVecNew = KinUtils::heliVec(mother_MotherRefRot_Axis, ref1Vec, mother1Vec, mother1Vec);
  Vector4<double> recoilMother1_MotherRefRot_HeliVecNew = KinUtils::heliVec(mother_MotherRefRot_Axis, ref1Vec, mother1Vec, recoilMother1Vec);
  Vector4<double> decA1_MotherRefRot_HeliVecNew = KinUtils::heliVec(mother_MotherRefRot_Axis, ref1Vec, mother1Vec, decA1Vec);
  Vector4<double> decA2_MotherRefRot_HeliVecNew = KinUtils::heliVec(mother_MotherRefRot_Axis, ref1Vec, mother1Vec, decA2Vec);

  Vector4<double> ref1Recoil_MotherRefRot_HeliVecNew = KinUtils::heliVec(mother_MotherRefRot_Axis,ref1Vec, mother1Vec, ref1Recoil_MotherRefRot_Vec);
  Vector4<double> motherref1_MotherRefRot_HeliVecNew = ref1_MotherRefRot_HeliVecNew+ref1Recoil_MotherRefRot_HeliVecNew;

  Vector4<double> ref1Heli_MotherRefRot_VecEvtGen = KinUtils::heliEvtGenVec(mother_MotherRefRot_Axis, ref1Vec, mother1Vec, ref1Vec);
  Vector4<double> mother1Heli_MotherRefRot_VecEvtGen = KinUtils::heliEvtGenVec(mother_MotherRefRot_Axis, ref1Vec, mother1Vec, mother1Vec);
  Vector4<double> recoilMother1_MotherRefRot_HeliVecEvtGen = KinUtils::heliEvtGenVec(mother_MotherRefRot_Axis, ref1Vec, mother1Vec, recoilMother1Vec);
  Vector4<double> decA1Heli_MotherRefRot_VecEvtGen = KinUtils::heliEvtGenVec(mother_MotherRefRot_Axis, ref1Vec, mother1Vec, decA1Vec);
  Vector4<double> decA2Heli_MotherRefRot_VecEvtGen = KinUtils::heliEvtGenVec(mother_MotherRefRot_Axis, ref1Vec, mother1Vec, decA2Vec);

  Vector4<double> ref1Recoil_MotherRefRot_HeliVecEvtGen = KinUtils::heliEvtGenVec(mother_MotherRefRot_Axis,ref1Vec, mother1Vec, ref1RecoilVec);
  Vector4<double> motherref1_MotherRefRot_HeliVecEvtGen = ref1HeliVecEvtGen+ref1RecoilHeliVecEvtGen;

  cout << endl;
  Print4Vec(ref1HeliVec, "ref1HeliVec");
  Print4Vec(ref1Heli_MotherRefRot_VecEvtGen, "ref1Heli_MotherRefRot_VecEvtGen");
  Print4Vec(ref1_MotherRefRot_HeliVecNew, "ref1_MotherRefRot_HeliVecNew");
  cout << endl;
  Print4Vec(mother1HeliVec, "mother1HeliVec");
  Print4Vec(mother1Heli_MotherRefRot_VecEvtGen, "mother1Heli_MotherRefRot_VecEvtGen");
  Print4Vec(mother1_MotherRefRot_HeliVecNew, "mother1_MotherRefRot_HeliVecNew");
  cout << endl;
  Print4Vec(recoilMother1HeliVec, "recoilMother1HeliVec: anti || to z-axis");
  Print4Vec(recoilMother1_MotherRefRot_HeliVecEvtGen, "recoilMother1_MotherRefRot_HeliVecEvtGen: anti || to z-axis");
  Print4Vec(recoilMother1_MotherRefRot_HeliVecNew, "recoilMother1_MotherRefRot_HeliVecNew: anti || to z-axis");
  cout << endl;
  Print4Vec(decA1HeliVec, "decA1HeliVec");
  Print4Vec(decA1Heli_MotherRefRot_VecEvtGen, "decA1Heli_MotherRefRot_VecEvtGen");
  Print4Vec(decA1_MotherRefRot_HeliVecNew, "decA1_MotherRefRot_HeliVecNew");
  cout << endl;
  Print4Vec(decA2HeliVec, "decA2HeliVec");
  Print4Vec(decA2Heli_MotherRefRot_VecEvtGen, "decA2Heli_MotherRefRot_VecEvtGen");
  Print4Vec(decA2_MotherRefRot_HeliVecNew, "decA2_MotherRefRot_HeliVecNew");

  Vector4<double> crossMotherRefMother_MotherRefRot_HeliNew=KinUtils::perpTo(ref1Recoil_MotherRefRot_HeliVecNew, refAxisPawian); //refAxisPawian = z-axis
  Vector4<double> crossNormMotherRefMother_Mother_MotherRefRot_HeliNew=KinUtils::perpTo(crossMotherRefMother_MotherRefRot_HeliNew, refAxisPawian);
  // Vector4<double> crossMotherRefMother_MotherRefRot_HeliNew=KinUtils::perpTo(recoilMother1_MotherRefRot_HeliVecNew, motherref1_MotherRefRot_HeliVecNew); 
  // Vector4<double> crossNormMotherRefMother_Mother_MotherRefRot_HeliNew=KinUtils::perpTo(crossMotherRefMother_MotherRefRot_HeliNew, motherref1_MotherRefRot_HeliVecNew);
  cout << "\ncrossMotherRefMother_MotherRefRot_HeliNew (y-direction): " << crossMotherRefMother_MotherRefRot_HeliNew << endl;
  cout << "crossNormMotherRefMother_Mother_MotherRefRot_HeliNew (x-direction): " << crossNormMotherRefMother_Mother_MotherRefRot_HeliNew << "\n" << endl;  
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
