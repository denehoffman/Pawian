#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/EtacToapi0Fit/EtacToapi0EventList.hh"
#include "ErrLogger/ErrLogger.hh"


EtacToapi0EventList::EtacToapi0EventList(int kindOfData) :
  _kindOfData(kindOfData),
  _nOfData(20000),
  _nOfMc(50000)
{
  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR");
  _mcPath=theSourcePath+std::string("/Examples/qft++/data/mc100Mgev.dat");

  if (kindOfData==0) _dataPath=theSourcePath+std::string("/Examples/qft++/data/dataEtacToA0Pi_100Mgev.dat");
  else if (kindOfData==2)   _dataPath=theSourcePath+std::string("/Examples/qft++/data/dataSpin2100MgevNew.dat");
  else{
    Alert <<"this kind of data: " << kindOfData << "  is not supported!!!!\n"
	  <<"initialze EtacToapi0EventList either with 0 or with 2 !!!!"
	  << endmsg; 
    exit(1);
  }

  read4Vecs(_dataPath, _nOfData, _data4Vecs);
  read4Vecs(_mcPath,_nOfMc, _mc4Vecs);

}



EtacToapi0EventList::~EtacToapi0EventList()
{
}

void EtacToapi0EventList::read4Vecs(std::string& path,  int nEvts, std::vector<evt4Vec>& the4Vecs)
{
  Info << "calculate 4Vecs and amplitudes for " << path << endmsg;
  std::ifstream inputStream(path.c_str(), std::ios::in);

   if (!inputStream.good()) 
   {
     std::cout << "Input file doesn't exsits !!! ";
     assert(0);
   }

   int counter=0; 
   while (!inputStream.eof() && counter<nEvts )
     {
       if ( counter%1000 == 0 ) Info << "event " << counter << endmsg;
       Vector4<double> pi14V,pi24V,eta4V; // 4-momenta
       get4Vecs(inputStream, pi14V);
       get4Vecs(inputStream, pi24V);
       get4Vecs(inputStream, eta4V);

       Vector4<double>  cm_4V(pi14V+pi24V+eta4V);

       Vector4<double> pi1_cm_4V(pi14V);
       pi1_cm_4V.Boost(cm_4V);
       Vector4<double>  pi2_cm_4V(pi24V);
       pi2_cm_4V.Boost(cm_4V);
       Vector4<double>  eta_cm_4V(eta4V);
       eta_cm_4V.Boost(cm_4V);

       Vector4<double> a2_0_cm_4Vec(pi1_cm_4V+eta_cm_4V);
       Vector4<double> a2_1_cm_4Vec(pi2_cm_4V+eta_cm_4V);
       Vector4<double>  cm_cm_4V(cm_4V);
       cm_cm_4V.Boost(cm_4V);

       evt4Vec theEvt4Vec;
       theEvt4Vec.pi0_0_4Vec=pi1_cm_4V;
       theEvt4Vec.pi0_1_4Vec=pi2_cm_4V;
       theEvt4Vec.eta_4Vec=eta_cm_4V;
       theEvt4Vec.a2_0_4Vec=a2_0_cm_4Vec;
       theEvt4Vec.a2_1_4Vec=a2_1_cm_4Vec;
       theEvt4Vec.cm_4Vec=cm_cm_4V;

       theEvt4Vec.spin2_0_amp = calcSpin2Amp(theEvt4Vec.pi0_0_4Vec, theEvt4Vec.eta_4Vec, theEvt4Vec.pi0_1_4Vec);

       theEvt4Vec.spin2_1_amp = calcSpin2Amp(theEvt4Vec.pi0_1_4Vec, theEvt4Vec.eta_4Vec, theEvt4Vec.pi0_0_4Vec);

       theEvt4Vec.spin1_0_amp = calcSpin1Amp(theEvt4Vec.pi0_0_4Vec, theEvt4Vec.eta_4Vec, theEvt4Vec.pi0_1_4Vec);

       theEvt4Vec.spin1_1_amp = calcSpin1Amp(theEvt4Vec.pi0_1_4Vec, theEvt4Vec.eta_4Vec, theEvt4Vec.pi0_0_4Vec);

       the4Vecs.push_back(theEvt4Vec);
       ++counter;
     } 
}


void EtacToapi0EventList::get4Vecs (std::ifstream& inStream, Vector4<double>& the4Vec)
{
  double tmpPx,tmpPy, tmpPz, tmpE;
  inStream >> tmpPx;
  inStream >> tmpPy;
  inStream >> tmpPz;
  inStream  >> tmpE;
  the4Vec.SetP4(tmpE, tmpPx, tmpPy, tmpPz);  
}


Tensor<complex<double> > EtacToapi0EventList::calcSpin2Amp(Vector4<double>& pi_cm, Vector4<double>& eta_cm, Vector4<double>& pi_recoil_cm)
{

  Tensor<complex<double> > result(0); //rank-0 Tensor for the amplitude

  Vector4<double> inter_cm=eta_cm + pi_cm;
  
  PolVector inter_Pol_cm(2);
  inter_Pol_cm.SetP4(inter_cm,inter_cm.M());

  OrbitalTensor orb_eta_To_inter_piRecoil(2);
  orb_eta_To_inter_piRecoil.SetP4(inter_cm, pi_recoil_cm);
  
  OrbitalTensor orb_inter_To_pieta(2);
  orb_inter_To_pieta.SetP4(pi_cm, eta_cm);

  result = (orb_inter_To_pieta) | ( inter_Pol_cm.Projector() | orb_eta_To_inter_piRecoil );
 
  return result;
}

Tensor<complex<double> > EtacToapi0EventList::calcSpin1Amp(Vector4<double>& pi_cm, Vector4<double>& eta_cm, Vector4<double>& pi_recoil_cm)
{


  Tensor<complex<double> > result(0); // rank-0 Tensor for the amplitude

  Vector4<double> inter_cm=eta_cm + pi_cm;
  
  PolVector inter_Pol_cm(1);
  inter_Pol_cm.SetP4(inter_cm,inter_cm.M());

  OrbitalTensor orb_eta_To_inter_piRecoil(1);
  orb_eta_To_inter_piRecoil.SetP4(inter_cm, pi_recoil_cm);
  
  OrbitalTensor orb_inter_To_pieta(1);
  orb_inter_To_pieta.SetP4(pi_cm, eta_cm);

  result = (orb_inter_To_pieta) | ( inter_Pol_cm.Projector() | orb_eta_To_inter_piRecoil );
 
  return result;
}

