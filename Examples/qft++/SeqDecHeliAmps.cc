// 
//
//_____________________________________________________________________________
#include "Examples/qft++/SeqDecHeliAmps.hh"
#include <getopt.h>
#include <cassert>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
//_____________________________________________________________________________
 

SeqDecHeliAmps::SeqDecHeliAmps(jpChain& motherChain, jpChain& daughterChain, std::vector<double>& motherMs, bool verbose):
  _verbose(verbose){

  assert(checkValidity(motherChain, daughterChain, motherMs));

  TFile* theTFile=new TFile("./seqDecHeli.root","recreate");
  TH2F* dalitzAsymWOHist= new TH2F("dalitzAsymWOHist","dalitzplot asym w/o weight",50, 0., 8.,50, 0., 6.8);
  TH2F* dalitzWOHist= new TH2F("dalitzWOHist","dalitzplot w/o weight",50, 0., 8.,50, 0., 8.);  
  TH2F* dalitzHist= new TH2F("dalitzHist","dalitzplot",50, 0., 8.,50, 0., 8.);
  TH1F* invpietaWOHist= new TH1F("invpietaWOHist","M(#pi #eta) w/o weight",120, 0.6, 3.0);
  TH1F* invpipiWOHist= new TH1F("invpipiWOHist","M(#pi #pi) w/o weight",120, 0.2, 2.6);
  TH1F* invpietaHist= new TH1F("invpietaHist","M(#pi #eta)",120, 0.6, 3.0);
  TH1F* invpipiHist= new TH1F("invpipiHist","M(#pi #pi)",120, 0.2, 2.6);
  TH1F* cosThetaMotherWOHist= new TH1F("cosThetaMotherWOHist","cos(#Theta) mother decay w/o weight",200, -1., 1.);
  TH1F* cosThetaMotherHist= new TH1F("cosThetaMotherHist","cos(#Theta) mother decay",200, -1., 1.);
  TH1F* cosThetaInterHist= new TH1F("cosThetaInterHist","cos(#Theta) intermediate resonance decay",200, -1., 1.);
  TH1F* cosThetaInterWOHist= new TH1F("cosThetaInterWOHist","cos(#Theta) intermediate resonance decay w/o weight",200, -1., 1.);

  /*______________________________Set 4-Momenta______________________________*/

  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR");
  std::string inputFile(theSourcePath+"/Examples/qft++/data/mc100Mgev.dat");
  std::ifstream inputStream(inputFile.c_str(), std::ios::in);

   if (!inputStream.good()) 
     {
       std::cout << "Input file: " << inputFile << " doesn't exsits !!! ";
       assert(0);
     }

  std::map<std::string,TH1F*> theMotherHistMap; 
  std::map<std::string,TH1F*> theDaughterHistMap;   
  int counter=0; 
 while (!inputStream.eof() && counter<30000 )
   {
     Vector4<double> pi14V,piRecoil4V,eta4V; // 4-momenta       

     get4Vecs(inputStream, pi14V);     
     get4Vecs(inputStream, piRecoil4V);
     get4Vecs(inputStream, eta4V);
     
     if (_verbose){
       std::cout << "\n"
		 << "*** event no " <<  counter << " ***" << "\n" 
		 << "pi0_1  Mass: " <<  pi14V.Mass() << " px: " << pi14V.Px() << " py: " << pi14V.Px() << " pz: " << pi14V.Pz() << " e: " << pi14V.E() << "\n"
		 << "eta Mass: " <<  eta4V.Mass() << " px: " << eta4V.Px() << " py: " << eta4V.Px() << " pz: " << eta4V.Pz() << " e: " << eta4V.E() <<"\n" 
		 << "piRecoil Mass: " <<  piRecoil4V.Mass() << " px: " << piRecoil4V.Px() << " py: " << piRecoil4V.Px() << " pz: " << piRecoil4V.Pz() << " e: " << piRecoil4V.E() << std::endl;
     }

     Vector4<double>  cm_4V(pi14V+piRecoil4V+eta4V); // 4-momenta
     Vector4<double>  inter_4V(pi14V+eta4V);

     Vector4<double> pi1_cm_4V(pi14V);
     pi1_cm_4V.Boost(cm_4V);
     Vector4<double>  piRecoil_cm_4V(piRecoil4V);
     piRecoil_cm_4V.Boost(cm_4V);
     Vector4<double>  eta_cm_4V(eta4V);
     eta_cm_4V.Boost(cm_4V);
     Vector4<double>  inter_cm_4V(inter_4V); // 4-momenta
     inter_cm_4V.Boost(cm_4V);
     Vector4<double>  cm_cm_4V(cm_4V);
     cm_cm_4V.Boost(cm_4V);

     Vector4<double> pi1_cminter_4V(pi14V);
     pi1_cminter_4V.Boost(inter_4V);

     Vector4<double> eta_cminter_4V(eta4V);
     eta_cminter_4V.Boost(inter_4V);


     Vector4<double> pi1_heli_inter= getHeliVec(cm_4V, inter_4V, pi14V);
     Vector4<double> eta_heli_inter = getHeliVec(cm_4V, inter_4V, eta4V);

     double cosEta_a1=EvtDecayAngle(cm_4V, inter_4V, eta4V);


     double intensity=0.;

     std::vector<heliAmp> theAmps;

     calcHeliAmps(theAmps, motherChain, motherMs, cm_4V, eta_cm_4V+pi1_cm_4V, piRecoil_cm_4V);
     std::vector<heliAmp>::iterator heliAmpIt;


     for ( heliAmpIt=theAmps.begin(); heliAmpIt!=theAmps.end(); ++heliAmpIt)
       {
	 if (_verbose){
	   std::cout << "***Mother amp***" << std::endl;        
	   printHeliAmp((*heliAmpIt));
	 }

         double tmpIntensity=norm((*heliAmpIt).amp*BreitWigner(pi1_cm_4V+eta_cm_4V, 1.320, 0.100));

	 std::string ampString=writeHeliAmpString((*heliAmpIt));
	 std::string theHistName="costMother_"+ampString;

	 std::map<std::string, TH1F*>::iterator itHist;
         itHist=theMotherHistMap.find(theHistName);

	 TH1F* theHist=0;
          if (itHist != theMotherHistMap.end()) theHist=theMotherHistMap[theHistName];
	  else{
	    TH1F* cosThetaHist=new TH1F(theHistName.c_str(),theHistName.c_str(),200, -1., 1.);
	    theMotherHistMap[theHistName]= cosThetaHist;
	    theHist= cosThetaHist;
	  }
	 
	  double cost=(eta_cm_4V+pi1_cm_4V).CosTheta();
	  theHist->Fill(cost,tmpIntensity);
       }




     std::vector<double> daughterMs;
     double maxDaughterM=maxAbsM(motherMs);

     for(double i=-maxDaughterM;i<=maxDaughterM; i++){
       daughterMs.push_back(i);
     }

     std::vector<heliAmp> theDaughterAmps;

     calcHeliAmps(theDaughterAmps, daughterChain, daughterMs, cm_cm_4V, eta_cm_4V, pi1_cm_4V);

     for ( heliAmpIt=theDaughterAmps.begin(); heliAmpIt!=theDaughterAmps.end(); ++heliAmpIt)
       {
	 if (_verbose){
	 std::cout << "***Inter amp***" << std::endl; 
	 printHeliAmp((*heliAmpIt));
	 }

         double tmpIntensity=norm((*heliAmpIt).amp*BreitWigner(pi1_cm_4V+eta_cm_4V, 1.320, 0.100));

	 std::string ampString=writeHeliAmpString((*heliAmpIt));
	 std::string theHistName="costInter_"+ampString;

	 std::map<std::string, TH1F*>::iterator itHist;
         itHist=theDaughterHistMap.find(theHistName);

	 TH1F* theHist=0;
          if (itHist != theDaughterHistMap.end()) theHist=theDaughterHistMap[theHistName];
	  else{
	    TH1F* cosThetaHist=new TH1F(theHistName.c_str(),theHistName.c_str(),200, -1., 1.);
	    theDaughterHistMap[theHistName]= cosThetaHist;
	    theHist= cosThetaHist;
	  }
	 
	  theHist->Fill(cosEta_a1,tmpIntensity);
       }


     complex<double> theSelAmp(0.,0.);
     std::vector<heliAmp>::iterator heliAmpIt1;
     std::vector<heliAmp>::iterator heliAmpIt2;

     for ( heliAmpIt1=theAmps.begin(); heliAmpIt1!=theAmps.end(); ++heliAmpIt1)
       {

	     complex<double> theSelAmpTmp=(*heliAmpIt1).amp;
	     complex<double> theDecayAmp(0.,0.);
	     for ( heliAmpIt2=theDaughterAmps.begin(); heliAmpIt2!=theDaughterAmps.end(); ++heliAmpIt2)
	       {
		 if ((*heliAmpIt2).M==(*heliAmpIt1).lambda1)
		   {
		     theDecayAmp+=(*heliAmpIt2).amp;
		   }                 
	       }
	     theSelAmpTmp*= theDecayAmp*BreitWigner(pi1_cm_4V+eta_cm_4V, 1.320, 0.100);
	     theSelAmp+=theSelAmpTmp;
	 
     
   }

     intensity=norm(theSelAmp);

 
     Vector4<double> pi1eta4V=pi14V + eta4V;
     Vector4<double> piRecoileta4V=piRecoil4V + eta4V;
     Vector4<double> pipi4V=pi14V + piRecoil4V;
     
     dalitzHist->Fill(pi1eta4V.M()*pi1eta4V.M(),piRecoileta4V.M()*piRecoileta4V.M(),intensity);
     dalitzHist->Fill(piRecoileta4V.M()*piRecoileta4V.M(),pi1eta4V.M()*pi1eta4V.M(),intensity);
     dalitzWOHist->Fill(pi1eta4V.M()*pi1eta4V.M(),piRecoileta4V.M()*piRecoileta4V.M());
     dalitzWOHist->Fill(piRecoileta4V.M()*piRecoileta4V.M(),pi1eta4V.M()*pi1eta4V.M());
     
     dalitzAsymWOHist->Fill(pi1eta4V.M()*pi1eta4V.M(),pipi4V.M()*pipi4V.M());
     dalitzAsymWOHist->Fill(piRecoileta4V.M()*piRecoileta4V.M(),pipi4V.M()*pipi4V.M());
     
     invpietaWOHist->Fill(pi1eta4V.M());
     invpietaWOHist->Fill(piRecoileta4V.M());
     invpipiWOHist->Fill(pipi4V.M());
     invpietaHist->Fill(pi1eta4V.M(),intensity);
     invpietaHist->Fill(piRecoileta4V.M(),intensity);
     invpipiHist->Fill(pipi4V.M(),intensity);
     
     
     if ( fabs(inter_cm_4V.M() -1.320) <.05 )
       {
	   cosThetaMotherHist->Fill(inter_cm_4V.CosTheta(),intensity);
	   cosThetaMotherWOHist->Fill(inter_cm_4V.CosTheta());
	   
	   cosThetaInterHist->Fill(eta_heli_inter.CosTheta(),intensity);
	   cosThetaInterWOHist->Fill(eta_heli_inter.CosTheta());
	   // 	   cosThetaMotherHist->Fill(cosEta_a1,intensity);
	   // 	   cosThetaMotherWOHist->Fill(cosEta_a1);
       }

     counter++;          
   }

 theTFile->Write();
 theTFile->Close();
}

SeqDecHeliAmps::~SeqDecHeliAmps()
{
}

//_____________________________________________________________________________

void SeqDecHeliAmps::get4Vecs (std::ifstream& inStream, Vector4<double>& the4Vec)
{
  double tmpPx,tmpPy, tmpPz, tmpE;
  inStream >> tmpPx;
  inStream >> tmpPy;
  inStream >> tmpPz;
  inStream  >> tmpE;
  the4Vec.SetP4(tmpE, tmpPx, tmpPy, tmpPz);  
}

double SeqDecHeliAmps::EvtDecayAngle(const Vector4<double>& parent,const Vector4<double>& resonance,
                     const Vector4<double>& daughter) {

  double pd=parent*daughter;
  double pr=parent*resonance;
  double rd=resonance*daughter;
  double mp2=parent.Mass2();
  double mr2=resonance.Mass2();
  double md2=daughter.Mass2();

  double cost=(pd*mr2-pr*rd)/sqrt((pr*pr-mr2*mp2)*(rd*rd-mr2*md2));
  return cost;

}

Vector4<double> SeqDecHeliAmps::getHeliVec(const Vector4<double>& cm, const Vector4<double>& mother, const Vector4<double>& daughter){

  Vector4<double> mother_cm(mother);
  mother_cm.Boost(cm);
  
  Vector4<double> mother_cm_rot(mother_cm);
  mother_cm_rot.RotateZ(-mother_cm.Phi());
  mother_cm_rot.RotateY(-mother_cm.Theta()); 

  Vector4<double> result(daughter);
  result.Boost(cm);
  result.RotateZ(-mother_cm.Phi());
  result.RotateY(-mother_cm.Theta()); 
  
  result.Boost(mother_cm_rot);

  return result;
}



void SeqDecHeliAmps::calcHeliAmps(std::vector<heliAmp>& theAmps, const jpChain& theChain,  const std::vector<double>& Ms, const Vector4<double>& reference4V, const Vector4<double>& daughter1_4V, const Vector4<double>& daughter2_4V){

  Vector4<double> mother_4V=daughter1_4V + daughter2_4V;
  Vector4<double> daughter1_heli_4V=getHeliVec(reference4V, mother_4V, daughter1_4V);
  Vector4<double> daughter2_heli_4V=getHeliVec(reference4V, mother_4V, daughter2_4V);

  Vector4<double> mother_heli_4V=daughter1_heli_4V+daughter2_heli_4V;

  // calculate mother amp
  complex<double> theMotherAmp(0.,0.);

  vector<LS> motherLSs=GetValidLS(theChain.mother.J, theChain.mother.P, 
				  theChain.daughter1.J,theChain.daughter1.P,
				  theChain.daughter2.J, theChain.daughter2.P);

  map<Spin,map<Spin,map<Spin,complex<double> > > > Df;    
  Wigner_D(theChain.mother.J, daughter1_heli_4V.Phi(), daughter1_heli_4V.Theta(),0,Df);

  std::vector<double>::const_iterator itm;
  for (itm=Ms.begin(); itm!=Ms.end(); ++itm){

  for(Spin lam1 = -theChain.daughter1.J; lam1 <= theChain.daughter1.J; lam1++){
    for(Spin lam2 = -theChain.daughter2.J; lam2 <= theChain.daughter2.J; lam2++){
      Spin lamj = lam1 - lam2;      
       if ( fabs(lamj) > theChain.mother.J) continue;

	complex<double>  theMotherAmp(0.,0.);
	int num_ls = (int) motherLSs.size();

	
	for(int ls = 0; ls < num_ls; ls++){
	  Spin L= motherLSs[ls].L;
	  Spin S= motherLSs[ls].S;

	  theMotherAmp += sqrt(2*L+1)*Clebsch(L,0,S,lamj, theChain.mother.J,lamj)
	    *Clebsch(theChain.daughter1.J,lam1, theChain.daughter2.J,-lam2,S,lamj)*conj(Df[theChain.mother.J][(*itm)][lamj]);
	}

	heliAmp tmpHeliAmp;
	tmpHeliAmp.J=theChain.mother.J;
	tmpHeliAmp.M=(*itm);
	tmpHeliAmp.lambda1=lam1;
	tmpHeliAmp.lambda2=lam2;
	tmpHeliAmp.amp=theMotherAmp;
	theAmps.push_back(tmpHeliAmp);
      
    }
  }
  }
  
  
}

std::string SeqDecHeliAmps::writeHeliAmpString(heliAmp& theHeliAmp){
  stringstream jStrString;
  jStrString  << theHeliAmp.J;
  stringstream mStrString;
  mStrString  << theHeliAmp.M;
  stringstream lam1StrString;
  lam1StrString  << theHeliAmp.lambda1;
  stringstream lam2StrString;
  lam2StrString  << theHeliAmp.lambda2;

  std::string result="J"+jStrString.str()+"_M"+mStrString.str()+"_l"+lam1StrString.str()+"_l"+lam2StrString.str();
  return result;
}

double SeqDecHeliAmps::maxAbsM(const std::vector<double>& Ms){
  double result=0;
  std::vector<double>::const_iterator itm;
  for (itm=Ms.begin(); itm!=Ms.end(); ++itm){
    if ( result < fabs(*itm) ) result=fabs(*itm);
  }
  return result;
}

void SeqDecHeliAmps::printHeliAmp(const heliAmp& theAmp){

  std::cout << "J:" << theAmp.J
	    << "  M:" << theAmp.M 
	    << "  lam1:" << theAmp.lambda1 
	    << "  lam2:" << theAmp.lambda2
	    << "  Amp:" <<  theAmp.amp
	    << std::endl;
}

bool SeqDecHeliAmps::checkValidity(const jpChain& motherChain, const jpChain& daughterChain, std::vector<double>& motherMs){
  bool result=true;
  
  // check wheather abs(maxM) <= Jmother
  double maxM=maxAbsM(motherMs);
  std::cout <<"maxM= " << maxM << std::endl; 
  if (motherChain.mother.J < maxM){
    std::cout << "ERROR: Jmother=" << motherChain.mother.J << " < max(M)=" <<  maxM << std::endl;
    return false;  
  }

  // check LS combinations of the mother decay
  vector<LS> motherLSs=GetValidLS(motherChain.mother.J, motherChain.mother.P, 
				  motherChain.daughter1.J, motherChain.daughter1.P,
				  motherChain.daughter2.J, motherChain.daughter2.P);

  int num_ls_mother = (int) motherLSs.size();
  if(num_ls_mother==0){
    std::cout << "ERROR: No LS combination possible  for the mother decay" <<  std::endl;
    return false;
  }
  else{
    std::cout <<"The LS combinations for the mother decay are: " << std::endl;

    for(int ls = 0; ls < num_ls_mother; ls++){
      std::cout << "L= " << motherLSs[ls].L << "  S= " << motherLSs[ls].S << std::endl;
    } 
    std::cout << std::endl; 
  }

  // check LS combinations of the intermediate decay
  vector<LS> interLSs=GetValidLS(daughterChain.mother.J, daughterChain.mother.P, 
				  daughterChain.daughter1.J, daughterChain.daughter1.P,
				  daughterChain.daughter2.J, daughterChain.daughter2.P);
  int num_ls_inter = (int) interLSs.size();
  if(num_ls_inter==0){
    std::cout << "ERROR: No LS combination possible for the intermediate decay" <<  std::endl;
    return false;
  }
  else{
    std::cout <<"The ls combinations for the intermediate decay are: " << std::endl;

    for(int ls = 0; ls < num_ls_inter; ls++){
      std::cout << "l= " << interLSs[ls].L << "  s= " << interLSs[ls].S << std::endl;
    } 
    std::cout << std::endl;  
  }

  return result;
}
