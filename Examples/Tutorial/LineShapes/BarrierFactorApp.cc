//************************************************************************//
//                                                                        //
//  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)                      //
//                 Julian Pychy (julian@ep1.rub.de)                       //
//                 - Ruhr-Universität Bochum                              //
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

#include "PwaDynamics/BarrierFactor.hh"
#include "PwaDynamics/BreitWignerFunction.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"

#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"

#include <sstream>

int main(int argc, char *argv[])
{
   double decayParticle1Mass = 0.494;
   double decayParticle2Mass = 0.494;
   double motherMass = 2.049; // pbarp @ 900 MeV/c
   double recoilMass = 0.135;
   double massMin = 0.0;
   double massMax = 0.0;
   double resonanceMass = 1.68;
   double resonanceWidth = 0.150;
   double normMass=-1;
   double qR = BarrierFactor::qRDefault;
   int lprod=-1;
   int ldec=-1;

   for(int i = 0; i < argc; i++){
      if((std::string(argv[i]).compare(std::string("-h")) == 0) ||
	 (std::string(argv[i]).compare(std::string("--help")) == 0))
      {
	 std::cout << "USAGE:" << std::endl;
	 std::cout << "Definition of decay tree: mother -> recoil + X (-> dec1 + dec2)\n" << std::endl;
	 std::cout << "-dec1:\t\tmass of decay particle 1 [GeV]" << std::endl;
	 std::cout << "-dec2:\t\tmass of decay particle 2 [GeV]" << std::endl;
	 std::cout << "-mother:\tmass of (dec1, dec2)'s mother system [GeV]" << std::endl;
	 std::cout << "-recoil:\tmass of recoil system [GeV]" << std::endl;
	 std::cout << "-min:\t\tmaximum mass in plots [GeV]" << std::endl;
	 std::cout << "-max:\t\tminimum mass in plots [GeV]" << std::endl;
	 std::cout << "-resmass:\tmass of Breit-Wigner resonance [GeV]" << std::endl;
	 std::cout << "-reswidth:\twidth of Breit-Wigner resonance [GeV]" << std::endl;
         std::cout << "-qR:\t\tqR value [GeV]" << std::endl;
         std::cout << "-lprod:\t\tfix l of production [GeV]" << std::endl;
         std::cout << "-ldec:\t\tfix l of decay [GeV]" << std::endl;
         return 0;
      }
      else if(std::string(argv[i]).compare(std::string("-dec1")) == 0)
      {
         std::stringstream strStream(std::string(argv[i+1]));
         strStream >> decayParticle1Mass;
      }
      else if(std::string(argv[i]).compare(std::string("-dec2")) == 0)
      {
         std::stringstream strStream(std::string(argv[i+1]));
         strStream >> decayParticle2Mass;
      }
      else if(std::string(argv[i]).compare(std::string("-mother")) == 0)
      {
         std::stringstream strStream(std::string(argv[i+1]));
         strStream >> motherMass;
      }
      else if(std::string(argv[i]).compare(std::string("-recoil")) == 0)
      {
         std::stringstream strStream(std::string(argv[i+1]));
         strStream >> recoilMass;
      }
      else if(std::string(argv[i]).compare(std::string("-max")) == 0)
      {
         std::stringstream strStream(std::string(argv[i+1]));
         strStream >> massMax;
      }
      else if(std::string(argv[i]).compare(std::string("-min")) == 0)
      {
         std::stringstream strStream(std::string(argv[i+1]));
         strStream >> massMin;
      }
      else if(std::string(argv[i]).compare(std::string("-resmass")) == 0)
      {
         std::stringstream strStream(std::string(argv[i+1]));
         strStream >> resonanceMass;
      }
      else if(std::string(argv[i]).compare(std::string("-reswidth")) == 0)
      {
         std::stringstream strStream(std::string(argv[i+1]));
         strStream >> resonanceWidth;
      }
      else if(std::string(argv[i]).compare(std::string("-ldec")) == 0)
      {
         std::stringstream strStream(std::string(argv[i+1]));
         strStream >> ldec;
      }
      else if(std::string(argv[i]).compare(std::string("-lprod")) == 0)
      {
         std::stringstream strStream(std::string(argv[i+1]));
         strStream >> lprod;
      }
      else if(std::string(argv[i]).compare(std::string("-norm")) == 0)
      {
         std::stringstream strStream(std::string(argv[i+1]));
         strStream >> normMass;
      }
      else if(std::string(argv[i]).compare(std::string("-qR")) == 0)
      {
         std::stringstream strStream(std::string(argv[i+1]));
         strStream >> qR;
      }
   }

   if(massMin == 0){
      massMin = decayParticle1Mass + decayParticle2Mass;
   }
   if(massMax == 0){
      massMax = motherMass - recoilMass;
   }

   std::cout << "Decay particle 1 mass: " << decayParticle1Mass << std::endl;
   std::cout << "Decay particle 2 mass: " << decayParticle2Mass << std::endl;
   std::cout << "Mother mass: " << motherMass << std::endl;
   std::cout << "Recoil mass: " << recoilMass << std::endl;
   std::cout << "Mass min: " << massMin << std::endl;
   std::cout << "Mass max: " << massMax << std::endl;
   std::cout << "Resonance mass: " << resonanceMass << std::endl;
   std::cout << "Resonance width: " << resonanceWidth << std::endl;
   std::cout << "qR: " << qR << std::endl;

   TFile* file = new TFile("BarrierFactorApp.root", "recreate");
 
   int lprodMin = 0;
   int lprodMax = 6;
   int ldecMin = 0;
   int ldecMax = 6;

   if(lprod >= 0){
      lprodMin = lprodMax = lprod;
   }
   if(ldec >= 0){
      ldecMin = ldecMax = ldec;
   }

   // Barrier factors for decay
   TCanvas* decOnlyCanvas = new TCanvas("decOnlyCanvas", "decOnlyCanvas", 800, 600);
   for(int l = 6; l>=0; l--){
      std::ostringstream histName;
      histName << "BW_dec_l" << l;

      TH1F* dHist = new TH1F(histName.str().c_str(), histName.str().c_str(), 500, massMin, massMax);
      dHist->GetXaxis()->SetTitle("m [GeV]");
      dHist->SetLineColor(l+1);
      dHist->SetLineWidth(2);

      for(int i=1; i<= dHist->GetNbinsX(); i++){
	 double q = breakupMomQ(dHist->GetBinCenter(i), decayParticle1Mass, decayParticle2Mass).real();
	 std::complex<double> value = BarrierFactor::BlattWeisskopf(l, q, qR);
   	 dHist->SetBinContent(i, std::abs(value));
      }
      dHist->Scale(1./dHist->GetBinContent(dHist->GetNbinsX()));
      dHist->Draw("same");
   }
   decOnlyCanvas->Write();


   // Barrier factors for production
   TCanvas* prodOnlyCanvas = new TCanvas("prodOnlyCanvas", "prodOnlyCanvas", 800, 600);
   for(int l = 6; l>=0; l--){
      std::ostringstream histName;
      histName << "BW_prod_l" << l;

      TH1F* dHist = new TH1F(histName.str().c_str(), histName.str().c_str(), 500, massMin, massMax);
      dHist->GetXaxis()->SetTitle("m [GeV]");
      dHist->SetLineColor(l+1);
      dHist->SetLineWidth(2);

      for(int i=1; i<= dHist->GetNbinsX(); i++){
	 double q = breakupMomQ(motherMass, recoilMass, dHist->GetBinCenter(i)).real();
	 std::complex<double> value = BarrierFactor::BlattWeisskopf(l, q, qR);
   	 dHist->SetBinContent(i, std::abs(value.real()));
      }
      dHist->Scale(1./dHist->GetBinContent(1));
      dHist->Draw("same");
   }
   prodOnlyCanvas->Write();


   // Barrier factors production*decay
   TCanvas* combinedCanvas = new TCanvas("combinedCanvas", "combinedCanvas", 800, 600);
   for(int lp = lprodMax; lp>=lprodMin; lp--){
      for(int ld = ldecMax; ld>=ldecMin; ld--){
	 std::ostringstream histName;
	 histName << "BW_combined_l" << ld << "lprod" << lp;

	 TH1F* dHist = new TH1F(histName.str().c_str(), histName.str().c_str(), 500, massMin, massMax);
	 dHist->GetXaxis()->SetTitle("m [GeV]");
	 dHist->SetLineColor(lp+ld*10+1);
	 dHist->SetLineWidth(2);

	 for(int i=1; i<= dHist->GetNbinsX(); i++){
	    double qDec = breakupMomQ(dHist->GetBinCenter(i), decayParticle1Mass, decayParticle2Mass).real();
	    double qProd = breakupMomQ(motherMass, recoilMass, dHist->GetBinCenter(i)).real();
	    std::complex<double> valueDec = BarrierFactor::BlattWeisskopf(ld, qDec, qR);
  	    std::complex<double> valueProd = BarrierFactor::BlattWeisskopf(lp, qProd, qR);
	    dHist->SetBinContent(i, std::abs(valueDec.real() * valueProd.real()));
	 }

	 if(normMass > 0){
	    dHist->Scale(1./dHist->GetBinContent(dHist->FindBin(normMass)));
	 }
	 else{
	    dHist->Scale(1./dHist->GetBinContent(250));
	 }
	 dHist->Draw("same");
      }
   }
   combinedCanvas->Write();


   // Barrier factors production*decay * Breit Wigner resonance
   TCanvas* resshapeCanvas = new TCanvas("resshapeCanvas", "resshapeCanvas", 800, 600);
   for(int lp = lprodMax; lp>=lprodMin; lp--){
      for(int ld = ldecMax; ld>=ldecMin; ld--){
	 std::ostringstream histName;
	 histName << "BW_resshape_ldec" << ld << "lprod" << lp;

	 TH1F* dHist = new TH1F(histName.str().c_str(), histName.str().c_str(), 500, massMin, massMax);
	 dHist->GetXaxis()->SetTitle("m [GeV]");
	 dHist->SetLineColor(lp+ld*10+1);
	 dHist->SetLineWidth(2);

	 for(int i=1; i<= dHist->GetNbinsX(); i++){
	    double qDec = breakupMomQ(dHist->GetBinCenter(i), decayParticle1Mass, decayParticle2Mass).real();
	    double qProd = breakupMomQ(motherMass, recoilMass, dHist->GetBinCenter(i)).real();
	    std::complex<double> valueDec = BarrierFactor::BlattWeisskopf(ld, qDec, qR);
  	    std::complex<double> valueProd = BarrierFactor::BlattWeisskopf(lp, qProd, qR);
	    std::complex<double> breitWigner = BreitWignerFunction::NonRel(dHist->GetBinCenter(i), resonanceMass, resonanceWidth);
	    dHist->SetBinContent(i, std::abs(breitWigner * valueDec * valueProd));
	 }

	 if(normMass > 0){
	    dHist->Scale(1./dHist->GetBinContent(dHist->FindBin(normMass)));
	 }
	 else{
	    dHist->Scale(1./dHist->GetBinContent(dHist->FindBin(resonanceMass)));
	 }
	 dHist->Draw("same");
      }
   }
   resshapeCanvas->Write();


   file->Close();
   return 0;
}
