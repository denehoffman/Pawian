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
#include "TLegend.h"

#include <sstream>
#include <vector>

int main(int argc, char *argv[])
{
  std::vector<TH1F*>::iterator it;
  
  double decayParticle1Mass = 0.135;
   double decayParticle2Mass = 0.135;
   double motherMass = 2.049; // pbarp @ 900 MeV/c
   double recoilMass = 0.547;
   double massMin = 0.0;
   double massMax = 0.0;
   double resonanceMass = 1.26;
   double resonanceWidth = 0.150;
   double normMass=-1;
   double qRdec = BarrierFactor::qRDefault;
   double qRprod = BarrierFactor::qRDefault;
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
         std::cout << "-qRdec:\t\tqRdec value [GeV]" << std::endl;
	 std::cout << "-qRprod:\t\tqRprod value [GeV]" << std::endl;
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
      else if(std::string(argv[i]).compare(std::string("-qRdec")) == 0)
      {
         std::stringstream strStream(std::string(argv[i+1]));
         strStream >> qRdec;
      }
      else if(std::string(argv[i]).compare(std::string("-qRprod")) == 0)
      {
         std::stringstream strStream(std::string(argv[i+1]));
         strStream >> qRprod;
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
   std::cout << "qRprod: " << qRprod << std::endl;
   std::cout << "qRdec: " << qRdec << std::endl;

   TFile* file = new TFile("BarrierFactorApp.root", "recreate");
 
   int lprodMin = 0;
   int lprodMax = 5;
   int ldecMin = 0;
   int ldecMax = 2;

   if(lprod >= 0){
      lprodMin = lprodMax = lprod;
   }
   if(ldec >= 0){
      ldecMin = ldecMax = ldec;
   }

   TCanvas* allCanvas = new TCanvas("allCanvas", "allCanvas", 800, 600);
   allCanvas->Divide(3,2);

   // q values for production and decay
   TCanvas* qValCanvas = new TCanvas("qValCanvas", "qValCanvas", 800, 600);
   qValCanvas->SetObjectStat(10);
 
   TH1F* dHistqValProd = new TH1F("qValProd", "qValProd", 500, massMin, massMax);
   dHistqValProd->SetStats(0);
   dHistqValProd->GetXaxis()->SetTitle("m [GeV]");
   dHistqValProd->SetLineColor(1);
   dHistqValProd->SetLineWidth(2);

   TH1F* dHistqValDec = new TH1F("qValDec", "qValDec", 500, massMin, massMax);
   dHistqValDec->SetStats(0);
   dHistqValDec->GetXaxis()->SetTitle("m [GeV]");
    dHistqValDec->SetLineColor(2);
   dHistqValDec->SetLineWidth(2);

   for(int i=1; i<= dHistqValProd->GetNbinsX(); i++){
     double qProd = breakupMomQ(motherMass, recoilMass, dHistqValProd->GetBinCenter(i)).real();
     double qDec = breakupMomQ(dHistqValDec->GetBinCenter(i), decayParticle1Mass, decayParticle2Mass).real();
     dHistqValProd->SetBinContent(i, qProd);
     dHistqValDec->SetBinContent(i, qDec);

   }
   dHistqValProd->Draw();
   dHistqValDec->Draw("same");

   TLegend *legend_qVal = new TLegend(0.8,0.75,0.96,0.98);
   legend_qVal->AddEntry(dHistqValProd,"q prod","l");
   legend_qVal->AddEntry(dHistqValDec,"q dec","l");
   legend_qVal->Draw();
   qValCanvas->Write();

   allCanvas->cd(5);
   qValCanvas->DrawClonePad();
   complex<double> q0Dec = breakupMomQ(resonanceMass, decayParticle1Mass, decayParticle2Mass).real();

   // Barrier factors for decay
   TCanvas* decOnlyCanvas = new TCanvas("decOnlyCanvas", "decOnlyCanvas", 800, 600);
   TLegend *legend_bdec = new TLegend(0.8,0.75,0.96,0.98);
   for(int l = ldecMax; l>=ldecMin; l--){
      std::ostringstream histName;
      histName << "BW_dec_l" << l;

      TH1F* dHist = new TH1F(histName.str().c_str(), histName.str().c_str(), 500, massMin, massMax);
      dHist->SetStats(0);
      dHist->GetXaxis()->SetTitle("m [GeV]");
      dHist->SetLineColor(l+1);
      dHist->SetLineWidth(2);

      for(int i=1; i<= dHist->GetNbinsX(); i++){
	complex<double> q = breakupMomQ(dHist->GetBinCenter(i), decayParticle1Mass, decayParticle2Mass).real();
	 //	 std::complex<double> value = BarrierFactor::BlattWeisskopf(l, q, qRdec);
	 std::complex<double> value = BarrierFactor::BlattWeisskopfRatio(l, q, q0Dec, qRdec);
	 // std::cout << "BarrierFactor::BlattWeisskopf(l, q, qRdec) = " << BarrierFactor::BlattWeisskopf(l, q, qRdec) << "\t BarrierFactor::BlattWeisskopfRatio(l, q, q0Dec, qRdec) = " << BarrierFactor::BlattWeisskopfRatio(l, q, q0Dec, qRdec) << std::endl;
   	 dHist->SetBinContent(i, std::abs(value));
      }
      dHist->Scale(1./dHist->GetBinContent(dHist->GetNbinsX()));
      dHist->Draw("same");
      legend_bdec->AddEntry(dHist,histName.str().c_str(),"l");
   }
   legend_bdec->Draw();
   decOnlyCanvas->Write();

   allCanvas->cd(1);
   decOnlyCanvas->DrawClonePad();


   // Barrier factors for production
   TCanvas* prodOnlyCanvas = new TCanvas("prodOnlyCanvas", "prodOnlyCanvas", 800, 600);
   TLegend *legend_bprod = new TLegend(0.8,0.75,0.96,0.98);
   for(int l = lprodMax; l>=lprodMin; l--){
      std::ostringstream histName;
      histName << "BW_prod_l" << l;

      TH1F* dHist = new TH1F(histName.str().c_str(), histName.str().c_str(), 500, massMin, massMax);
      dHist->SetStats(0);
      dHist->GetXaxis()->SetTitle("m [GeV]");
      dHist->SetLineColor(l+1);
      dHist->SetLineWidth(2);

      for(int i=1; i<= dHist->GetNbinsX(); i++){
	 std::complex<double> q = breakupMomQ(motherMass, recoilMass, dHist->GetBinCenter(i)).real();
	 std::complex<double> qNorm=breakupMomQ(motherMass, recoilMass, decayParticle1Mass+decayParticle2Mass).real();
	 std::complex<double> value = BarrierFactor::BlattWeisskopf(l, q, qRprod)/BarrierFactor::BlattWeisskopf(l, qNorm, qRprod);
   	 dHist->SetBinContent(i, std::abs(value.real()));
      }
      dHist->Scale(1./dHist->GetBinContent(1));
      dHist->Draw("same");
      legend_bprod->AddEntry(dHist,histName.str().c_str(),"l");
   }
   legend_bprod->Draw();
   prodOnlyCanvas->Write();

   allCanvas->cd(2);
   prodOnlyCanvas->DrawClonePad();


   // Barrier factors production*decay
   TCanvas* combinedCanvas = new TCanvas("combinedCanvas", "combinedCanvas", 800, 600);
   TLegend *legend_bcomb = new TLegend(0.8,0.75,0.96,0.98);
   std::vector<TH1F*> combinedHistVec;
   double maxValCombined=0.;
   for(int lp = lprodMax; lp>=lprodMin; lp--){
      for(int ld = ldecMax; ld>=ldecMin; ld--){
	 std::ostringstream histName;
	 histName << "BW_combined_l" << ld << "lprod" << lp;

	 TH1F* dHist = new TH1F(histName.str().c_str(), histName.str().c_str(), 500, massMin, massMax);
	 dHist->SetStats(0);
	 dHist->GetXaxis()->SetTitle("m [GeV]");
	 //	 dHist->SetLineColor(lp+ld*10+1);
	 dHist->SetLineColor(lp+1);
	 dHist->SetLineWidth(2);
	 combinedHistVec.push_back(dHist);

	 for(int i=1; i<= dHist->GetNbinsX(); i++){
	   complex<double> qDec = breakupMomQ(dHist->GetBinCenter(i), decayParticle1Mass, decayParticle2Mass).real();
	   complex<double> qProd = breakupMomQ(motherMass, recoilMass, dHist->GetBinCenter(i)).real();
	   complex<double> qNorm= breakupMomQ(motherMass, recoilMass, decayParticle1Mass+decayParticle2Mass).real();
	    std::complex<double> valueDec = BarrierFactor::BlattWeisskopfRatio(ld, qDec, q0Dec, qRdec);
  	    std::complex<double> valueProd = BarrierFactor::BlattWeisskopf(lp, qProd, qRprod)/BarrierFactor::BlattWeisskopf(lp, qNorm, qRprod);
	    //	    std::cout << "valueDec: " << valueDec << "\tvalueProd: " << valueProd << "\t(valueDec* valueProd).real(): " << (valueDec*valueProd).real() << std::endl; 
	    //	    dHist->SetBinContent(i, std::abs(valueDec.real() * valueProd.real()));
	    dHist->SetBinContent(i, (valueDec*valueProd).real());
	 }

	 // if(normMass > 0){
	 //    dHist->Scale(1./dHist->GetBinContent(dHist->FindBin(normMass)));
	 // }
	 // else{
	 //    dHist->Scale(1./dHist->GetBinContent(250));
	 // }
	 double currentMax = dHist->GetMaximum();
	 if(currentMax>maxValCombined) maxValCombined=currentMax;

	 legend_bcomb->AddEntry(dHist,histName.str().c_str(),"l");
      }
   }
 
   for (it=combinedHistVec.begin(); it != combinedHistVec.end(); ++it){
     (*it)->SetMaximum(1.05 * maxValCombined);
     (*it)->Draw("same");
    }

   legend_bcomb->Draw();
   combinedCanvas->Write();

   allCanvas->cd(3);
   combinedCanvas->DrawClonePad();

   // Barrier factors production*decay * Breit Wigner resonance
   TCanvas* resshapeCanvas = new TCanvas("resshapeCanvas", "resshapeCanvas", 800, 600);
   TLegend *legend_resshape = new TLegend(0.8,0.75,0.96,0.98);

   std::vector<TH1F*> resHistVec;
   double maxValRes=0.;
   for(int lp = lprodMax; lp>=lprodMin; lp--){
      for(int ld = ldecMax; ld>=ldecMin; ld--){
	 std::ostringstream histName;
	 histName << "BW_resshape_ldec" << ld << "lprod" << lp;

	 TH1F* dHist = new TH1F(histName.str().c_str(), histName.str().c_str(), 500, massMin, massMax);
	 dHist->SetStats(0);
	 dHist->GetXaxis()->SetTitle("m [GeV]");
	 //	 dHist->SetLineColor(lp+ld*10+1);
	 dHist->SetLineColor(lp+1);
	 dHist->SetLineWidth(2);
	 resHistVec.push_back(dHist);

	 for(int i=1; i<= dHist->GetNbinsX(); i++){
	    double qProd = breakupMomQ(motherMass, recoilMass, dHist->GetBinCenter(i)).real();
	    double qNorm=breakupMomQ(motherMass, recoilMass, decayParticle1Mass+decayParticle2Mass).real();
	    std::complex<double> valueProd = BarrierFactor::BlattWeisskopf(lp, qProd, qRprod)/BarrierFactor::BlattWeisskopf(lp, qNorm, qRprod);
	    std::complex<double> breitWigner = BreitWignerFunction::BlattWRel(ld, dHist->GetBinCenter(i), resonanceMass, 
									      resonanceWidth,decayParticle1Mass, decayParticle2Mass, qRdec );

	    dHist->SetBinContent(i, std::norm(breitWigner * valueProd));
	 }

	 if(normMass > 0){
	   double bcont=dHist->GetBinContent(dHist->FindBin(normMass));
	   if(bcont > 0.) dHist->Scale(1./bcont);
	 }
	 else{
	   double bcont=dHist->GetBinContent(dHist->FindBin(resonanceMass));
	   if(bcont > 0.) dHist->Scale(1./bcont);
	 }
	 double currentMax = dHist->GetMaximum();
	 if(currentMax>maxValRes) maxValRes=currentMax;

	 legend_resshape->AddEntry(dHist,histName.str().c_str(),"l");
      }
   }

   for (it=resHistVec.begin(); it != resHistVec.end(); ++it){
     (*it)->Scale(maxValRes/(*it)->GetMaximum());
     (*it)->SetMaximum(1.05 * maxValRes);
     (*it)->Draw("same");
    }

   legend_resshape->Draw();
   resshapeCanvas->Write();

   allCanvas->cd(4);
   resshapeCanvas->DrawClonePad();
   allCanvas->Write();
   file->Write();
   file->Close();
   return 0;
}
