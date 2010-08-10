{
  gROOT->Reset();
  gStyle->SetOptTitle(0);
  TCanvas *cc = new TCanvas("cc","cc",0,0,800,1200);
  cc->Divide(1,2);

  std::vector<long> generation;
  std::vector<double> evaluation0;

  std::vector<double> evaluation1;

  generation.push_back(0);
  evaluation0.push_back(1040.77);
  evaluation1.push_back(2781.61);

