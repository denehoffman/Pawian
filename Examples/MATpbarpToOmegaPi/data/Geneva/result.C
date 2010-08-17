{
  gROOT->Reset();
  gStyle->SetOptTitle(0);
  TCanvas *cc = new TCanvas("cc","cc",0,0,800,1200);
  cc->Divide(1,2);

  std::vector<long> generation;
  std::vector<double> evaluation0;

  std::vector<double> evaluation1;

  generation.push_back(0);
  evaluation0.push_back(367.042);
  evaluation1.push_back(382.39);

  generation.push_back(1);
  evaluation0.push_back(-55.6183);
  evaluation1.push_back(-53.5644);

  generation.push_back(2);
  evaluation0.push_back(-81.3641);
  evaluation1.push_back(-75.7804);

  generation.push_back(3);
  evaluation0.push_back(-93.5021);
  evaluation1.push_back(-82.5565);

  generation.push_back(4);
  evaluation0.push_back(-99.2031);
  evaluation1.push_back(-96.3605);

  generation.push_back(5);
  evaluation0.push_back(-129.743);
  evaluation1.push_back(-115.763);

  generation.push_back(6);
  evaluation0.push_back(-145.032);
  evaluation1.push_back(-136.543);

  generation.push_back(7);
  evaluation0.push_back(-155.554);
  evaluation1.push_back(-155.548);

  generation.push_back(8);
  evaluation0.push_back(-194.764);
  evaluation1.push_back(-170.865);

  generation.push_back(9);
  evaluation0.push_back(-212.492);
  evaluation1.push_back(-208.085);

  generation.push_back(10);
  evaluation0.push_back(-215.349);
  evaluation1.push_back(-215.305);

  generation.push_back(11);
  evaluation0.push_back(-218.729);
  evaluation1.push_back(-216.121);

  generation.push_back(12);
  evaluation0.push_back(-225.884);
  evaluation1.push_back(-224.056);

  generation.push_back(13);
  evaluation0.push_back(-230.274);
  evaluation1.push_back(-228.606);

  generation.push_back(14);
  evaluation0.push_back(-230.285);
  evaluation1.push_back(-230.274);

  generation.push_back(15);
  evaluation0.push_back(-231.923);
  evaluation1.push_back(-230.628);

  generation.push_back(16);
  evaluation0.push_back(-232.302);
  evaluation1.push_back(-231.923);

  generation.push_back(17);
  evaluation0.push_back(-232.36);
  evaluation1.push_back(-232.322);

  generation.push_back(18);
  evaluation0.push_back(-232.385);
  evaluation1.push_back(-232.36);

  generation.push_back(19);
  evaluation0.push_back(-232.753);
  evaluation1.push_back(-232.408);

  generation.push_back(20);
  evaluation0.push_back(-232.986);
  evaluation1.push_back(-232.839);

  generation.push_back(21);
  evaluation0.push_back(-233.07);
  evaluation1.push_back(-232.986);

  generation.push_back(22);
  evaluation0.push_back(-233.09);
  evaluation1.push_back(-233.074);

  generation.push_back(23);
  evaluation0.push_back(-233.107);
  evaluation1.push_back(-233.09);

  generation.push_back(24);
  evaluation0.push_back(-233.185);
  evaluation1.push_back(-233.118);

  generation.push_back(25);
  evaluation0.push_back(-233.567);
  evaluation1.push_back(-233.203);

  generation.push_back(26);
  evaluation0.push_back(-233.574);
  evaluation1.push_back(-233.569);

  generation.push_back(27);
  evaluation0.push_back(-233.619);
  evaluation1.push_back(-233.584);

  generation.push_back(28);
  evaluation0.push_back(-233.634);
  evaluation1.push_back(-233.628);

  generation.push_back(29);
  evaluation0.push_back(-233.636);
  evaluation1.push_back(-233.635);

  generation.push_back(30);
  evaluation0.push_back(-233.636);
  evaluation1.push_back(-233.636);

  generation.push_back(31);
  evaluation0.push_back(-233.636);
  evaluation1.push_back(-233.636);

  generation.push_back(32);
  evaluation0.push_back(-233.636);
  evaluation1.push_back(-233.636);

  generation.push_back(33);
  evaluation0.push_back(-233.664);
  evaluation1.push_back(-233.636);

  generation.push_back(34);
  evaluation0.push_back(-233.669);
  evaluation1.push_back(-233.664);

  generation.push_back(35);
  evaluation0.push_back(-233.67);
  evaluation1.push_back(-233.669);

  generation.push_back(36);
  evaluation0.push_back(-233.67);
  evaluation1.push_back(-233.67);

  generation.push_back(37);
  evaluation0.push_back(-233.673);
  evaluation1.push_back(-233.67);

  generation.push_back(38);
  evaluation0.push_back(-233.673);
  evaluation1.push_back(-233.673);

  generation.push_back(39);
  evaluation0.push_back(-233.673);
  evaluation1.push_back(-233.673);

  generation.push_back(40);
  evaluation0.push_back(-233.673);
  evaluation1.push_back(-233.673);

  generation.push_back(41);
  evaluation0.push_back(-233.673);
  evaluation1.push_back(-233.673);

  generation.push_back(42);
  evaluation0.push_back(-233.673);
  evaluation1.push_back(-233.673);

  generation.push_back(43);
  evaluation0.push_back(-233.673);
  evaluation1.push_back(-233.673);

  generation.push_back(44);
  evaluation0.push_back(-233.673);
  evaluation1.push_back(-233.673);

  generation.push_back(45);
  evaluation0.push_back(-233.674);
  evaluation1.push_back(-233.673);

  generation.push_back(46);
  evaluation0.push_back(-233.674);
  evaluation1.push_back(-233.674);

  generation.push_back(47);
  evaluation0.push_back(-233.674);
  evaluation1.push_back(-233.674);

  generation.push_back(48);
  evaluation0.push_back(-233.674);
  evaluation1.push_back(-233.674);

  generation.push_back(49);
  evaluation0.push_back(-233.674);
  evaluation1.push_back(-233.674);

  generation.push_back(50);
  evaluation0.push_back(-233.674);
  evaluation1.push_back(-233.674);

  // Transfer the vectors into arrays
  double generation_arr[generation.size()];
  double evaluation0_arr[evaluation0.size()];

  for(std::size_t i=0; i<generation.size(); i++) {
     generation_arr[i] = (double)generation[i];
     evaluation0_arr[i] = evaluation0[i];
  }

  // Create a TGraph object
  TGraph *evGraph0 = new TGraph(evaluation0.size(), generation_arr, evaluation0_arr);

  double evaluation1_arr[evaluation1.size()];

  for(std::size_t i=0; i<generation.size(); i++) {
     evaluation1_arr[i] = evaluation1[i];
  }

  // Create a TGraph object
  TGraph *evGraph1 = new TGraph(evaluation1.size(), generation_arr, evaluation1_arr);

  // Do the actual drawing
  cc->cd(1);
  evGraph0->Draw("AP");
  cc->cd(2);
  evGraph1->Draw("AP");
  cc->cd();
}
