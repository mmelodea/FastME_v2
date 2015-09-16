{
  ///To plot discriminant distributions
/*  
  //Sig XS (pb)
  //double sig_sx[] = {0.000245537, 0.000245537, 0.000430708};
  //Bkg XS (pb)
  //double bkg_sx[] = {0.00666065, 0.00666065, 0.0132134};
 
  //Int_t sig_fs, bkg_fs;
  //TFile *fosig = TFile::Open("NtuplesMadGraph/ggZZ4l_MadGraph.root");
  //TTree *tosig = (TTree*)fosig->Get("LHE_Tree");
  //tosig->SetBranchAddress("FinalState",&sig_fs);

  //TFile *fobkg = TFile::Open("NtuplesMadGraph/qqZZ4l_MadGraph.root");
  //TTree *tobkg = (TTree*)fobkg->Get("LHE_Tree");
  //tobkg->SetBranchAddress("FinalState",&bkg_fs);

  TFile *fsig = TFile::Open("ggZZ4l_FME_MinDist_scale_dPt1.root");
  TFile *fbkg = TFile::Open("qqZZ4l_FME_MinDist_scale_dPt1.root");  
  //TFile *fsig = TFile::Open("Official_4lResults/ggZZ4l_FME_Media.root");
  //TFile *fbkg = TFile::Open("Official_4lResults/qqZZ4l_FME_Media.root");  

  TTree *tsig = (TTree*)fsig->Get("FastME_Results");
  TTree *tbkg = (TTree*)fbkg->Get("FastME_Results");
  
  Double_t psbDSig, psbWSig, psbDBkg, psbWBkg;
  tsig->SetBranchAddress("PSB_Distance",&psbDSig);
  tsig->SetBranchAddress("PSB_Weight",&psbWSig);
  tbkg->SetBranchAddress("PSB_Distance",&psbDBkg);
  tbkg->SetBranchAddress("PSB_Weight",&psbWBkg);
  int nsig = tsig->GetEntries();
  int nbkg = tbkg->GetEntries();
  cout<<"nSig: "<<nsig<<"  nBkg: "<<nbkg<<endl;

  TH1D *sig_psbD = new TH1D("sig_psbD","",50,0,1);
  sig_psbD->SetLineWidth(2);
  TH1D *sig_psbW = new TH1D("sig_psbW","",50,0,1);
  sig_psbW->SetLineWidth(2);
  TH1D *bkg_psbD = new TH1D("bkg_psbD","",50,0,1); 
  bkg_psbD->SetLineWidth(2);
  bkg_psbD->SetLineColor(kRed);
  bkg_psbD->GetXaxis()->SetTitle("P_{SB}(Distance)");
  bkg_psbD->GetYaxis()->SetTitle("Events/0.02");
  TH1D *bkg_psbW = new TH1D("bkg_psbW","",50,0,1);
  bkg_psbW->SetLineWidth(2);
  bkg_psbW->SetLineColor(kRed);
  bkg_psbW->GetXaxis()->SetTitle("P_{SB}(Weight)");
  bkg_psbW->GetYaxis()->SetTitle("Events/0.02");

  tsig->Project("sig_psbD","PSB_Distance");
  tsig->Project("sig_psbW","PSB_Weight");
  tbkg->Project("bkg_psbD","PSB_Distance");
  tbkg->Project("bkg_psbW","PSB_Weight");

  for(int i=0; i<nsig; i++){
    tosig->GetEntry(i);
    tsig->GetEntry(i);
    sig_psbD->Fill(psbDSig,1);
    sig_psbW->Fill(psbWSig,1);
  }
  for(int i=0; i<nbkg; i++){
    tobkg->GetEntry(i);
    tbkg->GetEntry(i);
    bkg_psbD->Fill(psbDBkg,1);
    bkg_psbW->Fill(psbWBkg,1);
  }

    
  TLegend *leg = new TLegend(0.6,0.85,0.9,0.95);
  leg->AddEntry(sig_psbD,"gg #rightarrow ZZ #rightarrow 4l","l");
  leg->AddEntry(bkg_psbD,"qq #rightarrow ZZ #rightarrow 4l","l");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);

  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas();
  cv->Divide(2,1);
  cv->cd(1);
  bkg_psbD->Draw();
  sig_psbD->Draw("same");
  leg->Draw();
  cv->cd(2);
  bkg_psbW->Draw();
  sig_psbW->Draw("same");
  leg->Draw();
*/
 
///To plot Efficiency vs Purity
 
  TFile *fdata = TFile::Open("/home/sabayon/Developing/ggZZ4l_plus_qqZZ4l_FME_MinDist_scale_dPt70.root");  
  Int_t EvType;
  Double_t psbD_Data, psbW_Data;
  TTree *tdata = (TTree*)fdata->Get("FastME_Results");
  tdata->SetBranchAddress("PSB_Distance",&psbD_Data);
  tdata->SetBranchAddress("PSB_Weight",&psbW_Data);
  tdata->SetBranchAddress("EventType",&EvType);
  int ndata = tdata->GetEntries();
  
  TH2D *eff_pur_D = new TH2D("eff_pur_D","",1000,0,1,1000,0,1);
  eff_pur_D->SetMarkerColor(kBlack);
  TH2D *eff_pur_W = new TH2D("eff_pur_W","",1000,0,1,1000,0,1);
  eff_pur_W->SetMarkerColor(kRed);

  int npoints = 100, n1=0, n2=0;
  float purityD[npoints], efficiencyD[npoints], purityW[npoints], efficiencyW[npoints], x1[npoints], x2[npoints];
  for(int cut=0; cut<npoints; cut++){
    float vcut = cut/float(npoints);
    int N_SigFME_D=0, N_SigTrue_D=0, N_SigFME_W=0, N_SigTrue_W=0;
    for(int i=0; i<ndata; i++){
      tdata->GetEntry(i);
      
      if(psbD_Data<vcut) N_SigFME_D += 1;
      if(psbW_Data>vcut) N_SigFME_W += 1;
      ///EvType == 0 is Signal and EvType == 1 is Background
      if(psbD_Data<vcut && EvType == 0) N_SigTrue_D +=1;
      if(psbW_Data>vcut && EvType == 0) N_SigTrue_W +=1;
    }
    if(N_SigFME_D != 0){
      x1[n1] = vcut;
      efficiencyD[n1] = N_SigFME_D/float(ndata);
      purityD[n1] = N_SigTrue_D/float(N_SigFME_D);
      //Histogram with efficiency vs purity showing cuts
      if(purityD[n1]>0.51 && efficiencyD[n1]>0.05) eff_pur_D->Fill(purityD[n1],efficiencyD[n1],vcut);
      n1 += 1;
    }
    if(cut%3==0 && N_SigFME_W != 0){
      x2[n2] = vcut;
      efficiencyW[n2] = N_SigFME_W/float(ndata);
      purityW[n2] = N_SigTrue_W/float(N_SigFME_W);
      eff_pur_W->Fill(purityW[n2],efficiencyW[n2],vcut);
      n2 += 1;
    }
  }

  ///Efficiency and Purity vs Discriminant Cut
/*  
  TGraph *effD_cut = new TGraph(n1,x1,efficiencyD);
  effD_cut->SetMarkerColor(kBlack);
  effD_cut->SetMarkerStyle(20);
  TGraph *purD_cut = new TGraph(n1,x1,purityD);
  purD_cut->SetMarkerColor(kRed);
  purD_cut->SetMarkerStyle(20);
  
  TGraph *effW_cut = new TGraph(n2,x2,efficiencyW);
  effW_cut->SetMarkerColor(kViolet);
  effW_cut->SetMarkerStyle(20);
  TGraph *purW_cut = new TGraph(n2,x2,purityW);
  purW_cut->SetMarkerColor(kBlue);
  purW_cut->SetMarkerStyle(20);

  TMultiGraph DistAnalysis("DistAnalysis","");
  DistAnalysis.Add(effD_cut);
  DistAnalysis.Add(purD_cut);
  DistAnalysis.Add(effW_cut);
  DistAnalysis.Add(purW_cut);
  DistAnalysis.Draw("AP");
  
  TLegend *leg = new TLegend(0.6,0.85,0.9,0.95);
  leg->AddEntry(effD_cut,"P_{SB}(Distance) Efficiency","p");
  leg->AddEntry(purD_cut,"P_{SB}(Distance) Purity","p");
  leg->AddEntry(effW_cut,"P_{SB}(Weight) Efficiency","p");
  leg->AddEntry(purW_cut,"P_{SB}(Weight) Purity","p");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->Draw();
*/

 ///Efficiency vs Purity without show cuts
 /*
 TGraph *effD_vs_purD = new TGraph(n1,purityD,efficiencyD);
 effD_vs_purD->SetMarkerColor(kBlack);
 TGraph *effW_vs_purW = new TGraph(n2,purityW,efficiencyW);
 effW_vs_purW->SetMarkerColor(kRed);
 
 TMultiGraph DistAnalysis("DistAnalysis","");
 DistAnalysis.Add(effD_vs_purD);
 DistAnalysis.Add(effW_vs_purW);
 DistAnalysis.Draw("AP");
 */ 

 ///With cuts
 eff_pur_D->Draw("Text");
 eff_pur_W->Draw("Text,same");
 
 
 ///Plot Discriminant vs ZZ mass
/*
 //Sig XS (pb)
 double sig_sx[] = {1,1,1};//0.000245537, 0.000245537, 0.000430708};
 //Bkg XS (pb)
 double bkg_sx[] = {1,1,1};//0.00666065, 0.00666065, 0.0132134};
 
 Int_t sig_fs, bkg_fs;
 Double_t sig_ZZ_mass, sig_PsbD, sig_PsbW, bkg_ZZ_mass, bkg_PsbD, bkg_PsbW;
 TFile *fsig = TFile::Open("NtuplesMadGraph/ggZZ4l_MadGraph.root");
 TTree *tsig = (TTree*)fsig->Get("LHE_Tree");
 tsig->SetBranchAddress("ZZ_mass",&sig_ZZ_mass);
 tsig->SetBranchAddress("FinalState",&sig_fs);
 int nsig = tsig->GetEntries();

 TFile *sig_fme = TFile::Open("dPtScaleTests/ggZZ4l_FME_MinDist_scale_dPt70.root");
 TTree *tsig_fme = (TTree*)sig_fme->Get("FastME_Results");
 tsig_fme->SetBranchAddress("PSB_Distance",&sig_PsbD);
 tsig_fme->SetBranchAddress("PSB_Weight",&sig_PsbW);

 
 TFile *fbkg = TFile::Open("NtuplesMadGraph/qqZZ4l_MadGraph.root");
 TTree *tbkg = (TTree*)fbkg->Get("LHE_Tree");
 tbkg->SetBranchAddress("ZZ_mass",&bkg_ZZ_mass);
 tbkg->SetBranchAddress("FinalState",&bkg_fs);
 int nbkg = tbkg->GetEntries();

 TFile *bkg_fme = TFile::Open("dPtScaleTests/qqZZ4l_FME_MinDist_scale_dPt70.root");
 TTree *tbkg_fme = (TTree*)bkg_fme->Get("FastME_Results");
 tbkg_fme->SetBranchAddress("PSB_Distance",&bkg_PsbD);
 tbkg_fme->SetBranchAddress("PSB_Weight",&bkg_PsbW);
 
  
  TH2D *sig_psbD_ZZmass = new TH2D("sig_psbD_ZZmass","",75,50,800,50,0,1);
  sig_psbD_ZZmass->SetMarkerColor(kBlue);
  sig_psbD_ZZmass->GetXaxis()->SetTitle("ZZ Mass (GeV)");
  sig_psbD_ZZmass->GetYaxis()->SetTitle("P_{SB}(Distance)");
  sig_psbD_ZZmass->SetTitle("gg #rightarrow ZZ #rightarrow 4l");
  
  TH2D *sig_psbW_ZZmass = new TH2D("sig_psbW_ZZmass","",75,50,800,50,0,1);
  sig_psbW_ZZmass->SetMarkerColor(kBlue);
  sig_psbW_ZZmass->GetXaxis()->SetTitle("ZZ Mass (GeV)");
  sig_psbW_ZZmass->GetYaxis()->SetTitle("P_{SB}(Weight)");
  sig_psbW_ZZmass->SetTitle("gg #rightarrow ZZ #rightarrow 4l");
  
  TH2D *bkg_psbD_ZZmass = new TH2D("bkg_psbD_ZZmass","",75,50,800,50,0,1);
  bkg_psbD_ZZmass->SetMarkerColor(kRed);
  bkg_psbD_ZZmass->GetXaxis()->SetTitle("ZZ Mass (GeV)");
  bkg_psbD_ZZmass->GetYaxis()->SetTitle("P_{SB}(Distance)");
  bkg_psbD_ZZmass->SetTitle("qq #rightarrow ZZ #rightarrow 4l");

  TH2D *bkg_psbW_ZZmass = new TH2D("bkg_psbW_ZZmass","",75,50,800,50,0,1);
  bkg_psbW_ZZmass->SetMarkerColor(kRed);
  bkg_psbW_ZZmass->GetXaxis()->SetTitle("ZZ Mass (GeV)");
  bkg_psbW_ZZmass->GetYaxis()->SetTitle("P_{SB}(Weight)");
  bkg_psbW_ZZmass->SetTitle("qq #rightarrow ZZ #rightarrow 4l");

  
  for(int i=0; i<nsig; i++){
    tsig->GetEntry(i);
    tsig_fme->GetEntry(i);
    sig_psbD_ZZmass->Fill(sig_ZZ_mass,sig_PsbD,sig_sx[sig_fs]);
    sig_psbW_ZZmass->Fill(sig_ZZ_mass,sig_PsbW,sig_sx[sig_fs]);
  }
  for(int i=0; i<nbkg; i++){
    tbkg->GetEntry(i);
    tbkg_fme->GetEntry(i);
    bkg_psbD_ZZmass->Fill(bkg_ZZ_mass,bkg_PsbD,bkg_sx[bkg_fs]);
    bkg_psbW_ZZmass->Fill(bkg_ZZ_mass,bkg_PsbW,bkg_sx[bkg_fs]);    
  }
  
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas();
  cv->Divide(2,2);
  cv->cd(1);
  sig_psbD_ZZmass->Draw("Colz");
  cv->cd(2);
  sig_psbW_ZZmass->Draw("Colz");
  cv->cd(3);
  bkg_psbD_ZZmass->Draw("Colz");
  cv->cd(4);
  bkg_psbW_ZZmass->Draw("Colz");
*/
}
