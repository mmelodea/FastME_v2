{
  TFile *fsigo = TFile::Open("Ntuples/VBF_QED4_QCD0_SIG2.root");
  TFile *fbkgo = TFile::Open("Ntuples/VBF_QED4_QCD0_BKG2.root");

  Float_t sig_mass4l, bkg_mass4l;
  TTree *tsigo = (TTree*)fsigo->Get("Higgs");
  tsigo->SetBranchAddress("H_reco_mass",&sig_mass4l);
  TTree *tbkgo = (TTree*)fbkgo->Get("Higgs");
  tbkgo->SetBranchAddress("H_reco_mass",&bkg_mass4l);
  
  
  TFile *fsig = TFile::Open("SigVBF_FME_DrOrder.root");
  TFile *fbkg = TFile::Open("BkgVBF_FME_DrOrder.root");

  Double_t sig_psb, bkg_psb;
  TTree *tsig = (TTree*)fsig->Get("FastME_Results");
  tsig->SetBranchAddress("P_SB",&sig_psb);
  TTree *tbkg = (TTree*)fbkg->Get("FastME_Results");
  tbkg->SetBranchAddress("P_SB",&bkg_psb);

  TH1D *hsig1D = new TH1D("hsig1D","P_{SB}",120,0,1); 
  hsig1D->SetLineColor(kBlue);
  hsig1D->SetLineWidth(2);
  TH1D *hbkg1D= new TH1D("hbkg1D","P_{SB}",120,0,1);
  hbkg1D->SetLineColor(kRed);
  hbkg1D->SetLineWidth(2);
  
  TH2D *hsig2D = new TH2D("hsig2D","",70,100,800,100,0,1);
  TH2D *hbkg2D = new TH2D("hbkg2D","",70,100,800,100,0,1);
  hbkg2D->SetMarkerColor(kRed);

  int n = tsigo->GetEntries();
  for(int i=0; i<n; i++){
    ///Getting the Trees buffer
    tsigo->GetEntry(i);
    tbkgo->GetEntry(i);
    tsig->GetEntry(i);
    tbkg->GetEntry(i);

    ///P_SB 1D
    hsig1D->Fill(sig_psb,1.);
    hbkg1D->Fill(bkg_psb,1.);
    
    ///P_SB vs. mass4l
    hsig2D->Fill(sig_mass4l,sig_psb,1.);
    hbkg2D->Fill(bkg_mass4l,bkg_psb,1.);
    
  }

 
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cv","",0,0,1024,768);
  //cv->Divide(2,1);
  //cv->cd(1);
  hsig1D->Draw();
  hbkg1D->Draw("same");
  //cv->cd(2);
  //hsig2D->Draw("Colz");
  //hbkg2D->Draw("Colz,same");

}
