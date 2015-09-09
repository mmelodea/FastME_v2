{
  TFile *fsig = TFile::Open("ggZZ4l_MadGraph_FME_MinDist.root");
  TFile *fbkg = TFile::Open("qqZZ4l_MadGraph_FME_MinDist.root");  

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
  TH1D *sig_psbW = new TH1D("sig_psbW","",50,0,1);
  TH1D *bkg_psbD = new TH1D("bkg_psbD","",50,0,1); bkg_psbD->SetLineColor(kRed);
  TH1D *bkg_psbW = new TH1D("bkg_psbW","",50,0,1); bkg_psbW->SetLineColor(kRed);

  tsig->Project("sig_psbD","PSB_Distance");
  tsig->Project("sig_psbW","PSB_Weight");
  tbkg->Project("bkg_psbD","PSB_Distance");
  tbkg->Project("bkg_psbW","PSB_Weight");
  //TH2D *sig_comb = new TH2D("sig_comb","",50,0,1,50,0,1);
  //TH2D *bkg_comb = new TH2D("bkg_comb","",50,0,1,50,0,1);
    
  //TH2D *frac_sb_comb = new TH2D("frac_sb_comb","Fraction S/B - Distance and Weight Combined",20,0,1,20,0,1);
/*  
  for(int i=0; i<nsig; i++){
    tsig->GetEntry(i);
    tbkg->GetEntry(i);
    
    sig_comb->Fill(psbDSig,psbWSig,1);
    bkg_comb->Fill(psbDBkg,psbWBkg,1);
  }
*/  
  TCanvas *cv = new TCanvas();
  cv->Divide(2,1);
  cv->cd(1);
  sig_psbD->Draw();
  bkg_psbD->Draw("same");
  cv->cd(2);
  sig_psbW->Draw();
  bkg_psbW->Draw("same");

  //if(sig_sum == 0 && bkg_sum == 0){ sig_sum=1.; bkg_sum=1.; }
  //float ratio = (sig_sum-bkg_sum)/ntot;
  //float ratio = sig_sum/(sig_sum+bkg_sum);
  //float ratio = (sig_sum-bkg_sum)/(sig_sum+bkg_sum);
  //float dp = df+0.025, wp = wf+0.025;
  //cout<<"df: "<<dp<<"  wf: "<<wp<<"  Ssum: "<<sig_sum<<"  Bsum: "<<bkg_sum<<"  Ratio: "<<ratio<<endl;
  //frac_sb_comb->Fill(dp,wp,ratio);
  
/*
  int ntot = nsig;
  float sig_sum1, bkg_sum1, sig_sum2, bkg_sum2;
  float df;
  if(nsig != nbkg) cout<<"nSig != nBkg!!!"<<endl;
  for(int d=0; d<20; d++){
    df = d/20.;
    sig_sum1 = 0;
    bkg_sum1 = 0;
    sig_sum2 = 0;
    bkg_sum2 = 0;
  for(int i=0; i<ntot; i++){
    tsig->GetEntry(i);
    tbkg->GetEntry(i);
    
    if(psbDSig<df+0.05) sig_sum1 += 1;
    if(psbDBkg<df+0.05) bkg_sum1 += 1;
    if(psbWSig<df+0.05) sig_sum2 += 1;
    if(psbWBkg<df+0.05) bkg_sum2 += 1;
  }
  
  float ratio1 = (sig_sum1-bkg_sum1)/ntot;
  float ratio2 = (sig_sum2-bkg_sum2)/ntot;
  float dp = df+0.025;
  cout<<"df: "<<dp<<"  Ssum: "<<sig_sum2<<"  Bsum: "<<bkg_sum2<<"  Ratio: "<<ratio2<<endl;
  psbD_eff->Fill(dp,ratio1);
  psbW_eff->Fill(dp,ratio2);
}
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas();
  //cv->Divide(2,1);
  //cv->cd(1);
  //frac_sb_psbD->Draw();
  //sig_psbD_eff->Draw();
  //bkg_psbD_eff->Draw("same");
  //cv->cd(2);
  //frac_sb_psbW->Draw();
  //sig_psbW_eff->Draw();
  //bkg_psbW_eff->Draw("same");
  //cv->cd(1);
  //sig_comb->Draw("Colz");
  //cv->cd(2);
  //bkg_comb->Draw("Colz");
  //frac_sb_comb->Draw("Colz");
  psbD_eff->Draw();
  psbW_eff->Draw("same");
*/
  
/*  
  TH1D *hpsbD_sig = new TH1D("hpsbD_sig","P_{SB}(Distance)",60,-0.1,1.1);
  hpsbD_sig->SetLineColor(kRed);
  hpsbD_sig->SetLineWidth(2);
  TH1D *hpsbW_sig = new TH1D("hpsbW_sig","P_{SB}(Weight)",60,-0.1,1.1);
  hpsbW_sig->SetLineColor(kRed);
  hpsbW_sig->SetLineWidth(2);
  
  TH1D *hpsbD_bkg = new TH1D("hpsbD_bkg","P_{SB}(Distance)",60,-0.1,1.1);
  hpsbD_bkg->SetLineColor(kBlue);
  hpsbD_bkg->SetLineWidth(2);
  TH1D *hpsbW_bkg = new TH1D("hpsbW_bkg","P_{SB}(Weight)",60,-0.1,1.1);
  hpsbW_bkg->SetLineColor(kBlue);
  hpsbW_bkg->SetLineWidth(2);
 
  tsig->Project("hpsbD_sig","PSB_Distance");
  tsig->Project("hpsbW_sig","PSB_Weight");
  tbkg->Project("hpsbD_bkg","PSB_Distance");
  tbkg->Project("hpsbW_bkg","PSB_Weight");

  TLegend *leg = new TLegend(0.6,0.85,0.9,0.95);
  leg->AddEntry(hpsbD_sig,"gg #rightarrow 4l","l");
  leg->AddEntry(hpsbD_bkg,"qq #rightarrow 4l","l");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cv","",0,0,1024,768);
  cv->Divide(2,1);
  cv->cd(1);
  hpsbD_bkg->Draw();
  hpsbD_sig->Draw("same");
  leg->Draw();
  cv->cd(2);
  hpsbW_bkg->Draw();
  hpsbW_sig->Draw("same");
  leg->Draw();
*/
}
