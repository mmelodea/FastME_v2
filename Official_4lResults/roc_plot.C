{
  Double_t sPsbD, sPsbW, bPsbD, bPsbW;
  TFile *sig = TFile::Open("ggH4l_14TeV_madgraph_FMEe4MC_Scale_dEta5_dPt50.root");
  TTree *tsig = (TTree*)sig->Get("FastME_Results");
  tsig->SetBranchAddress("PSB_Distance",&sPsbD);
  tsig->SetBranchAddress("PSB_Weight",&sPsbW);
  int nsig = tsig->GetEntries();
  
  TFile *bkg = TFile::Open("ZZ4l_14TeV_madgraph_FMEe4MC_Scale_dEta5_dPt50.root");
  TTree *tbkg = (TTree*)bkg->Get("FastME_Results");
  tbkg->SetBranchAddress("PSB_Distance",&bPsbD);
  tbkg->SetBranchAddress("PSB_Weight",&bPsbW);
  int nbkg = tbkg->GetEntries();
  
  if(nsig != nbkg) cout<<">>>> Different event number ["<<abs(nsig-nbkg)<<"]"<<endl;
  float FPR[100], TPR[100];
  for(int w=0; w<100; w++){
    float wc = w/100.;
    float sacp=0, bacp=0;
    for(int y=0; y<nsig; y++){
      tsig->GetEntry(y);
      //cout<<"sPsbD: "<<sPsbD;
      if(sPsbW > wc) sacp += 1;
      
      tbkg->GetEntry(y);
      //cout<<"\tbPsbD: "<<bPsbD;
      if(bPsbW > wc) bacp += 1;
    }
    float TN = abs(nbkg-bacp);
    float FP = bacp;
    float FN = abs(nsig-sacp);
    float TP = sacp;
    
    if(TP == 0 && FN == 0){ TP = 1.; FN = 1.; }
    if(TN == 0 && FP == 0){ TN = 1.; FP = 1.; }
    FPR[w] = FP/(FP + TN);
    TPR[w] = TP/(TP + FN);
    //cout<<"sacp: "<<sacp<<",\tsig_total: "<<nsig<<",\tbacp: "<<bacp<<",\tbkg_total: "<<nbkg<<endl;
    //cout<<"TN: "<<TN<<",\tFP: "<<FP<<",\tFN: "<<FN<<",\tTP: "<<TP<<endl;
  }
  TGraph *roc = new TGraph(100,FPR,TPR);
  
  TCanvas *cv = new TCanvas();
  roc->Draw("A*");
  
}