{
  const int N = 16;
  TString sigfiles[N] = {"ggH4l_14TeV_powheg_FMEe3MC_noScale.root",
			 "ggH4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt10.root",
			 "ggH4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt20.root",
			 "ggH4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt30.root",
			 "ggH4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt40.root",
			 "ggH4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt50.root",
			 "ggH4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt60.root",
			 "ggH4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt70.root",
			 "ggH4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt80.root",
			 "ggH4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt90.root",
			 "ggH4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt100.root",
			 "ggH4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt110.root",
			 "ggH4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt120.root",
			 "ggH4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt130.root",
			 "ggH4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt150.root",
			 "ggH4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt200.root"};
  
  TString bkgfiles[N] = {"ZZ4l_14TeV_powheg_FMEe3MC_noScale.root",
			 "ZZ4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt10.root",
			 "ZZ4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt20.root",
			 "ZZ4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt30.root",
			 "ZZ4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt40.root",
			 "ZZ4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt50.root",
			 "ZZ4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt60.root",
			 "ZZ4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt70.root",
			 "ZZ4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt80.root",
			 "ZZ4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt90.root",
			 "ZZ4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt100.root",
			 "ZZ4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt110.root",
			 "ZZ4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt120.root",
			 "ZZ4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt130.root",
			 "ZZ4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt150.root",
			 "ZZ4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt200.root"};
  
  float nSigD[N], nSigW[N];
  Double_t psbD, psbW;
  for(int i=0; i<N; i++){
    //cout<<"loop "<<i<<endl;
    TFile *file = TFile::Open(sigfiles[i]);
    TTree *tfile = (TTree*)file->Get("FastME_Results");
    tfile->SetDirectory(0);
    tfile->SetBranchAddress("PSB_Distance",&psbD);
    tfile->SetBranchAddress("PSB_Weight",&psbW);
    int nevents = tfile->GetEntries();
    
    int acpD = 0, acpW = 0;
    for(int j=0; j<nevents; j++){
      tfile->GetEntry(j);
      if(psbD > 0.5) acpD += 1;
      if(psbW > 0.5) acpW += 1;
    }  
    nSigD[i] = acpD;
    nSigW[i] = acpW;
    
    file->Close();
    //tfile->Delete();
  }
  
  
  float nBkgD[N], nBkgW[N];
  for(int i=0; i<N; i++){
    //cout<<"loop "<<i<<endl;
    TFile *file = TFile::Open(bkgfiles[i]);
    TTree *tfile = (TTree*)file->Get("FastME_Results");
    tfile->SetDirectory(0);
    tfile->SetBranchAddress("PSB_Distance",&psbD);
    tfile->SetBranchAddress("PSB_Weight",&psbW);
    int nevents = tfile->GetEntries();
    
    int acpD = 0, acpW = 0;
    for(int j=0; j<nevents; j++){
      tfile->GetEntry(j);
      if(psbD > 0.5) acpD += 1;
      if(psbW > 0.5) acpW += 1;
    }  
    nBkgD[i] = acpD;
    nBkgW[i] = acpW;
    
    file->Close();
    //tfile->Delete();
  }
  
  float vD[N], vW[N];
  for(int l=0; l<N; l++){
    vD[l] = nBkgD[l]/nSigD[l];
    vW[l] = nBkgW[l]/nSigW[l];
  }
  float x[N] = {1,10,20,30,40,50,60,70,80,90,100,110,120,130,150,200};
  
  TGraph *gD = new TGraph(N,x,vD);
  TGraph *gW = new TGraph(N,x,vW);
  gW->SetMarkerColor(kRed);
  TMultiGraph *gfinal = new TMultiGraph("gfinal","");
  gfinal->Add(gD);
  gfinal->Add(gW);
  gfinal->Draw("AP");
}