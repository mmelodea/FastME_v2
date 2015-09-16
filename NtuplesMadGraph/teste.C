{
  Double_t ZZ_mass, PsbD;
  TFile *fdata = TFile::Open("/home/miqueias/Documentos/FastME/NtuplesMadGraph/ggZZ4l_plus_qqZZ4l_MadGraph.root");
  TTree *tdata = (TTree*)fdata->Get("LHE_Tree");
  tdata->SetBranchAddress("ZZ_mass",&ZZ_mass);
  int ndata = tdata->GetEntries();

  TFile *ffme = TFile::Open("/home/miqueias/Documentos/FastME/ggZZ4l_plus_qqZZ4l_MadGraphFME_MinDist.root");
  TTree *tfme = (TTree*)ffme->Get("FastME_Results");
  tfme->SetBranchAddress("PSB_Distance",&PsbD);
  
  TH1D *ZZmass = new TH1D("ZZmass","",75,50,800);
  for(int i=0; i<ndata; i++){
    tdata->GetEntry(i);
    tfme->GetEntry(i);
    
    if(PsbD < 0.5) ZZmass->Fill(ZZ_mass);
  }
  
  ZZmass->Draw();
}