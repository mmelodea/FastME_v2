///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::	       FAST MATRIX ELEMENT NEW VERSION		::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::          Author: Miquéias M. de Almeida 		::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#include "FS4l_DrComputers.cxx"
#include <iostream>
#include <string>
#include <ctime>
#include <exception>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TCanvas.h>

#define pedestal 	-99					///Reset Value to Variables
#define cut 	 	0.5					///Threshold to Separate Events (Ideal Cut 0.5 - MC #Sig and #Bkg Equals)

using namespace std;


int FastME(){
 
  
  ///Preparing Inputs
  TString Out_Name    = "ggZZ2e2mu_FME_DrMedio";
  TString Model       = "DrMedio";
  TString Data_Path   = "/home/sabayon/GitHub/FastME_v2/NtuplesSherpa/ggZZ2e2mu_unweighted.root";
  TString MC_Sig_Path = "/home/sabayon/GitHub/FastME_v2/NtuplesSherpa/ggH2e2mu_weighted.root";
  TString MC_Bkg_Path = "/home/sabayon/GitHub/FastME_v2/NtuplesSherpa/ggZZ2e2mu_weighted.root";
  TString Tree_Name   = "LHE_Tree";
  TString Objs_Branch = "RECO_PARTICLE";
  TString FS_Branch   = "FS_TYPE";
  TString EW_Branch   = "EVENT_WEIGHT";
  ///------------------------------------------------------------------------------------
  
  
  TFile *fData = TFile::Open(Data_Path);
  TFile *fMC_Sig = TFile::Open(MC_Sig_Path);
  TFile *fMC_Bkg = TFile::Open(MC_Bkg_Path);
  TTree *Data_Tree = (TTree*)fData->Get(Tree_Name);
  TTree *MC_Sig_Tree = (TTree*)fMC_Sig->Get(Tree_Name);
  TTree *MC_Bkg_Tree = (TTree*)fMC_Bkg->Get(Tree_Name);
        
  Int_t DataFS, MC_SIG_FS, MC_BKG_FS;
  Double_t MC_SIG_WEIGHT, MC_BKG_WEIGHT;
  Double_t Data[4][3][2], MC_SIG[4][3][2], MC_BKG[4][3][2];
  Data_Tree->SetBranchAddress(Objs_Branch,&Data);
  Data_Tree->SetBranchAddress(FS_Branch,&DataFS);
  int ndata = Data_Tree->GetEntries();
  MC_Sig_Tree->SetBranchAddress(Objs_Branch,&MC_SIG);
  MC_Sig_Tree->SetBranchAddress(FS_Branch,&MC_SIG_FS);
  MC_Sig_Tree->SetBranchAddress(EW_Branch,&MC_SIG_WEIGHT);
  int nsig = MC_Sig_Tree->GetEntries();
  MC_Bkg_Tree->SetBranchAddress(Objs_Branch,&MC_BKG);
  MC_Bkg_Tree->SetBranchAddress(FS_Branch,&MC_BKG_FS);
  MC_Bkg_Tree->SetBranchAddress(EW_Branch,&MC_BKG_WEIGHT);
  int nbkg = MC_Bkg_Tree->GetEntries();
      
  cout<<"\n::::::::::::::::::::::::::::::::::::::::"<<endl;
  cout<<":::::: Starting FastME Processing ::::::"<<endl;
  cout<<"::::::::::::::::::::::::::::::::::::::::"<<endl;
  cout<<":: Final State:    "<<"4l"<<endl;
  cout<<":: Model Chosen:   "<<Model<<endl;
  cout<<"::--------------------------------------"<<endl;
  cout<<":: #Data Events:   "<<ndata<<endl;
  cout<<":: #MC Sig Events: "<<nsig<<endl;
  cout<<":: #MC BKG Events: "<<nbkg<<endl;
  cout<<"::--------------------------------------"<<endl;
  
  ///Compute the final state yields in the ntuples
  int d4e=0, d4u=0, d2e2u=0;
  for(int fs=0; fs<ndata; fs++){
    Data_Tree->GetEntry(fs);
    if(DataFS==0) d4e    += 1;
    if(DataFS==1) d4u    += 1;
    if(DataFS==2) d2e2u  += 1;
  }
  int s4e=0, s4u=0, s2e2u=0;
  for(int fs=0; fs<nsig; fs++){
    MC_Sig_Tree->GetEntry(fs);
    if(MC_SIG_FS==0) s4e    += 1;
    if(MC_SIG_FS==1) s4u    += 1;
    if(MC_SIG_FS==2) s2e2u  += 1;
  }
  int b4e=0, b4u=0, b2e2u=0;
  for(int fs=0; fs<nbkg; fs++){
    MC_Bkg_Tree->GetEntry(fs);
    if(MC_BKG_FS==0) b4e    += 1;
    if(MC_BKG_FS==1) b4u    += 1;
    if(MC_BKG_FS==2) b2e2u  += 1;
  }
  cout<<":: Final State Yields in the Ntuples"<<endl;
  cout<<":: Ntuple   "<<"4e     "<<"4mu     "<<"2e2mu"<<endl;
  cout<<":: Data     "<<d4e<<"   "<<d4u<<"    "<<d2e2u<<endl;
  cout<<":: Sig      "<<s4e<<"   "<<s4u<<"    "<<s2e2u<<endl;
  cout<<":: Bkg      "<<b4e<<"   "<<b4u<<"    "<<b2e2u<<endl;
  cout<<"::______________________________________"<<endl;
  cout<<"::--------------------------------------"<<endl;
  
  Double_t dr_test, minDR_toSig, minDR_toBkg, psb_distance, sig_event_weight, bkg_event_weight, psb_weight;
  Double_t sig_frac1, bkg_frac1, sig_frac2, bkg_frac2, min_dr_sig, min_dr_bkg;
  int signal1 = 0, background1 = 0, signal2 = 0, background2 = 0;
  TTree *FME_out = new TTree("FastME_Results","Fast Matrix Element Results");
  FME_out->SetDirectory(0);
  FME_out->Branch("minDR_toSig",&minDR_toSig);
  FME_out->Branch("minDR_toBkg",&minDR_toBkg);
  FME_out->Branch("PSB_Distance",&psb_distance);
  FME_out->Branch("WSig_ToEvent",&sig_event_weight);
  FME_out->Branch("WBkg_ToEvent",&bkg_event_weight);
  FME_out->Branch("PSB_Weight",&psb_weight);

  ///For timming the process
  time_t start, stop;
  double seconds, elapsed_time;
  string unity = "s";
  time(&start);
  ///--------------------------
  
  for(int i=0; i<ndata; i++){
    if(i % (ndata/10) == 0 && i != 0) cout<<":: Remaining Data: "<<ndata-i<<endl;
    Data_Tree->GetEntry(i);
    
    ///Reseting Variables
    minDR_toSig   	= pedestal;
    minDR_toBkg   	= pedestal;
    psb_distance  	= pedestal;
    sig_event_weight 	= pedestal;
    bkg_event_weight 	= pedestal;
    psb_weight		= pedestal;
    min_dr_sig    	= 1.E15;
    min_dr_bkg    	= 1.E15;
        
    
    ///:::::::::::::  Using MC Sig  :::::::::::::::::::::::::::::::::::::::::::::::::::::
    for(int s=0; s<nsig; s++){
      MC_Sig_Tree->GetEntry(s);
      
      ///Checks Data-Sig Final State Compatibility
      if(DataFS != MC_SIG_FS) continue;
      
      dr_test = FS4l_DrComputers(Data,MC_SIG,Model);
      if(dr_test < min_dr_sig){
	min_dr_sig = dr_test;
	sig_event_weight = MC_SIG_WEIGHT;
      }
    }
    ///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    

    ///::::::::::  Using MC Bkg  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    for(int b=0; b<nbkg; b++){
      MC_Bkg_Tree->GetEntry(b);

      ///Checks Data-Sig Final State Compatibility
      if(DataFS != MC_BKG_FS) continue;
      
      dr_test = FS4l_DrComputers(Data,MC_BKG,Model);
      if(dr_test < min_dr_bkg){
	min_dr_bkg = dr_test;
	bkg_event_weight = MC_BKG_WEIGHT;
      }
    }    
    ///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    
    ///Getting the discriminant value based on phase space distance
    if(min_dr_sig == 0 && min_dr_bkg == 0){ 
      min_dr_sig = 1;
      min_dr_bkg = 1;
    }
    psb_distance = min_dr_sig/(min_dr_sig + min_dr_bkg);
    minDR_toSig  = min_dr_sig;
    minDR_toBkg  = min_dr_bkg;
    
    
    ///Getting data event weight
    if(sig_event_weight == 0 && bkg_event_weight == 0){
      sig_event_weight = 1.;
      bkg_event_weight = 1.;
    }
    psb_weight = sig_event_weight/(sig_event_weight + bkg_event_weight);
    
    ///Checking amount of kind choices
    if(psb_distance < cut) signal1 += 1;
    if(psb_distance > cut) background1 += 1;
    if(i == ndata-1){
      sig_frac1 = (signal1/float(ndata))*100;
      bkg_frac1 = (background1/float(ndata))*100;
    }    
    if(psb_weight < cut) signal2 += 1;
    if(psb_weight > cut) background2 += 1;
    if(i == ndata-1){
      sig_frac2 = (signal2/float(ndata))*100;
      bkg_frac2 = (background2/float(ndata))*100;
    }    
      
    ///Saving Current Results
    FME_out->Fill();
    
  }
  
  TFile *FME_Results = new TFile(Out_Name+".root","recreate");
  FME_out->Write();
  FME_Results->Close();
  
  time(&stop);
  cout<<":::::::::: Process Finished ::::::::::::"<<endl;
  cout<<"::------- Using PSB_Distance -----------"<<endl;
  cout<<":: SigFrac =      "<<sig_frac1<<"%"<<endl;
  cout<<":: BkgFrac =      "<<bkg_frac1<<"%"<<endl;
  cout<<"::-------- Using PSB_Weight -----------"<<endl;
  cout<<":: SigFrac =      "<<sig_frac2<<"%"<<endl;
  cout<<":: BkgFrac =      "<<bkg_frac2<<"%"<<endl;

      
  seconds = difftime(stop,start);
  if(seconds < 60) elapsed_time = seconds;
  if(seconds >= 60){
    elapsed_time = seconds/60.;
    unity = "min.";
  }
  if(seconds >= 3600){
    elapsed_time = seconds/3600.;
    unity = "h";
  }
  cout<<":: Time Consumed: "<<elapsed_time<<unity<<endl;
  cout<<"::--------------------------------------"<<endl;
  
  
  return 0; //if process well finished
}
///====================================================================================================
