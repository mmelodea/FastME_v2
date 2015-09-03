///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::	       FAST MATRIX ELEMENT NEW VERSION		::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::          Author: Miquéias M. de Almeida 		::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


//#include "FS4l_DrComputers.cxx"
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

#define pi		3.14159265358979312
#define pedestal 	-99					///Reset Value to Variables
#define cut 	 	0.5					///Threshold to Separate Events (Ideal Cut 0.5 - MC #Sig and #Bkg Equals)

///Sig XS
#define gg4e_XS		0.000251915
#define gg4mu_XS	0.000251915
#define gg2e2mu_XS	0.000435256

///Bkg XS
#define qq4e_XS	
#define qq4mu_XS	
#define qq2e2mu_XS	0.0521077

using namespace std;


///========== Compute Distance Between Leptons =============================================
Double_t ComputeDR(Int_t FS, Int_t DataParticleID[4], Double_t Data[4][3][2],
		 Int_t McParicleID[4], Double_t MC[4][3][2]){
  
  Double_t dPt=0, dEta=0, dPhi=0;
  Double_t sum_dPt2=0, sum_dEta2=0, sum_dPhi2=0, event_distance=-1;
  
  for(int idt=0; idt<4; idt++)
  for(int imc=0; imc<4; imc++){
    if(DataParticleID[idt] != McParicleID[imc]) continue;
  
    if(FS == 2){
      dPt  = (Data[idt][0][0]-MC[imc][0][0])/Data[idt][0][1];
      dEta = (Data[idt][1][0]-MC[imc][1][0])/Data[idt][1][1];
      dPhi = (Data[idt][2][0]-MC[imc][2][0])/Data[idt][2][1];
      if(fabs(dPhi) > pi)
	dPhi = (Data[idt][2][0]-MC[imc][2][0]-pi)/Data[idt][2][1];
    }
    else if(FS == 0 || FS == 1){
      dPt  = (Data[idt][0][0]-MC[imc][0][0])/(2*Data[idt][0][1]);
      dEta = (Data[idt][1][0]-MC[imc][1][0])/(2*Data[idt][1][1]);
      dPhi = (Data[idt][2][0]-MC[imc][2][0])/(2*Data[idt][2][1]);
      if(fabs(dPhi) > pi)
	dPhi = (Data[idt][2][0]-MC[imc][2][0]-pi)/2*Data[idt][2][1]);      
    }
    else{
      cout<<"[Error] Final State irregular passed to ComputeDR function!"<<endl;
      throw exception();
    }
    
    sum_dPt2  += dPt*dPt;
    sum_dEta2 += dEta*dEta;
    sum_dPhi2 += dPhi*dPhi;
  }      
  event_distance = sqrt(sum_dPt2 + sum_dEta2 + sum_dPhi2);
  
  
  if(event_distance == -1) throw exception();
  else return event_distance;
}
///================================================================================================


///========= Compute Discriminant Values ==========================================================
Double_t psbD(Double_t min_dr_sig, Double_t min_dr_bkg){
  return min_dr_sig/(min_dr_sig + min_dr_bkg);
}

Double_t psbW(Double_t sig_event_weight, Double_t bkg_event_weight, Int_t FS){
  
       if(FS == 0) return (sig_event_weight/gg4e_XS)/(sig_event_weight/gg4e_XS + bkg_event_weight/qq4e_XS);
  else if(FS == 1) return (sig_event_weight/gg4mu_XS)/(sig_event_weight/gg4mu_XS + bkg_event_weight/qq4mu_XS);
  else if(FS == 2) return (sig_event_weight/gg2e2mu_XS)/(sig_event_weight/gg2e2mu_XS + bkg_event_weight/qq2e2mu_XS);
  else return 1;
}
///================================================================================================

int FastME(){
 
  ///Preparing Inputs
  TString Out_Name    = "sig_plus_bkg";
  TString Data_Path   = "/home/sabayon/GitHub/FastME_v2/NtuplesMadGraph/Sig_plus_Bkg.root";
  TString MC_Sig_Path = "/home/sabayon/GitHub/FastME_v2/NtuplesSherpa/gg2e2mu_weighted.root";
  TString MC_Bkg_Path = "/home/sabayon/GitHub/FastME_v2/NtuplesSherpa/qq2e2mu_weighted.root";
  TString Tree_Name   = "LHE_Tree";
  TString Objs_Branch = "RecoParticle";
  TString FS_Branch   = "FinalState";
  TString EW_Branch   = "EventWeight";
  ///------------------------------------------------------------------------------------
  
  
  TFile *fData = TFile::Open(Data_Path);
  TFile *fMC_Sig = TFile::Open(MC_Sig_Path);
  TFile *fMC_Bkg = TFile::Open(MC_Bkg_Path);
  TTree *Data_Tree = (TTree*)fData->Get(Tree_Name);
  TTree *MC_Sig_Tree = (TTree*)fMC_Sig->Get(Tree_Name);
  TTree *MC_Bkg_Tree = (TTree*)fMC_Bkg->Get(Tree_Name);
        
  Int_t DataFS, DataID, MC_SIG_FS, MC_SIG_ID, MC_BKG_FS, MC_BKG_ID;
  Double_t MC_SIG_WEIGHT, MC_BKG_WEIGHT;
  Double_t Data[4][3][2], MC_SIG[4][3][2], MC_BKG[4][3][2];
  Data_Tree->SetBranchAddress(Objs_Branch,&Data);
  Data_Tree->SetBranchAddress(FS_Branch,&DataFS);
  Data_Tree->SetBranchAddress(ID_Branch,&DataID);
  int ndata = Data_Tree->GetEntries();
  MC_Sig_Tree->SetBranchAddress(Objs_Branch,&MC_SIG);
  MC_Sig_Tree->SetBranchAddress(FS_Branch,&MC_SIG_FS);
  MC_Sig_Tree->SetBranchAddress(ID_Branch,&MC_SIG_ID);
  MC_Sig_Tree->SetBranchAddress(EW_Branch,&MC_SIG_WEIGHT);
  int nsig = MC_Sig_Tree->GetEntries();
  MC_Bkg_Tree->SetBranchAddress(Objs_Branch,&MC_BKG);
  MC_Bkg_Tree->SetBranchAddress(FS_Branch,&MC_BKG_FS);
  MC_Bkg_Tree->SetBranchAddress(ID_Branch,&MC_BKG_ID);
  MC_Bkg_Tree->SetBranchAddress(EW_Branch,&MC_BKG_WEIGHT);
  int nbkg = MC_Bkg_Tree->GetEntries();
      
  cout<<"\n::::::::::::::::::::::::::::::::::::::::"<<endl;
  cout<<":::::: Starting FastME Processing ::::::"<<endl;
  cout<<"::::::::::::::::::::::::::::::::::::::::"<<endl;
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
  cout<<":: Data  "<<"4e: "<<d4e<<"  4mu: "<<d4u<<"  2e2mu: "<<d2e2u<<endl;
  cout<<":: Sig   "<<"4e: "<<s4e<<"  4mu: "<<s4u<<"  2e2mu: "<<s2e2u<<endl;
  cout<<":: Bkg   "<<"4e: "<<b4e<<"  4mu: "<<b4u<<"  2e2mu: "<<b2e2u<<endl;
  cout<<"::______________________________________"<<endl;
  cout<<"::--------------------------------------"<<endl;
  
  Double_t dr_test, minDR_toSig, minDR_toBkg, psb_distance, sig_event_weight, bkg_event_weight, psb_weight;
  Double_t min_dr_sig, min_dr_bkg;
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
      
      dr_test = ComputeDR(Data,MC_SIG);
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
      
      dr_test = ComputeDR(Data,MC_BKG);
      if(dr_test < min_dr_bkg){
	min_dr_bkg = dr_test;
	bkg_event_weight = MC_BKG_WEIGHT;
      }
    }
    ///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    
    ///Getting the discriminant value based on phase space distance
    if(min_dr_sig == 0 && min_dr_bkg == 0){ 
      min_dr_sig = 1.;
      min_dr_bkg = 1.;
    }
    psb_distance = psbD(min_dr_sig, min_dr_bkg);
    if(psb_distance<0 || psb_distance>1){
      cout<<"[Warning] PSB_Distance: "<<psb_distance<<"    SigD: "<<min_dr_sig<<"  BkgD: "<<min_dr_bkg<<endl;
    }
    minDR_toSig  = min_dr_sig;
    minDR_toBkg  = min_dr_bkg;
    
    
    ///Getting the discriminant value based on mc most close event weight
    if(sig_event_weight == 0 && bkg_event_weight == 0){
      sig_event_weight = 1.;
      bkg_event_weight = 1.;
    }
    psb_weight = psbW(sig_event_weight, bkg_event_weight, DataFS);
    if(psb_weight<0 || psb_weight>1){
      cout<<"[Warning] PSB_Weight: "<<psb_weight<<"    SigW: "<<sig_event_weight<<"  BkgW: "<<bkg_event_weight<<endl;
    }
          
    ///Saving Current Results
    FME_out->Fill();
    
  }
  
  TFile *FME_Results = new TFile(Out_Name+".root","recreate");
  FME_out->Write();
  FME_Results->Close();
  
  time(&stop);
  cout<<":::::::::: Process Finished ::::::::::::"<<endl;
      
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
  cout<<"::--------------------------------------\n"<<endl;
  
  
  return 0; //if process well finished
}
///====================================================================================================
