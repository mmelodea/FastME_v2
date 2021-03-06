///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::		    FAST MATRIX ELEMENT			::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::	    Code Author: Miquéias M. de Almeida		::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::															::
///:: This tool takes Data, MC Signal and Background tridimensional matrix containing leptons pT, eta and phi, and      ::
///:: compute the Data-MC event distance in a phase space defined by those leptons properties. In the 4e and 4mu final  ::
///:: state there's two choices, so the tool has two ways to make the comparison between the Data-MC leptons (minimum   ::
///:: distance between leptons and the media of the possible combinations). The tool use two branches to identify those ::
///:: special final states, one for indicate the final state and one for indicate what particle is being handle.	::
///:: The final states must be identify in the following way: 0 = 4e, 1 = 4mu and 2 = 2e2mu				::
///:: 															::
///:: To run the analysis make in Root (compilation mode is the fastest):						::
///:: root -l -b -q FastME.C+												::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


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


///Scale Factors to normalize Deltas
#define scale_dPt	50.
#define scale_dEta	5.
#define scale_dPhi	pi

///Flag to control dPhi use on event distance
#define use_dPhi	false

using namespace std;


///========== Compute Distance Between Events - Minimum Distance Method ===============
Double_t ComputeDR_MinDist(Int_t FS, Int_t DataParticleID[4], Double_t Data[4][3][2],
			   Int_t McParticleID[4], Double_t MC[4][3][2]){
  
  Double_t particles_distance, min_particles_distance;
  Double_t dPt=0, dEta=0, dPhi=0;
  Double_t sum_dPt2=0, sum_dEta2=0, sum_dPhi2=0, event_distance=-1;
  
  //cout<<"\nFinal State: "<<FS<<endl;
  ///2e2mu final state
  if(FS == 2){

    for(int idt=0; idt<4; idt++)
    for(int imc=0; imc<4; imc++){
      ///Avoid different Data-MC particles comparison
      if(DataParticleID[idt] != McParticleID[imc]) continue;
      //cout<<"DataPos: "<<idt<<"  ID: "<<DataParticleID[idt]<<"  MCPos: "<<imc<<"   ID: "<<McParticleID[imc]<<endl;
  
      dPt  = (Data[idt][0][0]-MC[imc][0][0])/(scale_dPt*Data[idt][0][1]);
      dEta = (Data[idt][1][0]-MC[imc][1][0])/(scale_dEta*Data[idt][1][1]);
      if(use_dPhi == true){
      dPhi = (Data[idt][2][0]-MC[imc][2][0]);
      if(fabs(dPhi) > pi)
	dPhi = (2*pi-fabs(dPhi))/(scale_dPhi*Data[idt][2][1]);
      else
	dPhi = dPhi/(scale_dPhi*Data[idt][2][1]);
      }
      sum_dPt2  += dPt*dPt;
      sum_dEta2 += dEta*dEta;
      if(use_dPhi == true) sum_dPhi2 += dPhi*dPhi;
    }
  }
    
  ///Takes the combination that gives minimum distance in 4e and 4mu final states
  else if(FS == 0 || FS == 1){
    
    int min_imc, vmin_imc[3] = {-1, -1, -1};
    for(int idt=0; idt<4; idt++){
      min_imc = -1, particles_distance = -1; min_particles_distance = 1.E15;
      for(int imc=0; imc<4; imc++){
	///Avoid different Data-MC particles comparison
	if(DataParticleID[idt] != McParticleID[imc]) continue;
	//cout<<"DataPos: "<<idt<<"  ID: "<<DataParticleID[idt]<<"  MCPos: "<<imc<<"   ID: "<<McParticleID[imc]<<endl;

	dPt  = (Data[idt][0][0]-MC[imc][0][0])/(scale_dPt*Data[idt][0][1]);
	dEta = (Data[idt][1][0]-MC[imc][1][0])/(scale_dEta*Data[idt][1][1]);
	if(use_dPhi == true){
	dPhi = (Data[idt][2][0]-MC[imc][2][0]);
	if(fabs(dPhi) > pi)
	  dPhi = (2*pi-fabs(dPhi))/(scale_dPhi*Data[idt][2][1]);
	else
	  dPhi = dPhi/(scale_dPhi*Data[idt][2][1]);
	}
	
	particles_distance = sqrt(dPt*dPt + dEta*dEta);
	if(use_dPhi == true) particles_distance = sqrt(dPt*dPt + dEta*dEta + dPhi*dPhi);
	if(particles_distance < min_particles_distance && imc != vmin_imc[0] && imc != vmin_imc[1] && imc != vmin_imc[2]){
	  min_imc = imc;
	  min_particles_distance = particles_distance;
	}
      }
      
      ///Monitor of chosen MCs to avoid object recounting
      vmin_imc[idt] = min_imc;
      //cout<<"DataPos: "<<idt<<"  ID: "<<DataParticleID[idt]<<"  MCPos: "<<min_imc<<"   ID: "<<McParticleID[min_imc]<<endl;
      
      dPt  = (Data[idt][0][0]-MC[min_imc][0][0])/(scale_dPt*Data[idt][0][1]);
      dEta = (Data[idt][1][0]-MC[min_imc][1][0])/(scale_dEta*Data[idt][1][1]);
      if(use_dPhi == true){
      dPhi = (Data[idt][2][0]-MC[min_imc][2][0]);
      if(fabs(dPhi) > pi)
	dPhi = (2*pi-fabs(dPhi))/(scale_dPhi*Data[idt][2][1]);
      else
	dPhi = dPhi/(scale_dPhi*Data[idt][2][1]);
      }
      sum_dPt2  += dPt*dPt;
      sum_dEta2 += dEta*dEta;
      if(use_dPhi == true) sum_dPhi2 += dPhi*dPhi;
    }
  }
  
  else{
    cout<<"[Error] Final State irregular passed to ComputeDR function!"<<endl;
    throw exception();
  }

  if(use_dPhi == false) event_distance = sqrt(sum_dPt2 + sum_dEta2);
  else event_distance = sqrt(sum_dPt2 + sum_dEta2 + sum_dPhi2);
  
  if(event_distance == -1) throw exception();
  else return event_distance;
}
///================================================================================================
///================================================================================================



///========== Compute Distance Between Data-MC Events - Media Method ++============================
Double_t ComputeDR_Media(Int_t FS, Int_t DataParticleID[4], Double_t Data[4][3][2],
			 Int_t McParticleID[4], Double_t MC[4][3][2]){
  
  Double_t dPt=0, dEta=0, dPhi=0;
  Double_t sum_dPt2=0, sum_dEta2=0, sum_dPhi2=0, event_distance=-1;
  
  //cout<<"\nFS: "<<FS<<endl;
  for(int idt=0; idt<4; idt++)
  for(int imc=0; imc<4; imc++){
    ///Avoid different Data-MC particles comparison
    if(DataParticleID[idt] != McParticleID[imc]) continue;
    //cout<<"DataPos: "<<idt<<"  ID: "<<DataParticleID[idt]<<"  MCPos: "<<imc<<"   ID: "<<McParticleID[imc]<<endl;
  
    if(FS == 2){
      dPt  = (Data[idt][0][0]-MC[imc][0][0])/(scale_dPt*Data[idt][0][1]);
      dEta = (Data[idt][1][0]-MC[imc][1][0])/(scale_dEta*Data[idt][1][1]);
      if(use_dPhi == true){
	dPhi = (Data[idt][2][0]-MC[imc][2][0]);
	if(fabs(dPhi) > pi)
	  dPhi = (2*pi-fabs(dPhi))/(scale_dPhi*Data[idt][2][1]);
	else
	  dPhi = dPhi/(scale_dPhi*Data[idt][2][1]);
      }
    }
    
    ///Takes the media of combinations in 4e and 4mu final states
    else if(FS == 0 || FS == 1){
      dPt  = (Data[idt][0][0]-MC[imc][0][0])/(2*scale_dPt*Data[idt][0][1]);
      dEta = (Data[idt][1][0]-MC[imc][1][0])/(2*scale_dEta*Data[idt][1][1]);
      if(use_dPhi == true){
	dPhi = (Data[idt][2][0]-MC[imc][2][0]);
	if(fabs(dPhi) > pi)
	  dPhi = (2*pi-fabs(dPhi))/(2*scale_dPhi*Data[idt][2][1]);
	else
	  dPhi = dPhi/(2*scale_dPhi*Data[idt][2][1]);
      }
    }

    else{
      cout<<"[Error] Final State irregular passed to ComputeDR function!"<<endl;
      throw exception();
    }
    
    sum_dPt2  += dPt*dPt;
    sum_dEta2 += dEta*dEta;
    if(use_dPhi == true)
      sum_dPhi2 += dPhi*dPhi;
  }      
  
  if(use_dPhi == false) event_distance = sqrt(sum_dPt2 + sum_dEta2);
  else event_distance = sqrt(sum_dPt2 + sum_dEta2 + sum_dPhi2);

  if(event_distance == -1) throw exception();
  else return event_distance;
}
///================================================================================================



///========= Compute Discriminant Values ==========================================================
///Based on Distance
Double_t PsbD(Double_t min_dr_sig, Double_t min_dr_bkg){
  Double_t DD = min_dr_bkg/(min_dr_sig + min_dr_bkg);
  return DD;
}

///Based on Weight
Double_t PsbW(Double_t sig_event_weight, Double_t sigXS, Double_t bkg_event_weight, Double_t bkgXS){
  ///Splits the discriminant calculus according to final state (the cross section are different)
  Double_t DW = (sig_event_weight/sigXS)/(sig_event_weight/sigXS + bkg_event_weight/bkgXS);
  return DW;
}
///================================================================================================

int FastME(){
   
  ///--------------	Preparing Inputs	-----------------------------------------
  //TString Out_Name    = "ggH4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt50_Media";
  //TString Data_Path   = "NtuplesPowhegPythia/ggH4l_14TeV_powheg-pythia8.root";
  TString Out_Name    = "ZZ4l_14TeV_powheg_FMEe3MC_Scale_dEta5_dPt50_Media";
  TString Data_Path   = "NtuplesPowhegPythia/ZZTo4l_14TeV_powheg.root";
  
  //TString Out_Name      = "ggH4l_14TeV_madgraph_FMEe4MC_Scale_dEta5_dPt50";
  //TString Data_Path     = "NtuplesMadGraph/gg4l_madgraph.root";
  //TString Out_Name      = "ZZ4l_14TeV_madgraph_FMEe4MC_Scale_dEta5_dPt50";
  //TString Data_Path     = "NtuplesMadGraph/qq4l_madgraph.root";
    
  /// ZZ events including photons
  TString MC_Sig_Path = "NtuplesSherpa/10e3Ev/gg4l_14TeV_weighted.root";
  TString MC_Bkg_Path = "NtuplesSherpa/10e3Ev/qq4l_14TeV_weighted.root";
  
  TString Tree_Name   = "LHE_Tree";
  TString ID_Branch   = "ParticleID";
  TString Objs_Branch = "RecoParticle";
  TString FS_Branch   = "FinalState";
  TString EW_Branch   = "EventWeight";
  TString Num_Trials  = "Ntrials";
  ///------------------------------------------------------------------------------------
  
  ///Flag to control the Data-MC comparisons to compute DR
  ///(0: minimum distance, 1: media of combinations)
  int DrMethod = 1;
  
  ///Flags to control final state usage
  ///(to use only one channel - combinations not implemented!)
  bool fs4e 	= false;
  bool fs4u 	= false;
  bool fs2e2u 	= false;  
  
  
  ///Opening the input files
  TFile *fData = TFile::Open(Data_Path);
  TFile *fMC_Sig = TFile::Open(MC_Sig_Path);
  TFile *fMC_Bkg = TFile::Open(MC_Bkg_Path);
  TTree *Data_Tree = (TTree*)fData->Get(Tree_Name);
  TTree *MC_Sig_Tree = (TTree*)fMC_Sig->Get(Tree_Name);
  TTree *MC_Bkg_Tree = (TTree*)fMC_Bkg->Get(Tree_Name);
  
  ///Creating the Data Tree
  Int_t DataFS, DataID[4]; 
  Double_t Data[4][3][2];
  Data_Tree->SetBranchAddress(Objs_Branch,&Data);
  Data_Tree->SetBranchAddress(FS_Branch,&DataFS);
  Data_Tree->SetBranchAddress(ID_Branch,&DataID);
  int ndata = Data_Tree->GetEntries();
  
  ///Creating the Signal Tree
  Int_t MC_SIG_FS, MC_SIG_ID[4], MC_SIG_NTRIALS;
  Double_t MC_SIG_WEIGHT, MC_SIG[4][3][2];
  MC_Sig_Tree->SetBranchAddress(Objs_Branch,&MC_SIG);
  MC_Sig_Tree->SetBranchAddress(FS_Branch,&MC_SIG_FS);
  MC_Sig_Tree->SetBranchAddress(ID_Branch,&MC_SIG_ID);
  MC_Sig_Tree->SetBranchAddress(EW_Branch,&MC_SIG_WEIGHT);
  MC_Sig_Tree->SetBranchAddress(Num_Trials,&MC_SIG_NTRIALS);  
  int nsig = MC_Sig_Tree->GetEntries();
  
  ///Creating the Background Tree
  Int_t MC_BKG_FS, MC_BKG_ID[4], MC_BKG_NTRIALS;
  Double_t MC_BKG_WEIGHT, MC_BKG[4][3][2]; 
  MC_Bkg_Tree->SetBranchAddress(Objs_Branch,&MC_BKG);
  MC_Bkg_Tree->SetBranchAddress(FS_Branch,&MC_BKG_FS);
  MC_Bkg_Tree->SetBranchAddress(ID_Branch,&MC_BKG_ID);
  MC_Bkg_Tree->SetBranchAddress(EW_Branch,&MC_BKG_WEIGHT);
  MC_Bkg_Tree->SetBranchAddress(Num_Trials,&MC_BKG_NTRIALS);  
  int nbkg = MC_Bkg_Tree->GetEntries();
  
  
  cout<<"\n::::::::::::::::::::::::::::::::::::::::"<<endl;
  cout<<":::::: Starting FastME Processing ::::::"<<endl;
  cout<<"::::::::::::::::::::::::::::::::::::::::"<<endl;
  cout<<":: DrComputation Mode: ";
       if(DrMethod == 0) cout<<"Minimum Distance"<<endl;
  else if(DrMethod == 1) cout<<"Media"<<endl;
  if(use_dPhi == false)
    cout<<":: No Using dPhi!"<<endl;
  cout<<"::--------------------------------------"<<endl;
  cout<<":: #Data Events:   "<<ndata<<endl;
  cout<<":: #MC Sig Events: "<<nsig<<endl;
  cout<<":: #MC BKG Events: "<<nbkg<<endl;
  cout<<"::--------------------------------------"<<endl;
  
  ///Compute the final states yields in the ntuples
  ///And the Cross Section to compute weights
  int d4e=0, d4u=0, d2e2u=0;
  for(int fs=0; fs<ndata; fs++){
    Data_Tree->GetEntry(fs);
    if(DataFS==0) d4e    += 1;
    if(DataFS==1) d4u    += 1;
    if(DataFS==2) d2e2u  += 1;
  }
  int s4e=0, s4u=0, s2e2u=0;
  Double_t sig_4eSum_weight=0, sig_4uSum_weight=0, sig_2e2uSum_weight=0;
  Double_t sig_4eSum_trials=0, sig_4uSum_trials=0, sig_2e2uSum_trials=0;
  for(int fs=0; fs<nsig; fs++){
    MC_Sig_Tree->GetEntry(fs);
    if(MC_SIG_FS==0){
      s4e    += 1;
      sig_4eSum_weight += MC_SIG_WEIGHT;
      sig_4eSum_trials += MC_SIG_NTRIALS;
    }
    if(MC_SIG_FS==1){
      s4u    += 1;
      sig_4uSum_weight += MC_SIG_WEIGHT;
      sig_4uSum_trials += MC_SIG_NTRIALS;
    }
    if(MC_SIG_FS==2){
      s2e2u  += 1;
      sig_2e2uSum_weight += MC_SIG_WEIGHT;
      sig_2e2uSum_trials += MC_SIG_NTRIALS;
    }      
  }
  int b4e=0, b4u=0, b2e2u=0;
  Double_t bkg_4eSum_weight=0, bkg_4uSum_weight=0, bkg_2e2uSum_weight=0;
  Double_t bkg_4eSum_trials=0, bkg_4uSum_trials=0, bkg_2e2uSum_trials=0;
  for(int fs=0; fs<nbkg; fs++){
    MC_Bkg_Tree->GetEntry(fs);
    if(MC_BKG_FS==0){
      b4e    += 1;
      bkg_4eSum_weight += MC_BKG_WEIGHT;
      bkg_4eSum_trials += MC_BKG_NTRIALS;
    }
    if(MC_BKG_FS==1){
      b4u    += 1;
      bkg_4uSum_weight += MC_BKG_WEIGHT;
      bkg_4uSum_trials += MC_BKG_NTRIALS;
    }
    if(MC_BKG_FS==2){
      b2e2u  += 1;
      bkg_2e2uSum_weight += MC_BKG_WEIGHT;
      bkg_2e2uSum_trials += MC_BKG_NTRIALS;
    }      
  }
  Double_t SigTrials[3] = { sig_4eSum_trials, sig_4uSum_trials, sig_2e2uSum_trials };
  Double_t BkgTrials[3] = { bkg_4eSum_trials, bkg_4uSum_trials, bkg_2e2uSum_trials };
  
  ///Computing and storing XS
  Double_t SigXS[3] = { sig_4eSum_weight/SigTrials[0], sig_4uSum_weight/SigTrials[1], sig_2e2uSum_weight/SigTrials[2] };
  Double_t BkgXS[3] = { bkg_4eSum_weight/BkgTrials[0], bkg_4uSum_weight/BkgTrials[1], bkg_2e2uSum_weight/BkgTrials[2] };
  ///Plot on PC Screen the XS
  cout<<":: XS Computed from the Ntuples"<<endl;
  cout<<":: Sig   "<<"4e: "<<SigXS[0]<<"\t4mu: "<<SigXS[1]<<"\t2e2mu: "<<SigXS[2]<<endl;
  cout<<":: Bkg   "<<"4e: "<<BkgXS[0]<<"\t4mu: "<<BkgXS[1]<<"\t2e2mu: "<<BkgXS[2]<<endl;
  cout<<"::--------------------------------------"<<endl;
  
  cout<<":: Final State Yields in the Ntuples"<<endl;
  cout<<":: Data  "<<"4e: "<<d4e<<"\t4mu: "<<d4u<<"\t2e2mu: "<<d2e2u<<endl;
  cout<<":: Sig   "<<"4e: "<<s4e<<"\t4mu: "<<s4u<<"\t2e2mu: "<<s2e2u<<endl;
  cout<<":: Bkg   "<<"4e: "<<b4e<<"\t4mu: "<<b4u<<"\t2e2mu: "<<b2e2u<<endl;
  cout<<"::______________________________________"<<endl;
  cout<<"::--------------------------------------"<<endl;
  
  
  ///Creating Tree to store Fast Matrix Element results
  Int_t FinalState;
  Double_t dr_test, minDR_toSig, minDR_toBkg, psb_distance, sig_event_weight;
  Double_t min_dr_sig, min_dr_bkg, bkg_event_weight, psb_weight;
  TTree *FME_out = new TTree("FastME_Results","Fast Matrix Element Results");
  FME_out->SetDirectory(0);
  FME_out->Branch("FinalState",&FinalState);
  FME_out->Branch("minDR_toSig",&minDR_toSig);
  FME_out->Branch("minDR_toBkg",&minDR_toBkg);
  FME_out->Branch("PSB_Distance",&psb_distance);
  FME_out->Branch("WSig_ToEvent",&sig_event_weight);
  FME_out->Branch("WBkg_ToEvent",&bkg_event_weight);
  FME_out->Branch("PSB_Weight",&psb_weight);

  
  ///For timming the process
  time_t start, stop, delaied;
  double seconds, elapsed_time;
  string unity = "s";
  time(&start);
  ///--------------------------
  
  
  ///Starts the analysis
  if(ndata > 10050) ndata = 10050;
  //ndata = 9202;
  for(int i=0; i<ndata; i++){
    time(&delaied);
    if(i % (ndata/10) == 0 && i != 0) cout<<":: Remaining Data: "<<ndata-i<<"\t||Elapsed: "<<difftime(delaied,start)<<"s"<<endl;
    Data_Tree->GetEntry(i);
    
    //cout<<"\n\nDataFS: "<<DataFS<<endl;
	 if(fs4e   == true && DataFS != 0) continue;
    else if(fs4u   == true && DataFS != 1) continue;
    else if(fs2e2u == true && DataFS != 2) continue;
    
    ///Reseting Variables
    //event_type          = Data_EvType;
    minDR_toSig   	= pedestal;
    minDR_toBkg   	= pedestal;
    psb_distance  	= pedestal;
    sig_event_weight 	= pedestal;
    bkg_event_weight 	= pedestal;
    psb_weight		= pedestal;
    min_dr_sig    	= 1.E15;
    min_dr_bkg    	= 1.E15;
        
    
    ///:::::::::::::  Tests MC Signal  ::::::::::::::::::::::::::::::::::::::::::::::::::
    //nsig = 10000;
    for(int s=0; s<nsig; s++){
      MC_Sig_Tree->GetEntry(s);
	   if(fs4e   == true && MC_SIG_FS != 0) continue;
      else if(fs4u   == true && MC_SIG_FS != 1) continue;
      else if(fs2e2u == true && MC_SIG_FS != 2) continue;
      
      ///Checks Data-Sig Final State Compatibility
      //cout<<"SigFS: "<<MC_SIG_FS<<endl;
      if(MC_SIG_FS != DataFS) continue;
      
	   if(DrMethod == 0) dr_test = ComputeDR_MinDist(DataFS,DataID,Data,MC_SIG_ID,MC_SIG);
      else if(DrMethod == 1) dr_test = ComputeDR_Media(DataFS,DataID,Data,MC_SIG_ID,MC_SIG);
      
      if(dr_test < min_dr_sig){
	min_dr_sig = dr_test;
	sig_event_weight = MC_SIG_WEIGHT;
      }
    }
    ///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    

    ///::::::::::  Using MC Background  :::::::::::::::::::::::::::::::::::::::::::::::::
    //nbkg = 10000;
    for(int b=0; b<nbkg; b++){
      MC_Bkg_Tree->GetEntry(b);
	   if(fs4e   == true && MC_BKG_FS != 0) continue;
      else if(fs4u   == true && MC_BKG_FS != 1) continue;
      else if(fs2e2u == true && MC_BKG_FS != 2) continue;

      ///Checks Data-Sig Final State Compatibility
      if(MC_BKG_FS != DataFS) continue;
      
	   if(DrMethod == 0) dr_test = ComputeDR_MinDist(DataFS,DataID,Data,MC_BKG_ID,MC_BKG);
      else if(DrMethod == 1) dr_test = ComputeDR_Media(DataFS,DataID,Data,MC_BKG_ID,MC_BKG);
      
      if(dr_test < min_dr_bkg){
	min_dr_bkg = dr_test;
	bkg_event_weight = MC_BKG_WEIGHT;
      }
    }
    ///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    
    ///Getting the discriminant value based on event distance
    if(min_dr_sig == 0 && min_dr_bkg == 0){ 
      min_dr_sig = 1.;
      min_dr_bkg = 1.;
    }
    psb_distance = PsbD(min_dr_sig, min_dr_bkg);
    if(psb_distance < 0 || psb_distance > 1){
      cout<<"[Warning] PSB_Distance: "<<psb_distance<<"    SigD: "<<min_dr_sig<<"  BkgD: "<<min_dr_bkg<<endl;
    }
    minDR_toSig  = min_dr_sig;
    minDR_toBkg  = min_dr_bkg;
    
    
    ///Getting the discriminant value based on most close MC event weight
    if(sig_event_weight == 0 && bkg_event_weight == 0){
      sig_event_weight = 1.;
      bkg_event_weight = 1.;
    }
    psb_weight = PsbW(sig_event_weight/SigTrials[DataFS], SigXS[DataFS],
		      bkg_event_weight/BkgTrials[DataFS], BkgXS[DataFS]);
    
    FinalState = DataFS;
    ///Stores in the Tree the results
    FME_out->Fill();
  }
  
  
  TFile *FME_Results = new TFile(Out_Name+".root","recreate");
  FME_out->Write();
  FME_Results->Close();
  
  
  ///Stoping timming
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
