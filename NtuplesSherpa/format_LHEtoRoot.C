#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#define pedestal -999
#define Z_mass	 91.188

using namespace std;

void format_LHEtoRoot(ifstream &Input, TString out_name){
  

  ///Checks if file given is ok
  int ncheck = 0;
  string info;
  do{
     Input >> info;
     ncheck += 1;
     //cout<<"ncheck: "<<ncheck<<endl;
  }while(info == "" && ncheck < 100);
  if(info == ""){
    cout<<"Where's the file?!"<<endl;
    throw exception();
  }
  ///--------------------------------
  
  TLorentzVector p1, p2, p3, p4;
  
  TH1D *Zon = new TH1D("Zon","Z On-Shell",40,40,120);
  TH1D *Zoff = new TH1D("Zoff","Z Off-Shell",54,12,120);
  TH1D *ZZ = new TH1D("ZZ","ZZ Invariant Mass",195,50,2000);
  
  Int_t FS_TYPE;
  Double_t Zon_mass, Zoff_mass, ZZ_mass, EVENT_WEIGHT, RECO_PARTICLE[4][3][2];
  TTree *lheTree = new TTree("LHE_Tree","LHE Tree formated to FastME");
  lheTree->Branch("FS_TYPE",&FS_TYPE,"FS_TYPE/I");
  lheTree->Branch("EVENT_WEIGHT",&EVENT_WEIGHT,"EVENT_WEIGHT/D");
  lheTree->Branch("RECO_PARTICLE",&RECO_PARTICLE,"RECO_PARTICLE[4][3][2]/D");
  lheTree->Branch("Zon_mass",&Zon_mass,"Zon_mass/D");
  lheTree->Branch("Zoff_mass",&Zoff_mass,"Zoff_mass/D");
  lheTree->Branch("ZZ_mass",&ZZ_mass,"ZZ_mass/D");

  ///Checks the begin of event list
  do{
    Input >> info;
    if(info == "") continue;
    //cout<<"Content: "<<info<<endl;
  }while(info != "</init>");
  
  int nevents = 0, n_ele, n_mu;
  double event_weight;
  do{
      ///Counting event number
      nevents += 1;
      
      ///Reseting tree variables
      n_ele = 0, n_mu = 0;
      FS_TYPE      = pedestal;
      EVENT_WEIGHT = pedestal;
      Zon_mass     = pedestal;
      Zoff_mass    = pedestal;
      ZZ_mass      = pedestal;
      
      ///Checks event by event
      do{
	  Input >> info;
	  if(info == "") continue; ///Dumps empity spaces

	  ///Searches for event weight
	  if(info == "6"){
	    Input >> info;
	    if(info == "1"){
	     Input >> event_weight;
	     //cout<<"Event Weight: "<<event_weight<<endl;
	    }
	  }
    
	  ///Searches for electrons
	  if(info == "11"){
	    //cout<<"PDG_ID: "<<info<<endl;
	    n_ele += 1;
	    Double_t px, py, pz, e;
	    for(int i=0; i<5; i++) Input >> info; ///Dumps not requested informations
     
	    Input >> px;     
	    Input >> py;
	    Input >> pz;
	    Input >> e;
	    p1.SetPxPyPzE(px, py, pz, e);
	    //cout<<"Px: "<<px<<"  Py: "<<py<<"  Pz: "<<pz<<"  E: "<<e<<endl;
	    //cout<<"Pt: "<<p1.Pt()<<"  Eta: "<<p1.Eta()<<"  Phi: "<<p1.Phi()<<"  E: "<<p1.E()<<endl;
	  }
    
	  ///Searches for positrons
	  if(info == "-11"){
	    //cout<<"PDG_ID: "<<info<<endl;
	    n_ele += 1;
	    Double_t px, py, pz, e;
	    for(int i=0; i<5; i++) Input >> info; ///Dumps not requested informations
     
	    Input >> px;     
	    Input >> py;
	    Input >> pz;
	    Input >> e;
	    p2.SetPxPyPzE(px, py, pz, e);
	    //cout<<"Px: "<<px<<"  Py: "<<py<<"  Pz: "<<pz<<"  E: "<<e<<endl;
	    //cout<<"Pt: "<<p2.Pt()<<"  Eta: "<<p2.Eta()<<"  Phi: "<<p2.Phi()<<"  E: "<<p2.E()<<endl;
	  }

	  ///Searches for muons
	  if(info == "13"){
	    //cout<<"PDG_ID: "<<info<<endl;
	    n_mu += 1;
	    Double_t px, py, pz, e;
	    for(int i=0; i<5; i++) Input >> info; ///Dumps not requested informations
     
	    Input >> px;     
	    Input >> py;
	    Input >> pz;
	    Input >> e;
	    p3.SetPxPyPzE(px, py, pz, e);
	    //cout<<"Px: "<<px<<"  Py: "<<py<<"  Pz: "<<pz<<"  E: "<<e<<endl;
	    //cout<<"Pt: "<<p3.Pt()<<"  Eta: "<<p3.Eta()<<"  Phi: "<<p3.Phi()<<"  E: "<<p3.E()<<endl;
	  }
	
	  ///Searches for anti-muons
	  if(info == "-13"){
	    //cout<<"PDG_ID: "<<info<<endl;
	    n_mu += 1;
	    Double_t px, py, pz, e;
	    for(int i=0; i<5; i++) Input >> info; ///Dumps not requested informations
     
	    Input >> px;     
	    Input >> py;
	    Input >> pz;
	    Input >> e;
	    p4.SetPxPyPzE(px, py, pz, e);
	    //cout<<"Px: "<<px<<"  Py: "<<py<<"  Pz: "<<pz<<"  E: "<<e<<endl;
	    //cout<<"Pt: "<<p4.Pt()<<"  Eta: "<<p4.Eta()<<"  Phi: "<<p4.Phi()<<"  E: "<<p4.E()<<endl;
	  }
      }while(info != "</event>"); ///End of reading process of an event
      
    ///Checks what is the Final State
    if(n_ele == 4 && n_mu == 0) FS_TYPE = 0;
    else if(n_ele == 0 && n_mu == 4) FS_TYPE = 1;
    else if(n_ele == 2 && n_mu == 2) FS_TYPE = 2;
    else{
      cout<<"Something worng happened!"<<endl;
      throw exception();
    }
    
    ///Store the event weight
    EVENT_WEIGHT = event_weight;
    
    ///Builds the matrix for FastME analysis
    RECO_PARTICLE[0][0][0] = p1.Pt();	RECO_PARTICLE[0][0][1] = 1.;
    RECO_PARTICLE[0][1][0] = p1.Eta();	RECO_PARTICLE[0][1][1] = 1.;
    RECO_PARTICLE[0][2][0] = p1.Phi();	RECO_PARTICLE[0][2][1] = 1.;
    
    RECO_PARTICLE[1][0][0] = p2.Pt();	RECO_PARTICLE[1][0][1] = 1.;
    RECO_PARTICLE[1][1][0] = p2.Eta();	RECO_PARTICLE[1][1][1] = 1.;
    RECO_PARTICLE[1][2][0] = p2.Phi();	RECO_PARTICLE[1][2][1] = 1.;
    
    RECO_PARTICLE[2][0][0] = p3.Pt();	RECO_PARTICLE[2][0][1] = 1.;
    RECO_PARTICLE[2][1][0] = p3.Eta();	RECO_PARTICLE[2][1][1] = 1.;
    RECO_PARTICLE[2][2][0] = p3.Phi();	RECO_PARTICLE[2][2][1] = 1.;
    
    RECO_PARTICLE[3][0][0] = p4.Pt();	RECO_PARTICLE[3][0][1] = 1.;
    RECO_PARTICLE[3][1][0] = p4.Eta();	RECO_PARTICLE[3][1][1] = 1.;
    RECO_PARTICLE[3][2][0] = p4.Phi();	RECO_PARTICLE[3][2][1] = 1.;
    ///-------------------------------------------------------------
    
    Double_t mpair1 = (p1+p2).M();
    Double_t mpair2 = (p3+p4).M();
    Zon_mass = (fabs(mpair1-Z_mass)<fabs(mpair2-Z_mass))? mpair1:mpair2;
    Zoff_mass = (Zon_mass == mpair1)? mpair2:mpair1; 
    ZZ_mass = (p1+p2+p3+p4).M();
    
    ///Stores all data in Tree
    lheTree->Fill();
    
    if(Zon_mass>40 && Zoff_mass<120) Zon->Fill(Zon_mass);
    if(Zon_mass>12 && Zoff_mass<120) Zoff->Fill(Zoff_mass);
    if(ZZ_mass>50  &&  ZZ_mass<2000) ZZ->Fill(ZZ_mass);

    Input >> info;
    //cout<<"nenvents: "<<nevents<<endl;
  }while(info != "</LesHouchesEvents>");
  cout<<"Events Generated: "<<nevents<<endl;

  TFile *lheRoot = new TFile(out_name+".root","recreate");
  lheTree->Write();
  Zon->Write();
  Zoff->Write();
  ZZ->Write();
  lheRoot->Close();

}
