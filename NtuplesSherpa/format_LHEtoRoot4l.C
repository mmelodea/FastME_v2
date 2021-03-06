#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <exception>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#define pedestal -999
#define Z_mass	 91.188

using namespace std;

void format_LHEtoRoot4l(){
  TString out_name;
  
  ///Signal
  //ifstream Input("/home/sabayon/SHERPA_v2_2_0/work_area/FastME_OfficialNtuples/10e4Ev/gg4e/gg4e_14TeV_weighted.lhe");
  //out_name = "gg4e_14TeV_weighted";
  //ifstream Input("/home/sabayon/SHERPA_v2_2_0/work_area/FastME_OfficialNtuples/10e4Ev/gg4mu/gg4mu_14TeV_weighted.lhe");
  //out_name = "gg4mu_14TeV_weighted";
  //ifstream Input("/home/sabayon/SHERPA_v2_2_0/work_area/FastME_OfficialNtuples/10e4Ev/gg2e2mu/gg2e2mu_14TeV_weighted.lhe");
  //out_name = "gg2e2mu_14TeV_weighted";

  ///Background
  //ifstream Input("/home/sabayon/SHERPA_v2_2_0/work_area/FastME_OfficialNtuples/10e4Ev/qq4e/qq4e_14TeV_weighted.lhe");
  //out_name = "qq4e_14TeV_weighted";
  //ifstream Input("/home/sabayon/SHERPA_v2_2_0/work_area/FastME_OfficialNtuples/10e4Ev/qq4mu/qq4mu_14TeV_weighted.lhe");
  //out_name = "qq4mu_14TeV_weighted";
  ifstream Input("/home/sabayon/SHERPA_v2_2_0/work_area/FastME_OfficialNtuples/10e4Ev/qq2e2mu/qq2e2mu_14TeV_weighted.lhe");
  out_name = "qq2e2mu_14TeV_weighted";

  string status;

  TH1D *Zon = new TH1D("Zon","Z On-Shell",40,40,120);
  TH1D *Zoff = new TH1D("Zoff","Z Off-Shell",54,12,120);
  TH1D *ZZ = new TH1D("ZZ","ZZ Invariant Mass",195,50,2000);  
  
  Int_t ParticleID[4], FinalState, EventType = 0, Ntrials;
  Double_t Zon_mass, Zoff_mass, ZZ_mass, EventWeight, RecoParticle[4][3][2];
  TTree *lheTree = new TTree("LHE_Tree","LHE Tree formated to FastME");
  lheTree->Branch("EventWeight",&EventWeight,"EventWeight/D");
  lheTree->Branch("Ntrials",&Ntrials,"Ntrials/I");
  lheTree->Branch("ParticleID",&ParticleID,"ParticleID[4]/I");
  lheTree->Branch("FinalState",&FinalState,"FinalState/I");
  lheTree->Branch("RecoParticle",&RecoParticle,"RecoParticle[4][3][2]/D");
  //lheTree->Branch("EventType",&EventType,"EventType/I");
  //lheTree->Branch("Zon_mass",&Zon_mass,"Zon_mass/D");
  //lheTree->Branch("Zoff_mass",&Zoff_mass,"Zoff_mass/D");
  //lheTree->Branch("ZZ_mass",&ZZ_mass,"ZZ_mass/D");

  std::vector<Int_t> partID;
  std::vector<TLorentzVector> part4p;
  TLorentzVector reset4p(-1,0,0,0), fv_tmp;
  
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
    
  ///Checks the begin of event list
  do{
    Input >> info;
    if(info == "") continue;
    //cout<<"Content: "<<info<<endl;
  }while(info != "</init>");
  
  int nevents = 0, trials, ntrials=0, n_ele, n_mu, negative;
  Double_t event_weight, Total_XS=0;
  do{
      ///Counting event number
      nevents += 1;
      
      ///Reseting tree variables
      fv_tmp	  = reset4p;
      EventWeight = pedestal;
      Zon_mass    = pedestal;
      Zoff_mass   = pedestal;
      ZZ_mass     = pedestal;
      
      ///Checks event by event
      n_ele = 0; n_mu = 0;
      do{
	  status = "";
	  Input >> info;
	  if(info == "") continue; ///Dumps empity spaces
	  //cout<<"Content: "<<info<<endl;
    
	  ///Searches for trials to weight event
	  if(info == "trials="){
	    Input >> trials;
	    ntrials += trials;
	    //cout<<"Content: "<<trials<<endl;
	  }
	  ///Searches for event weight
	  if(info == "6"){
	    Input >> info;
	    if(info == "1"){
	     Input >> event_weight;
	     //cout<<"Event Weight: "<<event_weight<<endl;
	     ///Computes the total Cross Section of process - have to fix.. the value not matches with given by Sherpa
	     Total_XS += event_weight;
	     //cout<<"Currently XS: "<<Total_XS<<endl;
	    }
	  }
          //if(event_weight<0) cout<<"Strange! Event["<<nevents<<"] negative weight!"<<endl;
	  
	  ///Searches for electrons/ positrons
	  if(info == "11" || info == "-11"){
   	    Input >> status;
	    if(status != "1") continue; ///Checks if particle is outgoing final state

	    n_ele += 1;
	       if(info ==  "11") partID.push_back(11);
	    else (info == "-11") partID.push_back(-11);
	    
	    Double_t px, py, pz, e;
	    for(int i=0; i<4; i++) Input >> info; ///Dumps not requested informations
     
	    Input >> px;     
	    Input >> py;
	    Input >> pz;
	    Input >> e;
	    fv_tmp.SetPxPyPzE(px, py, pz, e);
	    part4p.push_back(fv_tmp);	    
	    //if(nevents == 1) cout<<"Px: "<<px<<"  Py: "<<py<<"  Pz: "<<pz<<"  E: "<<e<<endl;
	  }
    
	  ///Searches for muons/ anti-muons
	  if(info == "13" || info == "-13"){
	    Input >> status;
	    if(status != "1") continue; ///Checks if particle is outgoing final state

	    n_mu += 1;
	       if(info ==  "13") partID.push_back(13);
	    else (info == "-13") partID.push_back(-13);

	    Double_t px, py, pz, e;
	    for(int i=0; i<4; i++) Input >> info; ///Dumps not requested informations
     
	    Input >> px;     
	    Input >> py;
	    Input >> pz;
	    Input >> e;
	    fv_tmp.SetPxPyPzE(px, py, pz, e);
	    part4p.push_back(fv_tmp);	    
	    //if(nevents == 1) cout<<"Px: "<<px<<"  Py: "<<py<<"  Pz: "<<pz<<"  E: "<<e<<endl;
	  }

      }while(info != "</event>"); ///End of reading process of an event
      if(event_weight < 0) negative += 1;  
      
    ///Determine the final state
	 if(n_ele == 4 && n_mu == 0) FinalState = 0;
    else if(n_ele == 0 && n_mu == 4) FinalState = 1;
    else if(n_ele == 2 && n_mu == 2) FinalState = 2;
    else{
      cout<<"[Error] Final State  nEle: "<<n_ele<<", nMu: "<<n_mu<<" not expected!"<<endl;
      throw exception();
    }
          
    ///Store the event weight
    EventWeight = event_weight;
    Ntrials	= trials;
    
    ///Builds the matrix for FastME analysis
    for(int i=0; i<4; i++){
      ParticleID[i] = partID[i];
      
      RecoParticle[i][0][0] = part4p[i].Pt();	RecoParticle[i][0][1] = 1.;
      RecoParticle[i][1][0] = part4p[i].Eta();	RecoParticle[i][1][1] = 1.;
      RecoParticle[i][2][0] = part4p[i].Phi();	RecoParticle[i][2][1] = 1.;
    }
    ///-------------------------------------------------------------
    
    ///Builds Z on and off-shell based on CMS method
    TLorentzVector pair1, pair2;
    int fc=-1, sc=-1;
    for(int f=0; f<4; f++)
    for(int s=0; s<4; s++){
      if((abs(ParticleID[f]) != abs(ParticleID[s])) || (ParticleID[f] == ParticleID[s])) continue;
      pair1 = part4p[f] + part4p[s];
      fc = f;
      sc = s;
    }
    for(int f=0; f<4; f++)
    for(int s=0; s<4; s++){
      if((abs(ParticleID[f]) != abs(ParticleID[s])) || (ParticleID[f] == ParticleID[s]) || f==fc || f==sc || s==fc || s==sc) continue;
      pair2 = part4p[f] + part4p[s];
    }
    Double_t mpair1 = pair1.M();
    Double_t mpair2 = pair2.M();
    Zon_mass = (fabs(mpair1-Z_mass)<fabs(mpair2-Z_mass))? mpair1:mpair2;
    Zoff_mass = (Zon_mass == mpair1)? mpair2:mpair1; 
    ZZ_mass = (pair1 + pair2).M();
    
    ///Save tree
    lheTree->Fill();
    
    if(Zon_mass>40 && Zon_mass<120) Zon->Fill(Zon_mass);
    if(Zoff_mass>12 && Zoff_mass<120) Zoff->Fill(Zoff_mass);
    if(ZZ_mass>50  &&  ZZ_mass<2000) ZZ->Fill(ZZ_mass);

    ///Clear vectors
    partID.clear();
    part4p.clear();

    Input >> info;
    //cout<<"nenvents: "<<nevents<<endl;
  }while(info != "</LesHouchesEvents>");
   //while(nevents < 10000);
  cout<<"Events Generated: "<<nevents<<endl;
  cout<<"Fraction of Negative Weights: "<<(negative*100)/float(nevents)<<"%"<<endl;
  cout<<"Total XS:         "<<Total_XS/double(ntrials)<<" pb"<<endl;

  
  TFile *lheRoot = new TFile(out_name+".root","recreate");
  lheTree->Write();
  Zon->Write();
  Zoff->Write();
  ZZ->Write();
  lheRoot->Close();

}
