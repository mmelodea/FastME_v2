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

void format_MadGraphLHEtoRoot4l(){
  
  ifstream Input("");
  TString out_name = "";

  //ifstream Input("/home/sabayon/MG5_aMC_v2_2_3/work_area/qq2e2mu/Events/run_02/unweighted_events.lhe");
  //TString out_name = "qq2e2mu_MadGraph";
  
  string status;
  Double_t Total_XS;

  TH1D *Zon = new TH1D("Zon","Z On-Shell",40,40,120);
  TH1D *Zoff = new TH1D("Zoff","Z Off-Shell",54,12,120);
  TH1D *ZZ = new TH1D("ZZ","ZZ Invariant Mass",195,50,2000);  
  
  Int_t ParticleID[4], FinalState, EventType = 0;
  Double_t Zon_mass, Zoff_mass, ZZ_mass, EventWeight, RecoParticle[4][3][2];
  TTree *lheTree = new TTree("LHE_Tree","LHE Tree formated to FastME");
  //lheTree->Branch("EventType",&EventType,"EventType/I");
  lheTree->Branch("ParticleID",&ParticleID,"ParticleID[4]/I");
  lheTree->Branch("FinalState",&FinalState,"FinalState/I");
  lheTree->Branch("RecoParticle",&RecoParticle,"RecoParticle[4][3][2]/D");
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
  
  int nevents = 0, n_ele, n_mu, negative;
  Double_t event_weight;
  do{
      ///Counting event number
      nevents += 1;
      
      ///Reseting tree variables
      fv_tmp	  = reset4p;
      Zon_mass    = pedestal;
      Zoff_mass   = pedestal;
      ZZ_mass     = pedestal;
      
      ///Checks event by event
      n_ele = 0; n_mu = 0;
      do{
	  status = "";
	  Input >> info;
	  if(info == "") continue; ///Dumps empity spaces
	  
	  ///Searches for electrons/ positrons
	  if(info == "11" || info == "-11"){
	    n_ele += 1;
	         if(info ==  "11") partID.push_back(11);
	    else if(info == "-11") partID.push_back(-11);
	    
	    Input >> status;
	    if(status != "1") continue; ///Checks if particle is outgoing final state
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
	    n_mu += 1;
	         if(info ==  "13") partID.push_back(13);
	    else if(info == "-13") partID.push_back(-13);

	    Input >> status;
	    if(status != "1") continue; ///Checks if particle is outgoing final state
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
      
    ///Determine the final state
	 if(n_ele == 4 && n_mu == 0) FinalState = 0;
    else if(n_ele == 0 && n_mu == 4) FinalState = 1;
    else if(n_ele == 2 && n_mu == 2) FinalState = 2;
    else{
      cout<<"[Error] Final State  nEle: "<<n_ele<<", nMu: "<<n_mu<<" not expected!"<<endl;
      throw exception();
    }
          
    
    ///Builds the matrix for FastME analysis
    for(int i=0; i<4; i++){
      ParticleID[i] = partID[i];
      
      RecoParticle[i][0][0] = part4p[i].Pt();	RecoParticle[i][0][1] = 1.;
      RecoParticle[i][1][0] = part4p[i].Eta();	RecoParticle[i][1][1] = 1.;
      RecoParticle[i][2][0] = part4p[i].Phi();	RecoParticle[i][2][1] = 1.;
    }
    ///-------------------------------------------------------------
    //cout<<"\nEvent: "<<nevents<<endl;
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
    //cout<<"fc: "<<fc<<"  sc: "<<sc<<endl;
    for(int f=0; f<4; f++)
    for(int s=0; s<4; s++){
      if((abs(ParticleID[f]) != abs(ParticleID[s])) || (ParticleID[f] == ParticleID[s]) || f==fc || f==sc || s==fc || s==sc) continue;
      //cout<<"f: "<<f<<"  s: "<<s<<endl;
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
  cout<<"Events Generated: "<<nevents<<endl;
  //cout<<"Total XS:         "<<Total_XS<<" pb"<<endl;
  
  TFile *lheRoot = new TFile(out_name+".root","recreate");
  lheTree->Write();
  Zon->Write();
  Zoff->Write();
  ZZ->Write();
  lheRoot->Close();

}
