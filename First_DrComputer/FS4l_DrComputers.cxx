///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::		PHASE ESPACE DATA-MC DISTANCE COMPUTERS	LIBRARY		::::::::::
///::::::		       Author: Miquéias M. de Almeida			::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::: 4 LEPTONS FINAL STATE CASE :::::::::::::::::::::::::::::::

#define pi 3.14159265358979312
#include "FS4l_DrComputers.h"
#include <TString.h>
#include <iostream>


///DISTANCE ORDERING MEDIUM + RESSONANCE CRITERY
///TAKES THE MEDIUM IN dPt, dEta, dPhi FOR EACH PAIR
double FS4l_DrComputers(Double_t Data[4][3][2], Double_t MC[4][3][2], TString Model){
  int SetModel = -1;
       if(Model == "DrOrder") SetModel = 0;
  else if(Model == "DrMedio") SetModel = 1;
  else cout<<"Model '"<<Model<<"' not defined!"<<endl;
  
  switch(SetModel){
    case 0:
      return FS4l_DrOrder(Data,MC);
      break;

    case 1:
      return FS4l_DrMedio(Data,MC);
      break;

    default:
      return FS4l_DrOrder(Data,MC);
      break;
  }
}


///DR ORDERING
double FS4l_DrOrder(Double_t Data[4][3][2], Double_t MC[4][3][2]){
  double fdpt2, fdeta2, fdphi2, sdpt2, sdeta2, sdphi2, sum_dr1, sum_dr2;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int start = 0;
  
  do{
      fdpt2   = 0;
      fdeta2  = 0;
      fdphi2  = 0;
      sdpt2   = 0;
      sdeta2  = 0;
      sdphi2  = 0;
      sum_dr1 = 0;
      sum_dr2 = 0;
	
      fdpt2  = pow( (Data[start][0][0]-MC[start][0][0])/Data[start][0][1] ,2 );
      fdeta2 = pow( (Data[start][1][0]-MC[start][1][0])/Data[start][1][1] ,2 );
      
      if(Data[start][2][0]-MC[start][2][0] > pi)
	fdphi2 = pow( (Data[start][2][0]-MC[start][2][0]-pi)/Data[start][2][1] ,2 );
      else
	fdphi2 = pow( (Data[start][2][0]-MC[start][2][0])/Data[start][2][1] ,2 );
      
      sdpt2  = pow( (Data[start][0][0]-MC[start+1][0][0])/Data[start][0][1] ,2 );
      sdeta2 = pow( (Data[start][1][0]-MC[start+1][1][0])/Data[start][1][1] ,2 );
      
      if(Data[start][2][0]-MC[start+1][2][0] > pi)
	sdphi2 = pow( (Data[start][2][0]-MC[start+1][2][0]-pi)/Data[start][2][1] ,2 );
      else
	sdphi2 = pow( (Data[start][2][0]-MC[start+1][2][0])/Data[start][2][1] ,2 );
      
      sum_dr1 = sqrt(fdpt2 + fdeta2 + fdphi2);
      sum_dr2 = sqrt(sdpt2 + sdeta2 + sdphi2);
      
      if(sum_dr1 < sum_dr2){
	sum_dpt2  += fdpt2;
	sum_deta2 += fdeta2;
	sum_dphi2 += fdphi2;
	
	sum_dpt2  += pow( (Data[start+1][0][0]-MC[start+1][0][0])/Data[start+1][0][1] ,2 );
	sum_deta2 += pow( (Data[start+1][1][0]-MC[start+1][1][0])/Data[start+1][1][1] ,2 );
	if(Data[start+1][2][0]-MC[start+1][2][0] > pi)
	  sum_dphi2 += pow( (Data[start+1][2][0]-MC[start+1][2][0]-pi)/Data[start+1][2][1] ,2 );
	else
	  sum_dphi2 += pow( (Data[start+1][2][0]-MC[start+1][2][0])/Data[start+1][2][1] ,2 );
      }
      
      if(sum_dr2 < sum_dr1){
	sum_dpt2  += sdpt2;
	sum_deta2 += sdeta2;
	sum_dphi2 += sdphi2;
	
	sum_dpt2  += pow( (Data[start+1][0][0]-MC[start][0][0])/Data[start+1][0][1] ,2 );
	sum_deta2 += pow( (Data[start+1][1][0]-MC[start][1][0])/Data[start+1][1][1] ,2 );
	if(Data[start+1][2][0]-MC[start][2][0] > pi)
	  sum_dphi2 += pow( (Data[start+1][2][0]-MC[start][2][0]-pi)/Data[start+1][2][1] ,2 );
	else
	  sum_dphi2 += pow( (Data[start+1][2][0]-MC[start][2][0])/Data[start+1][2][1] ,2 );
      }
      
      start += 2;
  }while(start < 4);
  
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}


///MEDIA
double FS4l_DrMedio(Double_t Data[4][3][2], Double_t MC[4][3][2]){
  int idt, imc;
  double sum_dPt2, sum_dEta2, sum_dPhi2, event_distance;
  double sum_dPt2_med = 0, sum_dEta2_med = 0, sum_dPhi2_med = 0;
  
  for(idt=0; idt<4; idt++){
    sum_dPt2  = 0;
    sum_dEta2 = 0;
    sum_dPhi2 = 0;
    
    for(imc=0; imc<4; imc++){
      if( (idt<2 && imc>1) || (idt>1 && imc<2) ) continue; ///To avoid pair mixture
	
      sum_dPt2  += pow( (Data[idt][0][0]-MC[imc][0][0])/Data[idt][0][1] ,2 );
      sum_dEta2 += pow( (Data[idt][1][0]-MC[imc][1][0])/Data[idt][1][1] ,2 );
      
      if( fabs(Data[idt][2][0]-MC[imc][2][0]) > pi)
	sum_dPhi2 += pow( (Data[idt][2][0]-MC[imc][2][0]-pi)/Data[idt][2][1] ,2 );
      else
	sum_dPhi2 += pow( (Data[idt][2][0]-MC[imc][2][0])/Data[idt][2][1] ,2 );
    }
    
    ///Computing the media of pair combinations
    sum_dPt2_med  += sum_dPt2/2.;
    sum_dEta2_med += sum_dEta2/2.;
    sum_dPhi2_med += sum_dPhi2/2.;
  }
      
  event_distance = sqrt(sum_dPt2_med + sum_dEta2_med + sum_dPhi2_med);
  
  return event_distance;
}


///========================================================================================