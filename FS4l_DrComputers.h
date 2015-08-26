///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::		PHASE ESPACE DATA-MC DISTANCE COMPUTERS	LIBRARY		::::::::::
///::::::		       Author: Miquéias M. de Almeida			::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifndef FS4l_DrComputers_h
#define FS4l_DrComputers_h

#include <TString.h>

///FINAL STATE 4 LEPTONS
double	FS4l_DrComputers(Float_t Data[4][3][2], Float_t MC[4][3][2], TString Model);
  double	FS4l_DrOrder(Float_t Data[4][3][2], Float_t MC[4][3][2]);
  double	FS4l_DrMedio(Float_t Data[4][3][2], Float_t MC[4][3][2]);
  
#endif
