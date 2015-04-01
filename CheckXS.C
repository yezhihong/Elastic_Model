//---------------------------------------------------------------------
//
//      SUBROUTINE HE4_ELASTIC
//
//      AUTHOR: D.W.Higinbotham
//      DATE:   February 2000
//      PURPOSE:
//              Calculate 4He cross sections for elastic scattering.
//              The form factor  calculations use the fit of  C. R. 
//              Ottermann et  al.  Nucl. Phys. A436 (1985) 688-698.
//              The Coulomb correction is made by using Q effective.  
//
//      Update: Convert into C++ by Zhihong Ye,   Nov/15/2014
//---------------------------------------------------------------------

/*Include{{{*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <errno.h>
#include <sstream>
/*ROOT Includes{{{*/
#include <TSystem.h>
#include <TString.h>
#include <TStyle.h>
#include <Riostream.h>
#include "TObjString.h"
#include <TNamed.h>
#include <TPRegexp.h>
#include <TObjArray.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TError.h>
#include <TVirtualFitter.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TCut.h>
#include <TMultiGraph.h>
#include <TCutG.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
//#include <TMatrix.h>
/*}}}*/

/*}}}*/
using namespace::std;
using namespace::TMath;
#include "Constants.h"
#include "He4_Elastic.h"
#include "Deut_Elastic.h"
#include "Moller.h"
#include "XEMC.h"
//To make this code work, please go to the file 
//"./QE_XEMC/SRC/XEM_Target.h"
//and, make sure the "TARGET_TABLE" is pointed to the file:
//"target.table", which is in "./QE_XEM/SRC/target.table"
//
int main(){
	//void CheckXS(){

	double E0 = 1100.;//MeV
	double Ep,Ep_E,Q2,Q2_E,qfm,qfm_E;
	double Theta, Phi, XS_QE, XS_E, XS_M;

	TString Target = ""; cerr<<"--- Target (Deut or He4) "; cin >> Target;
	cerr<<"=== E0 = (MeV) "; cin >> E0;
	cerr<<"=== Ep = (MeV) "; cin >> Ep;
	cerr<<"=== Theta = (Deg) "; cin >> Theta;

	/*Setup for QE XS{{{*/
	int A=4; 
	int Z=2; 
	TString Target_Input = "./QE_XEMC/He4_Input.dat";	

	if(Target=="Deut"){
		A = 2; Z = 1; 
		Target_Input = "./QE_XEMC/H2_Input.dat";	
	}

	//Define a event to calculate radiated cross section
	XEMC* Event = new XEMC(); 
	Event->Init(Target_Input.Data());
	/*}}}*/

	double Theta_Rad = Theta *PI/180.0;//Rad
	double SinThSQ= pow(sin(Theta_Rad/2),2);
	double ETA= 1./(1.+ 2.*E0/He4_Mass*SinThSQ);
	Q2= 4.*E0*Ep*SinThSQ;
	qfm=sqrt(Q2)/HBARC;           //q4 in fm**-1

	Ep_E= ETA*E0;
	Q2_E= 4.*E0*Ep_E*SinThSQ;
	qfm_E=sqrt(Q2_E)/HBARC;           //q4 in fm**-1

	//Elastic Scattering
	if(Target=="Deut")
		XS_E = Deut_Elastic(E0,Theta);//fm2/sr
	else if(Target=="He4")
		XS_E = He4_Elastic(E0,Theta);//fm2/sr

	XS_E *= FM2NB;//nb/sr

	//Quasi-Elastic Scattering
	//                          GeV  GeV  Degree
	int err = Event->Process(E0/1000.,Ep/1000.,Theta,A,Z,0.0);	
	if(err>=0){
		//Go to H2_Input.dat, Line#23 set to one to enable RadC
		//Otherwise, cs_rad = cs_born
		//XS_QE = Event->XS_Rad(); //Radiated, in nb/sr/MeV
		XS_QE = Event->XS_Born(); //Non-Radiated, in nb/sr/MeV
	}

	//Moller Scattering
	XS_M = Z * Moller(E0, Theta);//fm2/sr
	XS_M *= FM2NB;//nb/sr

	cerr<<Form("--- E0=%6.2f, Ep=%6.2f(%6.2f), Theta=%6.2f, q2fm=%f(%f), Q2=%12.3e (%e) (GeV-2) ETA=%f",E0,Ep,Ep_E,Theta, qfm,qfm_E, Q2/1000./1000.,Q2_E/1000./1000., ETA)<<endl;
	cerr<<Form("&&& Elastic=%12.3e nb/sr, QE=%12.3e nb/sr/MeV, Moller=%10.4e nb/sr ",XS_E, XS_QE, XS_M)<<endl;
	return 0;
}
