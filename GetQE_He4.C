
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
#include "XEMC.h"
//To make this code work, please go to the file 
//"./QE_XEMC/SRC/XEM_Target.h"
//and, make sure the "TARGET_TABLE" is pointed to the file:
//"target.table", which is in "./QE_XEM/SRC/target.table"
//

const double Delta_E = 27.0;//MeV
const double Delta_X = 0.02; //0.02 cm from GEM
const double Delta_Y = 0.02; //0.02 cm from GEM
const double Length_HyCal = 500; //cm from target to Hycal

int main(){

	double E0 = 0.;//MeV
	cerr<<"--- Beam Energy (MeV) E0 = "; cin >> E0;

	double Ep, Q2,qfm,xbj;
	double Theta, Phi, XS, XS_QE, XS_DIS;

	TFile *file = new TFile(Form("XS_He4_QE_%d_F1F2QE09.root",(int) (E0)),"recreate");
	TTree* T = new TTree("T","a new tree");

	T->Branch("E0", &E0, "E0/D");
	T->Branch("Theta", &Theta , "Theta/D");
	T->Branch("Phi", &Phi, "Phi/D");
	T->Branch("Ep", &Ep, "Ep/D");
	T->Branch("Q2", &Q2, "Q2/D");
	T->Branch("qfm", &qfm, "qfm/D");
	T->Branch("XS", &XS, "XS/D");
	T->Branch("XS_QE", &XS_QE, "XS_QE/D");
	T->Branch("XS_DIS", &XS_DIS, "XS_DIS/D");
	T->Branch("xbj", &xbj, "xbj/D");

	double Q2_gen,Ep_gen, Theta_gen, Phi_gen,xbj_gen, qfm_gen;
	T->Branch("Phi_gen", &Phi_gen, "Phi_gen/D");
	T->Branch("Theta_gen", &Theta_gen , "Theta_gen/D");
	T->Branch("Ep_gen", &Ep_gen, "Ep_gen/D");
	T->Branch("Q2_gen", &Q2_gen, "Q2_gen/D");
	T->Branch("xbj_gen", &xbj_gen, "xbj_gen/D");
	T->Branch("qfm_gen", &qfm_gen, "qfm_gen/D");
	
	/*Setup for QE XS{{{*/
	const int A=4; 
	const int Z=2; 
	TString Target_Input = "./QE_XEMC/He4_Input.dat";	

	//Define a event to calculate radiated cross section
	XEMC* Event = new XEMC(); 
	Event->Init(Target_Input.Data());
	/*}}}*/

	const int N = 20000000;
	gRandom->SetSeed(0);
	for(int i=0;i<N;i++){
		cout<<"--- N="<< i <<"\r";
		Theta_gen = gRandom->Uniform(0.4,8.0);	//degree
		Phi_gen = gRandom->Uniform(0.,360.0);	//degree
		Ep_gen = gRandom->Uniform(1., E0*1.);//MeV

		double Theta_Rad = Theta_gen *PI/180.0;//Rad
		double SinThSQ= pow(sin(Theta_Rad/2),2);
		Q2_gen= 4.*E0*Ep_gen*SinThSQ;	
		qfm_gen=sqrt(Q2_gen)/HBARC;           //q4 in fm**-1
		xbj_gen = Q2_gen/(2.*Proton_Mass*(E0-Ep_gen));

		//Quasi-Elastic Scattering
		int err = -1000;
		if(xbj_gen<=4){
			//                          GeV  GeV  Degree
			err =Event->Process(E0/1000.,Ep_gen/1000.,Theta_gen,A,Z,0.0);	
			if(err>=0){
				//nb/sr
				//Go to H2_Input.dat, Line#23 set to one to enable RadC
				//Otherwise, cs_rad = cs_born
				XS_QE = Event->XS_QE(); //Non-Radiated, in nb/sr/MeV
				XS_DIS = Event->XS_DIS(); //Non-Radiated, in nb/sr/MeV
				//XS = Event->XS_Rad(); //Radiated, in nb/sr/MeV
				XS = Event->XS_Born(); //Non-Radiated, in nb/sr/MeV
			}else{
				XS = 1e-20;
			}
			if(isnan(XS))
				XS = 1e-20;//Fail to calculate XS
		}else{
			XS = 1e-20;
		}

		/////////////////////////////////////////////////
		//Now Smearing with the resolution:
		/////////////////////////////////////////////////
         Ep = gRandom->Gaus(Ep_gen, Delta_E);

		double R_L = Length_HyCal / cos(Theta_gen *PI/180.0);
		double x_gen = R_L * sin(Theta_gen *PI/180.0) * cos(Phi_gen *PI/180.0); 
		double y_gen = R_L * sin(Theta_gen *PI/180.0) * sin(Phi_gen *PI/180.0); 
        double x = gRandom->Gaus(x_gen, Delta_X);
        double y = gRandom->Gaus(y_gen, Delta_Y);
		Theta = asin( sqrt( (x*x+y*y)/R_L/R_L ) )*180/PI;
		Phi = atan(y/x)*180/PI;
		if(x<0. && y<0.) Phi+=180.0;
		if(x>0. && y<0.) Phi+=270.0;
        ////////////////////////////////////////////////

		Theta_Rad = Theta *PI/180.0;//Rad
		SinThSQ= pow(sin(Theta_Rad/2),2);
		Q2= 4.*E0*Ep*SinThSQ;
		qfm=sqrt(Q2)/HBARC;           //q4 in fm**-1
		xbj = Q2/(2.*Proton_Mass*(E0-Ep));

		T->Fill();

	}

	file->cd(); T->Write();file->Close();
	delete Event;
	return 0;
}
