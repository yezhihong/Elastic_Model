
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
#include "Moller.h"
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
	cerr<<"--- Beam Energy (MeV), E0 = "; cin >> E0;

	double Ep, Q2,qfm,xbj;
	double Theta, Phi, XS;

	TFile *file = new TFile(Form("XSoller_%d_Smear.root",(int) (E0) ),"recreate");
	TTree* T = new TTree("T","a new tree");

	T->Branch("E0", &E0, "E0/D");
	T->Branch("Theta", &Theta , "Theta/D");
	T->Branch("Phi", &Phi, "Phi/D");
	T->Branch("Ep", &Ep, "Ep/D");
	T->Branch("Q2", &Q2, "Q2/D");
	T->Branch("qfm", &qfm, "qfm/D");
	T->Branch("XS", &XS, "XS/D");
	T->Branch("xbj", &xbj, "xbj/D");

	double Q2_gen, Ep_gen, Theta_gen, Phi_gen;
	T->Branch("Phi_gen", &Phi_gen, "Phi_gen/D");
	T->Branch("Theta_gen", &Theta_gen , "Theta_gen/D");
	T->Branch("Ep_gen", &Ep_gen, "Ep_gen/D");
	T->Branch("Q2_gen", &Q2_gen, "Q2_gen/D");

	const int N = 20000000;
	gRandom->SetSeed(0);

	for(int i=0;i<N;i++){
		cout<<"--- N="<< i <<"\r";
		Theta_gen = gRandom->Uniform(0.4,8.0);	//degree
		Phi_gen = gRandom->Uniform(0.,360.0);	//degree

		double Theta_Rad = Theta_gen *PI/180.0;//Rad
		double SinThSQ= pow(sin(Theta_Rad/2),2);
		Q2_gen= 4.*E0*Ep_gen*SinThSQ;

		double P0 = sqrt( E0*E0 - Electron_Mass*Electron_Mass);
		double Theta_CM = 2.0 * atan( tan(Theta_Rad) * sqrt( (P0+Electron_Mass)/(2.0*Electron_Mass) ) );
		double CosTh_CM = cos(Theta_CM);
        Ep_gen = E0/2. * (1.+CosTh_CM);

		XS =  Moller(E0, Theta_gen);
        XS *= FM2NB;

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

	file->cd(); T->Write(); file->Close();
	return 0;
}
