////////////////////////////////////////////////////////////////////////
//  Deuteron Elastic Cross Section Model
//   --- Obtained from D. Higinbotham's fortran code in MCEEP
//   --- Zhihong Ye, Nov. 15 2014
////////////////////////////////////////////////////////////////////////

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

#include "../Constants.h"
#include "../Deut_Elastic.h"

int main(){
	double E0 = 0;//MeV
	double Theta = 0;//degree
	double SinThSQ, ETA,Ep,Q2,Q_FM,XS_M;
	double XS_EX,XS_EX_Err,A,A_Err,B,qfm, q2fm, xbj, A_M, B_M;

	TFile *file = new TFile("deut_simon.root","recreate");
	TTree *T = new TTree("T","A new Tree");
	T->Branch("E0", &E0, "E0/D");
	T->Branch("Ep", &Ep, "Ep/D");
	T->Branch("Theta", &Theta, "Theta/D");
	T->Branch("Q2", &Q2, "Q2/D");
	T->Branch("Q_FM", &Q_FM, "Q_FM/D");
	T->Branch("XS_M", &XS_M, "XS_M/D");
	T->Branch("XS_EX", &XS_EX, "XS_EX/D");
	T->Branch("XS_EX_Err", &XS_EX_Err, "XS_EX_Err/D");
	T->Branch("qfm", &qfm, "qfm/D");
	T->Branch("A", &A, "A/D");
	T->Branch("A_M", &A_M, "A_M/D");
	T->Branch("A_Err", &A_Err, "A_Err/D");
	T->Branch("B", &B, "B/D");
	T->Branch("B_M", &B_M, "B_M/D");
	T->Branch("xbj", &xbj, "xbj/D");

	double Q2_log[100], vXS_EX[100],vXS_MC[100],vXS_EX_Err[100];
	double Q2_Min = 100000, Q2_Max = -1000000, XS_Min = 1e23, XS_Max = -1e23;
    int N =0;

    //E0 Theta Q2 XS(cm2/sr*1e-30) stat_err (A+B*Tan) A
	ifstream inf("Simon.dat");
	TString Com;
	inf >> Com >> Com >> Com >> Com >> Com >> Com >> Com;

	while(inf>> E0 >> Theta >> q2fm >>  XS_EX >> XS_EX_Err >> B >> A){ 

		XS_EX *=  1e3;    //nb/sr, it was cm/sr*1e-30, CM2NB=1e33
		XS_EX_Err *= 1e3; //nb/sr, it was cm/sr*1e-30, CM2NB=1e33

		Theta *= PI/180.0;//rad
		SinThSQ= pow(sin(Theta/2),2);
		ETA= 1./(1.+ 2.*E0/Deut_Mass*SinThSQ);
		Ep= ETA*E0;
		Q2= 4.*E0*Ep*SinThSQ;
		xbj = Q2/2.0/(E0-Ep)/Deut_Mass;
		Q_FM  = sqrt(Q2)/HBARC;

		//                     MeV  Degree
		XS_M = Deut_Elastic(E0, Theta*180/PI)*FM2NB;//fm^2/sr

		Q2_log[N] = log10(Q2);
		vXS_MC[N] = log10(XS_M);
		vXS_EX[N] = log10(XS_EX);
		vXS_EX_Err[N] = log10(XS_EX_Err);
        if(Q2_Min > Q2_log[N]) Q2_Min = Q2_log[N];
        if(Q2_Max < Q2_log[N]) Q2_Max = Q2_log[N];
        if(XS_Min > vXS_EX[N]) XS_Min = vXS_EX[N];
        if(XS_Max < vXS_EX[N]) XS_Max = vXS_EX[N];
		if(XS_Min > vXS_MC[N]) XS_Min = vXS_MC[N];
		if(XS_Max < vXS_MC[N]) XS_Max = vXS_MC[N];
        N++;

		cerr<<Form("----E=%5.2f, Theta=%5.2f, Q_FM=%f(%f),  XS_EX=%e, XS_M = %e, Diff=%f", E0,Theta/PI*180,Q_FM,sqrt(q2fm), XS_EX, XS_M, (XS_EX-XS_M)/XS_EX*100. )<<endl; 

		T->Fill();
	}	

	file->cd(); T->Write(); file->Close();

	TCanvas *c1 = new TCanvas("c1","c1",800,600);
    TH2F *h1 = new TH2F("h1","Elastic Deuteron XS (Simon et.al) vs SOG model ", 100, Q2_Min-1,Q2_Max+1,100,XS_Min/10,XS_Max*10);
    h1->SetXTitle("log_{10}(Q^{2}) (MeV^{2})");
	h1->SetYTitle("log_{10}(#sigma_{el}) (fm^{2}/sr)");
    h1->Draw();

	TGraphErrors* gr = new TGraphErrors(N, Q2_log, vXS_EX, 0, vXS_EX_Err);
	gr->SetMarkerColor(2);
	gr->SetMarkerStyle(20);
	gr->Draw("p");

	TGraphErrors* gr1 = new TGraphErrors(N, Q2_log, vXS_MC, 0, 0);
	gr1->SetMarkerColor(4);
	gr1->SetMarkerStyle(21);
	gr1->Draw("p");

	c1->Print("simon.png");

	return 0;
}
