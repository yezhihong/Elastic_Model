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

	//E0(MeV) Q2(fm-2) angle(degree) XS (mb/sr) XS_Err 10-N B(Q2) 10-N A(Q2) A_Err 10-N
	ifstream inf("Platchkov_deut.txt");
	TString Com;
	int N1, N2,N3;

	/*Platchkov_deut{{{*/
		while(inf>> E0 >> qfm >> Theta >> XS_EX >> XS_EX_Err >> N1 >> B >> N2 >> A >> A_Err >> N3){ 

	//	cerr<<Form("--- E0=%f, Theta=%f, XS_EX=%f, N1=%d, N2=%d, N3=%d", E0,Theta, XS_EX, N1,N2,N3)<<endl;

	XS_EX *= pow(10,-N1)/FM2MB; //fm^2/sr
	XS_EX_Err *= pow(10,-N1)/FM2MB; //fm^2/sr
	B *= pow(10,-N2);
	A *= pow(10,-N3);
	A_Err *= pow(10,-N3);
	Theta *= PI/180.0;//rad

	SinThSQ= pow(sin(Theta/2),2);
	ETA= 1./(1.+ 2.*E0/Deut_Mass*SinThSQ);
	Ep= ETA*E0;
	Q2= 4.*E0*Ep*SinThSQ;
	xbj = Q2/2.0/(E0-Ep)/Deut_Mass;
	Q_FM  = sqrt(Q2)/HBARC;

	//                     MeV  Degree
	XS_M = Deut_Elastic(E0, Theta*180/PI);//fm^2/sr

	cerr<<Form("----E=%5.2f, Theta=%5.2f, Q_FM=%f,  XS_EX=%e (%e)  XS_M = %e, Diff = %f", E0,Theta/PI*180,Q_FM, XS_EX,XS_EX_Err, XS_M, (XS_EX-XS_M)/XS_EX*100. )<<endl; 

	T->Fill();
	}
	/*}}}*/

	file->cd(); T->Write(); file->Close();

	return 0;
}
