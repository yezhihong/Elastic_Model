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
#include "../Constants.h"
#include "../He4_Elastic.h"

void Test_He4(){
  gStyle->SetOptStat(0);

  double E0 = 0.0;     //cerr<<" --- E0 (GeV) = "; cin>> E0;
  double Theta = 0.0;  //cerr<<" --- Theta (Degree) = "; cin>> Theta;
  double q2, XS, XS_Err, G_E, G_E_Err;
  double Ep, Q2,q_fm, XS_M;
  int N;

  ifstream infile("Frosch_He4_XS.dat");
  ofstream outfile("Calc_Frosch_He4_XS.dat");
  
  TFile *file = new TFile("He4_Elastic_XS.root","recreate");
  TTree* T = new TTree("T","a new tree");

  T->Branch("E0", &E0, "E0/D");
  T->Branch("Theta", &Theta , "Theta/D");
  T->Branch("Ep", &Ep, "Ep/D");
  T->Branch("q2", &q2, "q2/D");
  T->Branch("Q2", &Q2, "Q2/D");
  T->Branch("q_fm", &q_fm, "q_fm/D");
  T->Branch("XS", &XS, "XS/D");
  T->Branch("XS_Err", &XS_Err, "XS_Err/D");
  T->Branch("XS_M", &XS_M, "XS_M/D");
  T->Branch("G_E", &G_E, "G_E/D");
  T->Branch("G_E_Err", &G_E_Err, "G_E_Err/D");


	double Q2_log[100], vXS_EX[100],vXS_MC[100],vXS_EX_Err[100];
	double vXS_Diff[100],vXS_Diff_Err[100];
	double Q2_Min = 100000, Q2_Max = -1000000, XS_Min = 1e23, XS_Max = -1e23;
    int C =0;


  TString Com;
//E(MeV) Angle q2(fm-2) xs(cm2/sr) xs_err power(10^-N) G_E    G_E err
  while(!infile.eof()){
      infile >> E0 >> Theta >> q2 >> XS >> XS_Err >> N >> G_E >> G_E_Err;

	  cerr<<Form(" -- E0=%f, Theta=%f", E0, Theta)<<endl;

	  XS *= pow(10,-N) * CM2NB;
	  XS_Err *= pow(10,-N) * CM2NB;

	  Theta *=PI/180.0;
	  double SinThSQ= pow(sin(Theta/2),2);
	  double ETA= 1./(1.+ 2.*E0/He4_Mass*SinThSQ);
	  Ep= ETA*E0;
	  Q2= 4.*E0*Ep*SinThSQ;
	  q_fm=sqrt(Q2)/HBARC;           //q4 in fm**-1

	  XS_M = He4_Elastic(E0,Theta*180/PI) * FM2NB;
	  cerr<<Form("--- Q2 = %f (GeV^2), q = %f (fm-1), XS_M = %e nb/sr, XS_EX=%e fm^2/sr ", Q2/GeV2MeV/GeV2MeV, q_fm, XS_M, XS)<<endl;

	  Q2_log[C] = log10(Q2);
	  vXS_MC[C] = log10(XS_M);
	  vXS_EX[C] = log10(XS);
	  vXS_EX_Err[C] = log10(XS_Err);
	  vXS_Diff[C] = XS_M/XS;
	  vXS_Diff_Err[C] = vXS_Diff[C] * XS_Err/XS;

	  if(Q2_Min > Q2_log[C]) Q2_Min = Q2_log[C];
	  if(Q2_Max < Q2_log[C]) Q2_Max = Q2_log[C];
	  if(XS_Min > vXS_EX[C]) XS_Min = vXS_EX[C];
	  if(XS_Max < vXS_EX[C]) XS_Max = vXS_EX[C];
	  if(XS_Min > vXS_MC[C]) XS_Min = vXS_MC[C];
	  if(XS_Max < vXS_MC[C]) XS_Max = vXS_MC[C];
	  C++;

	 
	  outfile << Form("%8.3f %8.3f %8.3f %14.4e %14.4e %14.4e %14.4e %14.4e", E0, Theta, q2, XS_M, XS, XS_Err) <<endl;

	  T->Fill();
  }
  infile.close();
  outfile.close();
  file->cd(); T->Write(); file->Close();

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->Divide(1,2);
  c1->cd(1);
  TH2F *h1 = new TH2F("h1","Elastic He4 XS (Frosh et.al) vs SOG model ", 100, Q2_Min-1,Q2_Max+1,100,XS_Min-1,XS_Max+2);
  h1->SetXTitle("log_{10}(Q^{2}) (MeV^{2})");
  h1->SetYTitle("log_{10}(#sigma_{el}) (fm^{2}/sr)");
  h1->Draw();

  TGraphErrors* gr = new TGraphErrors(C, Q2_log, vXS_EX, 0, vXS_EX_Err);
  gr->SetMarkerColor(2);
  gr->SetLineColor(2);
  gr->SetMarkerStyle(20);
  gr->Draw("p");

  TGraphErrors* gr1 = new TGraphErrors(C, Q2_log, vXS_MC, 0, 0);
  gr1->SetMarkerColor(4);
  gr1->SetMarkerStyle(21);
  gr1->Draw("p");

  c1->cd(2);
  TH2F *h2 = new TH2F("h2","Elastic He4 XS (Frosch et.al) vs SOG model ", 100, Q2_Min-1,Q2_Max+1,100,0,2);
  h2->SetXTitle("log_{10}(Q^{2}) (MeV^{2})");
  h2->SetYTitle("(#sigma_{el}^{MC}/#sigma_{el}^{EX}");
  h2->Draw();

  TGraphErrors* gr2= new TGraphErrors(C, Q2_log, vXS_Diff, 0, vXS_Diff_Err);
  gr2->SetMarkerColor(1);
  gr2->SetLineColor(1);
  gr2->SetMarkerStyle(20);
  gr2->Draw("p");



  c1->Print("simon.png");



}
