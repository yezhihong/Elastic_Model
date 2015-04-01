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
#include "Elastic_Main.h"
#include "He4_Elastic.h"
#include "QuasiElastic.h"

int main(){

  double E0 = 0.0;     //cerr<<" --- E0 (GeV) = "; cin>> E0;
  double Theta = 0.0;  //cerr<<" --- Theta (Degree) = "; cin>> Theta;
  double q2, XS, XS_Err, G_E, G_E_Err;
  double Ep, Q2,q_fm, Sig;
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
  T->Branch("Sig", &Sig, "Sig/D");
  T->Branch("Sig", &Sig, "Sig/D");
  T->Branch("G_E", &G_E, "G_E/D");
  T->Branch("G_E_Err", &G_E_Err, "G_E_Err/D");

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

	  Sig = he4_elastic(E0,Theta) * FM2NB;
	  cerr<<Form("--- Q2 = %f (GeV^2), q = %f (fm-1), Sig = %e nb/sr, XS_EX=%e fm^2/sr ", Q2/GeV2MeV/GeV2MeV, q_fm, Sig, XS)<<endl;
 
	  outfile << Form("%8.3f %8.3f %8.3f %14.4e %14.4e %14.4e %14.4e %14.4e", E0, Theta, q2, Sig, XS, XS_Err) <<endl;

	  T->Fill();
  }
  infile.close();
  outfile.close();
  file->cd(); T->Write(); file->Close();
  return 0;
}
