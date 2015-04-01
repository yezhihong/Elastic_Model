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
#include <TLegend.h>
//#include <TMatrix.h>
/*}}}*/

/*}}}*/
using namespace::std;
using namespace::TMath;

const double Day=24*60*60;//Sec
const double PI = 3.1415926;
int main(){

	TString Target;
    cerr<<"--- Which Target (Deut,He4) = "; cin >> Target;

	int E0 = 0;
    cerr<<"--- Which Beam Energy (MeV) = "; cin >> E0;

	double Beam_Time = 1.0*Day; //Day->Sec
	if(E0==1100)
		Beam_Time = 1.0*Day; //Day->Sec
	else if(E0==2200)
		Beam_Time = 1.0*Day; //Day->Sec
	else if(E0==3300)
		Beam_Time = 1.0*Day; //Day->Sec
	else if(E0==4400)
		Beam_Time = 1.0*Day; //Day->Sec
    else{
	    cerr<<"*ERROR*, I don't know this beam energy"<<endl;
		return -1;
	}

    /*Pre-Define{{{*/
	const double Deut_Mass = 1875.630;       //Deuteron mass in MeV
	const double He4_Mass=3727.40841;         //4He Nuclear Mass of 3He           
	const double FM2NB = 1e+7; //from ub to fm^2, 1fm^2 = 10mb, 1b=1e-28 m^2, 1fm = 1e-15m
	const double NB2CM = 1e-33; //nb/sr->cm2/sr

    const double Qe = 1.6022 * 1e-19; //C/e	
	const double Beam_Current = 10 * 1e-9; //nA->A
	const double Beam_Electron = (Beam_Current * Beam_Time)/Qe;

//	const double Target_Denstiy = 0.; // g/cm^3
//	const double Target_D = 4.0/10.; //4mm->0.4cm;
//	const double Target_Area = PI*pow(0.5*Target_D,2); //cm2
	const double Deut_Target_Atom = 0.95*1e18; //atoms/cm2
	const double He4_Target_Atom  = 0.45*1e18; //atoms/cm2

	const double Gen_Theta_Min = 0.5 * PI/180.;//Min. value in the MC
	const double Gen_Theta_Max = 8.0 * PI/180.;//Min. value in the MC
	const double dOmega = 2*PI * (Gen_Theta_Max-Gen_Theta_Min);//I will divide it by the N_gen_MC later
    const double dEp = (E0 - 0); //Ep phase space for QE process only
	
	const double Det_Eff = 0.90;//assuming 90% for all efficiency, including HyCal and GEM
	/*}}}*/

    /*Input,Output,Histograms{{{*/
	ofstream outfile(Form("%s_Bin_%d.dat",Target.Data(),E0));
	outfile<<Form("%8s %16s %16s %16s %16s %16s", "E0 (MeV)","Theta (Degree)", "Ep (MeV)", "Q2 (MeV2)","XS (nb/sr)", "XS_Err (nb/sr)")<<endl;

	TChain* TE = new TChain("T");
	TE->Add(Form("../XS_%s_Elastic_%d.root",Target.Data(),E0));
	const int N_E_Gen = TE->GetEntries();
	
	TChain* TQ = new TChain("T");
	TQ->Add(Form("../XS_%s_QE_%d.root",Target.Data(),E0));
	const int N_Q_Gen = TQ->GetEntries();

	TChain* TM = new TChain("T");
	TM->Add(Form("../XS_Moller_%d.root",E0));
	const int N_M_Gen = TM->GetEntries();
    /*}}}*/

	double Target_Atom = 0.0, Target_Mass = 0.0, Z = 0;
	if(Target=="Deut"){ Target_Atom = Deut_Target_Atom; Target_Mass = Deut_Mass; Z = 1;}
	if(Target=="He4") { Target_Atom =  He4_Target_Atom; Target_Mass =  He4_Mass; Z = 2;}

	const double Norm_Factor_E = NB2CM * Det_Eff * Beam_Electron * Target_Atom * dOmega / N_E_Gen;//I divide it by the N_gen_MC here 
	const double Norm_Factor_Q = NB2CM * Det_Eff * Beam_Electron * Target_Atom * dOmega / N_Q_Gen * dEp;
	const double Norm_Factor_M = NB2CM * Det_Eff * Beam_Electron * Target_Atom * dOmega / N_M_Gen * Z;
    //cerr<<Form("---- Eff=%f, Ne = %e, Ntg = %e, dOMega=%f", Det_Eff, Beam_Electron, Target_Atom, dOmega)<<endl;	
	
	cerr<<Form("---     EL Raw Events = %d for E0= %d, Norn_Factor = %e",N_E_Gen,E0, Norm_Factor_E)<<endl;
	cerr<<Form("---     QE Raw Events = %d for E0= %d, Norn_Factor = %e",N_Q_Gen,E0, Norm_Factor_Q)<<endl;
	cerr<<Form("--- Moller Raw Events = %d for E0= %d, Norn_Factor = %e",N_M_Gen,E0, Norm_Factor_M)<<endl;

	/*Get Bins for EL only{{{*/
	const int Bin = 1000;
	TH1D* h_theta = new TH1D("h_theta",Form("%s Rate E0=%d MeV",Target.Data(),E0),Bin, -0.41,0.91);
	h_theta->SetXTitle("log10(#theta) (Degree)");
	h_theta->SetYTitle("Counts");

	TCanvas *c1 = new TCanvas("c1","c1",800,600);
	TE->Draw("log10(Theta)>>h_theta",(TCut)(Form("XS*%e", Norm_Factor_E)),"goff");
	c1->Print(Form("%s_Theta_ELNorm_%d.png", Target.Data(),E0));

	const int Theta_Bin=200;//Make sure it is large enough
	double Theta[Theta_Bin] = {Theta_Bin*0.0};
	double Theta_Min[Theta_Bin] = {Theta_Bin*0.0};
	double Theta_Max[Theta_Bin] = {Theta_Bin*0.0};
	double Theta_Min_Log[Theta_Bin] = {Theta_Bin*0.0};
	double Theta_Max_Log[Theta_Bin] = {Theta_Bin*0.0};
  
	//const double Stat_Bin = 0.01;//Assume the statistical error for each bin is 1%
	const double Stat_Bin_Count = 1e4;//Assume the statistical error for each bin is 1%
  	int Theta_Min_Bin = 20;
  
	double N = 0.; 
	int B = 0, Bin_Temp = 0;
	Theta_Min[0] = pow(10.0, h_theta->GetBinCenter(0));	
	Theta_Min_Log[0] = h_theta->GetBinCenter(0);	
	for(int j=1;j<=Bin;j++){
		if( N<Stat_Bin_Count||Bin_Temp<Theta_Min_Bin){	
			N += (double) (h_theta->GetBinContent(j));
			Bin_Temp++;
			//cerr<<" Histo-Bin = "<< j  <<", N = "<< N <<endl;
			continue;
		}
		else{
			Theta_Max[B] = pow(10., h_theta->GetBinCenter(j));//Degree
			Theta_Max_Log[B] = h_theta->GetBinCenter(j);//Degree

			cerr<<Form("---E=%d MeV, Bin#%d, Theta_Min=%6.4f, Theta_Max=%6.4f, N=%e",E0, B, Theta_Min[B], Theta_Max[B], N)<<endl;
			B++;
			if(B<Theta_Bin){
				Theta_Min[B] = pow(10.,h_theta->GetBinCenter(j));//Degree
				Theta_Min_Log[B] = h_theta->GetBinCenter(j);//Degree
			}

			N = 0.; 
			Bin_Temp = 0;
			if(j<Bin){
				//j--;
				continue;
			}else
				break;
		}
	}
	delete h_theta;
	/*}}}*/
	
	/*Binning{{{*/
	const int Theta_Bin_Final = B;
	double Ep[Theta_Bin] = {Theta_Bin*0.0};
	double Q2[Theta_Bin] = {Theta_Bin*0.0};

	double XS_E_Mean[Theta_Bin] = {Theta_Bin*0.0};
	double XS_E_RMS[Theta_Bin] = {Theta_Bin*0.0};
	double NE_Bin[Theta_Bin] = {Theta_Bin*0.0};
	double NE_Bin_Err[Theta_Bin] = {Theta_Bin*0.0};

	double XS_Q_Mean[Theta_Bin] = {Theta_Bin*0.0};
	double XS_Q_RMS[Theta_Bin] = {Theta_Bin*0.0};
	double NQ_Bin[Theta_Bin] = {Theta_Bin*0.0};
	double NQ_Bin_Err[Theta_Bin] = {Theta_Bin*0.0};

	double XS_M_Mean[Theta_Bin] = {Theta_Bin*0.0};
	double XS_M_RMS[Theta_Bin] = {Theta_Bin*0.0};
	double NM_Bin[Theta_Bin] = {Theta_Bin*0.0};
	double NM_Bin_Err[Theta_Bin] = {Theta_Bin*0.0};

	double Q2_Min = 1e99, Q2_Max = -1e99;
	double XS_Min = 1e99, XS_Max = -1e99;
	double NBin_Min = 1e99, NBin_Max = -1e99;

	TH1D* htemp = new TH1D;
	//TString Phase_Space_Cut = Form("Ep>=%f",0.90*E0);//Allow some cuts on Ep & Theta to clean up some events other than EL electrons
	TString Phase_Space_Cut = Form("1");//Allow some cuts on Ep & Theta to clean up some events other than EL electrons

	for(int i=0;i<Theta_Bin_Final;i++){

		TE->Draw("Theta", Form("XS*%e*(Theta>%e && Theta<=%e && %s)", Norm_Factor_E, Theta_Min[i], Theta_Max[i], Phase_Space_Cut.Data()),"goff");
		htemp = (TH1D*) gROOT->FindObject("htemp");
		Theta[i] = htemp->GetMean();
		NE_Bin[i] = htemp->Integral(); 
		NE_Bin_Err[i] = sqrt(NE_Bin[i]);
		htemp->Reset();

		TQ->Draw("Theta", Form("XS_QE*%e*(Theta>%e && Theta<=%e && %s)", Norm_Factor_Q, Theta_Min[i], Theta_Max[i], Phase_Space_Cut.Data()),"goff");
		htemp = (TH1D*) gROOT->FindObject("htemp");
		NQ_Bin[i] = htemp->Integral();
		NQ_Bin_Err[i] = sqrt(NQ_Bin[i]);
		htemp->Reset();

		TM->Draw("Theta", Form("XS_M*%e*(Theta>%e && Theta<=%e && %s)", Norm_Factor_M, Theta_Min[i], Theta_Max[i], Phase_Space_Cut.Data()),"goff");
		htemp = (TH1D*) gROOT->FindObject("htemp");
		NM_Bin[i] = htemp->Integral();
		NM_Bin_Err[i] = sqrt(NM_Bin[i]);
		htemp->Reset();
			
		double Theta_Rad = Theta[i] *PI/180.0;//Rad
		double SinThSQ= pow(sin(Theta_Rad/2),2);
		double ETA= 1./(1.+ 2.*E0/Target_Mass*SinThSQ);
		Ep[i]= ETA*E0;
		Q2[i]= 4.*E0*Ep[i]*SinThSQ;

		TE->Draw("XS", Form("Theta>%e && Theta<=%e && %s", Theta_Min[i], Theta_Max[i], Phase_Space_Cut.Data()),"goff");
		htemp = (TH1D*) gROOT->FindObject("htemp");
		XS_E_Mean[i] = htemp->GetMean();
		XS_E_RMS[i] = htemp->GetRMS();
		htemp->Reset();

		TQ->Draw("XS_QE", Form("Theta>%e && Theta<=%e && %s", Theta_Min[i], Theta_Max[i], Phase_Space_Cut.Data()),"goff");
		htemp = (TH1D*) gROOT->FindObject("htemp");
		XS_Q_Mean[i] = htemp->GetMean();
		XS_Q_RMS[i] = htemp->GetRMS();
		htemp->Reset();

		TM->Draw("XS_M", Form("Theta>%e && Theta<=%e && %s", Theta_Min[i], Theta_Max[i], Phase_Space_Cut.Data()),"goff");
		htemp = (TH1D*) gROOT->FindObject("htemp");
		XS_M_Mean[i] = htemp->GetMean();
		XS_M_RMS[i] = htemp->GetRMS();
		htemp->Reset();

		outfile<<Form("%5d %16.4f %16.4f %16.4f %16.4e %16.4e %16.4e %16.4e %16.4e", E0, Theta[i],Ep[i],Q2[i], XS_E_Mean[i], XS_E_RMS[i], NE_Bin[i], NQ_Bin[i], NM_Bin[i])<<endl;

		cerr<<Form("E0=%d, Theta_Min=%6.4f, Theta_Max=%6.4f, Theta[%d]=%6.4f, XS=%e, RMS=%e", E0, Theta_Min[i],Theta_Max[i],i, Theta[i], XS_E_Mean[i], XS_E_RMS[i])<<endl;
		cerr<<Form("  N_EL = %e, N_QE = %e, N_M = %e", NE_Bin[i], NQ_Bin[i], NM_Bin[i])<<endl;
		cerr<<Form("---------------------")<<endl;
		
		if(Q2[i]<Q2_Min) Q2_Min = Q2[i];
		if(Q2[i]>Q2_Max) Q2_Max = Q2[i];

		if(XS_E_Mean[i]-XS_E_RMS[i]<XS_Min) XS_Min = XS_E_Mean[i]-XS_E_RMS[i];
		if(XS_E_Mean[i]+XS_E_RMS[i]>XS_Max) XS_Max = XS_E_Mean[i]+XS_E_RMS[i];

		if(NE_Bin[i]+NE_Bin_Err[i]>NBin_Max) NBin_Max = NE_Bin[i]+NE_Bin_Err[i];
		if(NE_Bin[i]-NE_Bin_Err[i]<NBin_Min) NBin_Min = NE_Bin[i]-NE_Bin_Err[i];
		if(NQ_Bin[i]+NQ_Bin_Err[i]>NBin_Max) NBin_Max = NQ_Bin[i]+NQ_Bin_Err[i];
		if(NQ_Bin[i]-NQ_Bin_Err[i]<NBin_Min) NBin_Min = NQ_Bin[i]-NQ_Bin_Err[i];
		if(NM_Bin[i]+NM_Bin_Err[i]>NBin_Max) NBin_Max = NM_Bin[i]+NM_Bin_Err[i];
		if(NM_Bin[i]-NM_Bin_Err[i]<NBin_Min) NBin_Min = NM_Bin[i]-NM_Bin_Err[i];
	}
	/*}}}*/
	delete htemp;

    gStyle->SetOptStat(0);	
	TH2D* h2 = new TH2D("h2",Form("%s Cross Section Binning at E0=%d MeV",Target.Data(), E0), 200, Q2_Min*0.95, Q2_Max*1.05, 200, XS_Min*0.1, XS_Max*10);
	h2->SetXTitle("Q^{2} (MeV)");
	h2->SetYTitle("#sigma (nb/sr)");
	gPad->SetLogy(1);
	h2->Draw();
	TGraphErrors *gr = new TGraphErrors(Theta_Bin_Final, Q2, XS_E_Mean, 0, XS_E_RMS);
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(20);
    gr->Draw("P");
	
	c1->Print(Form("%s_XS_Bin_%d.png",Target.Data(), E0));

	TH2D* h3 = new TH2D("h3",Form("%s Counts for Elastic, QE and Moller Processes at E0=%d MeV",Target.Data(),E0), 200, Q2_Min*0.95, Q2_Max*1.05, 200, NBin_Min*0.1, NBin_Max*10);
	h3->SetXTitle("Q^{2} (MeV)");
	h3->SetYTitle("Count");
	gPad->SetLogy(1);
	h3->Draw();

	TGraphErrors *gr1 = new TGraphErrors(Theta_Bin_Final, Q2, NE_Bin, 0, NE_Bin_Err);
    gr1->SetMarkerColor(2);
    gr1->SetMarkerStyle(20);
    gr1->Draw("P");

	TGraphErrors *gr2 = new TGraphErrors(Theta_Bin_Final, Q2, NQ_Bin, 0, NQ_Bin_Err);
    gr2->SetMarkerColor(4);
    gr2->SetMarkerStyle(21);
    gr2->Draw("P");

	TGraphErrors *gr3 = new TGraphErrors(Theta_Bin_Final, Q2, NM_Bin, 0, NM_Bin_Err);
    gr3->SetMarkerColor(6);
    gr3->SetMarkerStyle(21);
    gr3->Draw("P");

	TLegend *t1 = new TLegend(0.65,0.60,0.80,0.94,"Counts per bin:");
	t1->SetBorderSize(0);
	t1->SetTextSize(0.03);
	t1->SetTextFont(22);
	t1->AddEntry(gr1,"Elastic","p");
	t1->AddEntry(gr2,"Quasi-Elastic","p");
	t1->AddEntry(gr3,"Moller","p");
	t1->Draw();
	c1->Print(Form("%s_Bining_%d.png",Target.Data(), E0));

    delete h3;	
    delete gr;	
    delete h2;	
    delete gr1;	
    delete gr2;	
    delete gr3;	
	delete c1;
	delete TE;
	delete TQ;
	delete TM;
}
