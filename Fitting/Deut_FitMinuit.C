//////////////////////////////////////
//   Deut SOG Parameters Fitting    //
//     TMinuit Minimization         //
//    ---  Zhihong Ye 02/11/2015    //
//////////////////////////////////////

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
#include <TMinuit.h>
//#include <TMatrix.h>
/*}}}*/

/*}}}*/
using namespace::std;
using namespace::TMath;
#include "../../Constants.h"
static const Int_t N = 24;
static const Int_t iParaNum=4*N;
ofstream outlog("fitting.log");

void output(Double_t *inPar, Double_t *outPar,Double_t *outErr);
void LoadData();
void Init(double *kR, double *kAC, double *kAQ, double* kAM);
int GetFF_SOG(const double aQ_FM, double* aGE, double *aGQ, double *aGM, const double* R, const double* Q_E, const double* Q_Q, const double* Q_M);
double Deut_Elastic( double E0, double Theta, double *R, double* AC, double* AQ, double* AM);
void myFcn(Int_t& npar, Double_t* deriv, Double_t& f, Double_t *par, Int_t flag);
void DoMinuit(double *R, double* AC, double* AQ, double* AM);
void Test_Model(double *R, double* AC, double* AQ, double* AM);

vector<double> vE0, vTheta, vEp, vQ2, vXS, vXS_Err;

/*Main{{{*/
Int_t main()
{ 
	gStyle->SetOptFit(1);  
	gStyle->SetOptStat(0);

	double R[N] = {N*0};
	double AC[N] = {N*0};
	double AQ[N] = {N*0};
	double AM[N] = {N*0};

	//Initialize the SOG parameters
	Init(R, AC, AQ, AM);

	//Load the data points
	LoadData();

	//Test the XS model first
	Test_Model(R, AC, AQ, AM);

	/////////////////////////////
	//TMinuit Minimization
	/////////////////////////////
	DoMinuit(R, AC, AQ, AM);

	return 0;
}
/*}}}*/

/*DoMunuit{{{*/
void DoMinuit(double *R, double* AC, double* AQ, double* AM)
{
	cout << "      *************************************************" << endl; 
	cout << "      *          Minimization Fitting for Deuteron    *" << endl;
	cout << "      *              Elastic FF SOG Function          *" << endl;
	cout << "      *                 Z. Ye 02/11/2015              *" << endl;
	cout << "      *************************************************" << endl; 
	cout << endl;
	gSystem->Load("libMinuit");

	TMinuit *pMinuit = new TMinuit(iParaNum);  //initialize TMinuit with a maximum of iParaNum params
	pMinuit->SetFCN(myFcn);

	Double_t arglist[10];
	Int_t ierflg = 0;

	/// set print level
	/*  SET PRIntout  <level>
		Sets the print level, determining how much output will be
		produced. Allowed values and their meanings are displayed
		after a SHOw PRInt command, and are currently <level>=:
		[-1]  no output except from SHOW commands
		[0]  minimum output
		[1]  default value, normal output
		[2]  additional output giving intermediate results.
		[3]  maximum output, showing progress of minimizations. */
	arglist[0] = 1;
	pMinuit->mnexcm("SET PRIntout",arglist,1,ierflg);

	/*
	   SET NOWarnings
	   Supresses Minuit warning messages.
	   SET WARnings
	   Instructs Minuit to output warning messages when suspicious
	   conditions arise which may indicate unreliable results.
	   This is the default.
	   */
	arglist[0] = 1;
	pMinuit->mnexcm("SET NOWarnings",arglist,0,ierflg);

	//Set Name of Parameters
	TString namePar[iParaNum];
	for(int i=0;i<N;i++){
		namePar[i] = Form("R%d",i);
		namePar[N+i] = Form("AC%d",i);
		namePar[N+N+i] = Form("AQ%d",i);
		namePar[N+N+N+i] = Form("AM%d",i);
	}

	// Double_t inPar[iParaNum]={Target->f0,   
	// 			    Target->B,   
	// 			    Target->alpha,
	// 			    Target->a,
	// 			    Target->b};	
	Double_t inPar[iParaNum]={iParaNum*0.};   
	for(int i=0;i<N;i++){
		inPar[i] = R[i];
		inPar[N+i] = AC[i];
		inPar[N+N+i] = AQ[i];
		inPar[N+N+N+i] = AM[i];
	}

	//Set Stepsize of Parameters
	Double_t step[iParaNum]={ iParaNum*0.00000001};

	//Set Min of Parameters, value==0 for No-Bound
	Double_t minVal[iParaNum]={ iParaNum*0.0};

	//Set Max of Parameters, value==0 for No-Bound
	Double_t maxVal[iParaNum]={ iParaNum*0.0};

	//Initializing Parameters
	for(int ii=0;ii<iParaNum;ii++)
	{
		step[ii] = 1e-7;
		minVal[ii] = 0.0;
		maxVal[ii] = 0.0;

		pMinuit->DefineParameter(ii,namePar[ii], 
				inPar[ii], step[ii], minVal[ii], maxVal[ii]);
	}

	//Fix parameters
	//arglist[0] = 3;   pMinuit->mnexcm("FIX ", arglist ,1,ierflg);
	//	for(int i=2*N;i<3*N;i++)
	//		arglist[0] = 3;   pMinuit->mnexcm("FIX ", arglist ,1,ierflg);
	//	for(int i=3*N;i<4*N;i++)
	//		arglist[0] = 3;   pMinuit->mnexcm("FIX ", arglist ,1,ierflg);


	arglist[0] = 1;
	pMinuit->mnexcm("SET ERR", arglist,1,ierflg);

	/*
	//SET LIMits  0->[parno]  1->[lolim]  2->[uplim]
	//	  arglist[0] = 0;	   //this means set all parameters with the same limit
	//	  arglist[0] = 1;	   //this means set 1st parameters with the specified limit
	//	  arglist[0] = 2;	   //this means set 2nd parameters with the sepecified limit
	arglist[0] = 0;
	arglist[1] = 0.;
	arglist[2] = 0.;
	pMinuit->mnexcm("SET LIMits",arglist,3,ierflg);
	*/

	//Set Iteration Times and Use Migrad to start minimization
	arglist[0] = 500;
	arglist[1] = 0.00001;
	pMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

	//Get Results
	Double_t outPar[iParaNum]={iParaNum*0.};   
	Double_t err[iParaNum]={iParaNum*0.};  
	for(int ii = 0; ii < iParaNum; ii++)
	{
		pMinuit->GetParameter(ii,outPar[ii],err[ii]);    
	}
	//Put the results into a file
	output(inPar, outPar, err);

}
/*}}}*/

/*myFcn{{{*/
void myFcn(Int_t& npar, Double_t* deriv, Double_t& f, Double_t *par, Int_t flag)
{
	//Set Parameters
	Double_t R[iParaNum];
	Double_t AC[iParaNum];
	Double_t AQ[iParaNum];
	Double_t AM[iParaNum];

	for(int i = 0; i < N; i++){
		R[i]  = par[i];      
		AC[i] = par[N+i];   
		AQ[i] = par[N+N+i];
		AM[i] = par[N+N+N+i]; 
	}

	Double_t kElastic = 0.0;
	Double_t aDelta = 0.0, aSum=0.0;
	for ( unsigned int ii=0; ii< vE0.size(); ii++ )
	{
		kElastic = Deut_Elastic(vE0.at(ii), vTheta.at(ii), R, AC, AQ, AM);
		aDelta = (vXS.at(ii)-kElastic) /vXS.at(ii);
		aSum+=aDelta*aDelta;
	}

	f =  sqrt(aSum)/vE0.size();
	//f = aSum;
	//   cerr << " >>> Chisq = " << f <<endl;
}
/*}}}*/

/*Test_Model{{{*/
void Test_Model(double *R, double* AC, double* AQ, double* AM){

	Double_t kElastic = 0.0;
	for ( unsigned int ii=0; ii< 10; ii++ )
	{
		kElastic = Deut_Elastic(vE0.at(ii), vTheta.at(ii), R, AC, AQ, AM);
		outlog<<Form("---- E0 = %f, Theta = %f, XS_EX = %e,  XS_Calc = %e",vE0.at(ii),vTheta.at(ii), vXS.at(ii), kElastic)<<endl;	
	}
}
/*}}}*/

/*Deut_Elastic{{{*/
double Deut_Elastic( double E0, double Theta, double *R, double* AC, double* AQ, double* AM){

	double Theta_Rad = Theta *PI/180.0;//Rad
	double SinThSQ= pow(sin(Theta_Rad/2),2);
	double TanTh = tan(Theta_Rad/2.0);
	double ETA= 1./(1.+ 2.*E0/Deut_Mass*SinThSQ);
	double Ep= ETA*E0;
	double Q2= 4.*E0*Ep*SinThSQ;
	double TAO = Q2 / (4.0*Deut_Mass*Deut_Mass);

	double Q2_FM  = Q2/pow(HBARC,2);
	double Q_FM  = sqrt(Q2_FM);
	double G_C = 0.0, G_Q = 0.0, G_M = 0.0;

	int err = GetFF_SOG(Q_FM, &G_C, &G_Q, &G_M, R, AC, AQ, AM);
	if(err<0)
		return 1e-33;

	double GC2 = pow(G_C,2);
	double GM2 = pow(G_M,2);
	double GQ2 = pow(G_Q,2);
	double A   = GC2 + (2.0*TAO/3.0) * GM2 + (8.0*TAO*TAO/9.0) * GQ2;
	double B   = (4.0*TAO/3.0) * (1.0+TAO) * GM2;
	double FF  = A + B*TanTh*TanTh;

	//---------------------------------------------------------------------
	//    Get Mott cross section and calculate unpolarized Deuteron
	//    elastic cross section.
	//---------------------------------------------------------------------
	double Sig_Mott = 4.0*pow(HBARC,2)*pow((ALPHA*cos(Theta/2.0)*Ep)/Q2,2);//in fm^2/sr
	double REC = Ep/E0;                   //recoil factor
	double Sig_Deut = Sig_Mott * FF * REC;     // deuteron cross section in fm^2/sr
	////Debug
	//cerr<<Form("-- Sig_Mot=%f, GC=%f, GM=%f, GQ= %f, A=%f, B=%f, FF=%f, TAO=%f",
	//			Sig_Mott, GC2, GM2, GQ2, A, B, FF,TAO)<<endl;

	double Q_FM_Max = 0.88183;//PRAD Max Q_FM, at E0=2200MeV,Theta=4 Deg
	double CC = 0.01*(1- Q_FM/Q_FM_Max);
	if(CC<0.0) CC = 0.0;
	Sig_Deut /= (1+CC);//The reverse process, correct for the CC effect

	return Sig_Deut*FM2NB;//return nb/sr
}
/*}}}*/

/*SOG{{{*/
int GetFF_SOG(const double aQ_FM, double* aGE, double *aGQ, double *aGM, const double* R, const double* Q_E, const double* Q_Q, const double* Q_M){

	const double Gamma = sqrt(2./3.)*0.600; //fm

	//I. Sick, Prog. Part. Nucl. Phys. 47(2001) page#268
	double Exp_Term = TMath::Exp(-1./4.* aQ_FM*aQ_FM * Gamma*Gamma);
	double Sum_GE = Q_E[0], Sum_GQ = Q_Q[0], Sum_GM = Q_M[0];

	for(int i=1;i<N;i++){
		Sum_GE += (Q_E[i]/(1.+2.*R[i]*R[i]/Gamma/Gamma))
			* ( cos(aQ_FM*R[i]) + (2.*R[i]*R[i]/Gamma/Gamma)*sin(aQ_FM*R[i])/aQ_FM/R[i] );

		Sum_GQ += (Q_Q[i]/(1.+2.*R[i]*R[i]/Gamma/Gamma))
			* ( cos(aQ_FM*R[i]) + (2.*R[i]*R[i]/Gamma/Gamma)*sin(aQ_FM*R[i])/aQ_FM/R[i] );

		Sum_GM += (Q_M[i]/(1.+2.*R[i]*R[i]/Gamma/Gamma))
			* ( cos(aQ_FM*R[i]) + (2.*R[i]*R[i]/Gamma/Gamma)*sin(aQ_FM*R[i])/aQ_FM/R[i] );
	}

	aGE[0] = Exp_Term * Sum_GE;//I. Sick, Prog. Part. Nucl. Phys. 47(2001) page#267
	//  aGQ[0] = Exp_Term * Sum_GQ *( (Deut_Mass*Deut_Mass)   * Deut_Quad_Moment);//I. Sick, Prog. Part. Nucl. Phys. 47(2001) page#267 
	//  aGM[0] = Exp_Term * Sum_GM *( (Deut_Mass/Proton_Mass) * Deut_Mag_Moment);//I. Sick, Prog. Part. Nucl. Phys. 47(2001) page#267
	aGQ[0] = Exp_Term * Sum_GQ *(25.83);// experimental value,(Md/Mp*ud = 1.714)
	aGM[0] = Exp_Term * Sum_GM *(1.714 );// experimental value, (Md*Md*Qd = 25.83)

	if(aGE[0]>1e-13 && aGQ[0]>1e-33 && aGM[0]>1e-33)
		return 1;
	else
		return -1;
}
/*}}}*/

/*Init{{{*/
void Init(double *kR, double *kAC, double *kAQ, double* kAM){

	/*SOG Parameterization by Ingo Sick 2012 result{{{*/
	const double R[N] = { 0.000000, 0.400000, 0.800000, 1.200000, 1.600000,
		2.000000, 2.400000, 2.800000, 3.200000, 3.600000,
		4.000000, 4.400000, 4.800000, 5.200000, 5.600000,
		6.000000, 6.500000, 7.000000, 7.500000, 8.000000,
		8.500000, 9.000000, 9.500000, 10.00000};
	/*	const double Q_C[N] = { 0.001586, 0.153416, 0.177804, 0.198969, 0.145431,
		0.097231, 0.067972, 0.050103, 0.025602, 0.025764,
		0.017142, 0.012376, 0.008021, 0.006033, 0.003642,
		0.003163, 0.002188, 0.001342, 0.000851, 0.000538,
		0.000333, 0.000217, 0.000126, 0.000100};
		*/
	const double Q_Q[N] = { 0.046513, 0.138087, 0.181425, 0.174011, 0.139929,
		0.091150, 0.070100, 0.047832, 0.031397, 0.025241,
		0.016786, 0.012042, 0.007777, 0.005819, 0.003497,
		0.003023, 0.002082, 0.001270, 0.000802, 0.000503,
		0.000311, 0.000200, 0.000116, 0.000090};
	const double Q_M[N] = { 0.043842, 0.159966, 0.182563, 0.177163, 0.120109,
		0.103085, 0.067460, 0.038359, 0.031553, 0.023661,
		0.015861, 0.011489, 0.007484, 0.005640, 0.003418,
		0.002971, 0.002063, 0.001265, 0.000805, 0.000508,
		0.000316, 0.000204, 0.000121, 0.000093};
	/*}}}*/  

	const double Q_C[N] = { 0.000000, 0.153416, 0.177804, 0.198969, 0.145431,
		0.097231, 0.067972, 0.050103, 0.025602, 0.025764,
		0.017142, 0.012376, 0.008021, 0.006033, 0.003642,
		0.003163, 0.002188, 0.001342, 0.000851, 0.000538,
		0.000333, 0.000217, 0.000126, 0.000100};

	for(int i=0;i<N;i++){
		kR[i] = R[i];
		kAC[i] = Q_C[i];
		kAQ[i] = Q_Q[i];
		kAM[i] = Q_M[i];
	}
}

/*}}}*/  

/*LoadData{{{*/
void LoadData(){

	ifstream infile1("Deut_Bin_1100.dat");
	ifstream infile2("Deut_Bin_2200.dat");
	ifstream infile3("Deut_Bin_3300.dat");
	ifstream infile4("Deut_Bin_4400.dat");
	TString Com;
	infile1 >> Com >> Com >> Com >> Com >> Com >> Com;
	infile1 >> Com >> Com >> Com >> Com >> Com >> Com;
	infile2 >> Com >> Com >> Com >> Com >> Com >> Com;
	infile2 >> Com >> Com >> Com >> Com >> Com >> Com;
	infile3 >> Com >> Com >> Com >> Com >> Com >> Com;
	infile3 >> Com >> Com >> Com >> Com >> Com >> Com;
	infile4 >> Com >> Com >> Com >> Com >> Com >> Com;
	infile4 >> Com >> Com >> Com >> Com >> Com >> Com;

	double E0,Ep,Theta, Q2, XS, XS_Err;
	while(!(infile1.eof()) ){
		infile1 >> E0  >>Theta>> Ep  >> Q2  >> XS  >> XS_Err;
		vE0.push_back(E0);
		vEp.push_back(Ep);
		vTheta.push_back(Theta);
		vQ2.push_back(Q2);
		vXS.push_back(XS);
		vXS_Err.push_back(XS_Err);
	}
	infile1.close();
	while(!(infile2.eof()) ){
		infile2 >> E0  >>Theta>> Ep  >> Q2  >> XS  >> XS_Err;
		vE0.push_back(E0);
		vEp.push_back(Ep);
		vTheta.push_back(Theta);
		vQ2.push_back(Q2);
		vXS.push_back(XS);
		vXS_Err.push_back(XS_Err);
	}
	infile2.close();
	while(!(infile3.eof()) ){
		infile3 >> E0  >>Theta>> Ep  >> Q2  >> XS  >> XS_Err;
		vE0.push_back(E0);
		vEp.push_back(Ep);
		vTheta.push_back(Theta);
		vQ2.push_back(Q2);
		vXS.push_back(XS);
		vXS_Err.push_back(XS_Err);
	}
	infile3.close();
	while(!(infile4.eof()) ){
		infile4 >> E0  >>Theta>> Ep  >> Q2  >> XS  >> XS_Err;
		vE0.push_back(E0);
		vEp.push_back(Ep);
		vTheta.push_back(Theta);
		vQ2.push_back(Q2);
		vXS.push_back(XS);
		vXS_Err.push_back(XS_Err);
	}
	infile4.close();

	//for(int i=0;i<vQ2.size();i++)
	//	cerr<<Form("--- Q2 = %f,  XS = %e", vQ2.at(i), vXS.at(i))<<endl;
}
/*}}}*/

/*output{{{*/
//Output the results into a file
void output(Double_t *inPar, Double_t *outPar,Double_t *outErr)
{

	ofstream outfile("fitting.dat");
	Double_t R_Old[iParaNum],  R_Out[iParaNum], R_Out_Err[iParaNum];	
	Double_t AC_Old[iParaNum], AC_Out[iParaNum], AC_Out_Err[iParaNum];	
	Double_t AQ_Old[iParaNum], AQ_Out[iParaNum], AQ_Out_Err[iParaNum];	
	Double_t AM_Old[iParaNum], AM_Out[iParaNum], AM_Out_Err[iParaNum];	

	outfile<<"----------------------------------------------"<<endl;
	outfile<<Form("%10s %10s %10s","R_Old","R","R_Err")<<endl;
	outfile<<"----------------------------------------------"<<endl;
	for(int i = 0; i < N; i++){
		R_Old[i]  = inPar[i];       R_Out[i]  = outPar[i];       R_Out_Err[i]  = outErr[i];

		outfile<<Form("%10.4e %10.4e %10.4e",
				R_Old[i], R_Out[i],  R_Out_Err[i])<<endl;
	}

	outfile<<endl;
	outfile<<"----------------------------------------------"<<endl;
	outfile<<Form("%10s %10s %10s","AC_Old", "AC","AC_Err")<<endl;
	outfile<<"----------------------------------------------"<<endl;
	for(int i = 0; i < N; i++){
		AC_Old[i] = inPar[N+i];     AC_Out[i] = outPar[N+i];     AC_Out_Err[i] = outErr[N+i];
		outfile<<Form("%10.4e %10.4e %10.4e", AC_Old[i], AC_Out[i], AC_Out_Err[i])<<endl;
	}

	outfile<<endl;
	outfile<<"----------------------------------------------"<<endl;
	outfile<<Form("%10s %10s %10s","AQ_Old", "AQ","AQ_Err")<<endl;
	outfile<<"----------------------------------------------"<<endl;
	for(int i = 0; i < N; i++){
		AQ_Old[i] = inPar[N+N+i];   AQ_Out[i] = outPar[N+N+i];   AQ_Out_Err[i] = outErr[N+N+i];
		outfile<<Form("%10.4e %10.4e %10.4e", AQ_Old[i], AQ_Out[i], AQ_Out_Err[i])<<endl;
	}

	outfile<<endl;
	outfile<<"----------------------------------------------"<<endl;
	outfile<<Form("%10s %10s %10s","AM_Old", "AM","AM_Err")<<endl;
	outfile<<"----------------------------------------------"<<endl;
	for(int i = 0; i < N; i++){
		AM_Old[i] = inPar[N+N+N+i]; AM_Out[i] = outPar[N+N+N+i]; AM_Out_Err[i] = outErr[N+N+N+i];

		outfile<<Form("%10.4e %10.4e %10.4e", AM_Old[i], AM_Out[i], AM_Out_Err[i])<<endl;
	}


	outfile<<endl;
	outfile<<"----------------------------------------------"<<endl;
	outfile<<Form("Save the results as C/C++ arrays:")<<endl;
	outfile<<"----------------------------------------------"<<endl;

	outfile <<Form("R[%d] = { ",N);
	for(int i = 0; i < N-1; i++){
		outfile<<Form("%10.4e ,", R_Out[i]);
		if(!(i%6)) outfile<<endl;
	}
	outfile<<Form("%10.4e", R_Out[N-1]) <<"};"<<endl;

	outfile <<Form("A_C[%d] = { ",N);
	for(int i = 0; i < N-1; i++){
		outfile<<Form("%10.4e ,", AC_Out[i]);
		if(!(i%6)) outfile<<endl;
	}
	outfile<<Form("%10.4e", AC_Out[N-1]) <<"};"<<endl;

	outfile <<Form("A_Q[%d] = { ",N);
	for(int i = 0; i < N-1; i++){
		outfile<<Form("%10.4e ,", AQ_Out[i]);
		if(!(i%6)) outfile<<endl;
	}
	outfile<<Form("%10.4e", AQ_Out[N-1]) <<"};"<<endl;

	outfile <<Form("A_M[%d] = { ",N);
	for(int i = 0; i < N-1; i++){
		outfile<<Form("%10.4e ,", AM_Out[i]);
		if(!(i%6)) outfile<<endl;
	}
	outfile<<Form("%10.4e", AM_Out[N-1]) <<"};"<<endl;

}
/*}}}*/
