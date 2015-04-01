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
//      Update: Convert into C++ by Zhihong Ye,   Nov/15/ 2014
//---------------------------------------------------------------------

//                   MeV           Degree   --> returning fm^2/sr
double He4_Elastic(double E0, double Theta_Degree){
	double Sig_Mott,CosThSQ,SinThSQ,TanThSQ;
	double Ep,ETA,FACC,FACM,FACS,FCH,FMAG,Q,Q4SQ;
	double SIG,TAU,Qeff;

	//---------------------------------------------------------------------
	//    Calculation of Kinematic Quantities
	//---------------------------------------------------------------------
    const int A = 4, Z = 2;
	double Theta = Theta_Degree * PI/180.;
	SinThSQ= pow(sin(Theta/2),2);
	CosThSQ= pow(cos(Theta/2),2);
	TanThSQ= SinThSQ/CosThSQ;
	ETA= 1./(1.+ 2.*E0/He4_Mass*SinThSQ);
	Ep= ETA*E0;
	Q4SQ= 4.*E0*Ep*SinThSQ;
	Q=sqrt(Q4SQ)/HBARC;           //q4 in fm**-1
	Sig_Mott= 10000.*pow((HBARC*ALPHA/2./E0),2)* CosThSQ/pow(SinThSQ,2);

	//---------------------------------------------------------------------
	//    Simple Calculation Of The Effective Q  
	//---------------------------------------------------------------------
	Qeff=Q*( 1. + 1.5*Z*ALPHA*HBARC/E0/(1.12*pow(A,0.33333)));

	//---------------------------------------------------------------------
	//    Calculation Of The Elastic Form Factors  
	//    C.R.Ottermann et al., Nucl. Phys A436(1985)688-698.   
	//---------------------------------------------------------------------
	const double FactA=0.316;
	const double FactB=0.675;
	FMAG=0.0; 
	FCH=(1-pow(FactA*FactA*Q*Q,6))*exp(-1*FactB*FactB*Q*Q);

	//---------------------------------------------------------------------
	//    Approxiamte Calculation The DWBA Cross Section  
	//---------------------------------------------------------------------
	TAU= Q4SQ/(4.*He4_Mass*He4_Mass);
	FACC= 1./((1.+TAU));
	double Mag_MomentA= 3./Z*He4_Mag_Moment;
	FACM= TAU*pow(Mag_MomentA,2)* (FACC+ 2.*TanThSQ);
	FACS= FCH*FCH*FACC+ FMAG*FMAG*FACM;
	SIG= Z*Z*Sig_Mott*ETA*FACS; // microbarns/sr'

	double Sig_He4 =SIG/10000.;      //He4 cross section in fm^2/sr
	return Sig_He4;
}

double GetFF_He4(double aQ2){

	//---------------------------------------------------------------------
	//    Calculation Of The Elastic Form Factors  
	//    C.R.Ottermann et al., Nucl. Phys A436(1985)688-698.   
	//---------------------------------------------------------------------
	const double FactA=0.316;
	const double FactB=0.675;
	double Q=sqrt(aQ2)/HBARC;           //q4 in fm**-1
	double FCH=(1-pow(FactA*FactA*Q*Q,6))*exp(-1*FactB*FactB*Q*Q);

	return FCH;
}

