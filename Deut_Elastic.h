
int GetFF_SOG(const double aQ_FM, double* aGE, double *aGQ, double *aGM);
double RINTERPQ(const double x0, double* x, double* y, const int N);
int GetFF(double Q2_New, double *FC_New, double *FQ_New, double *FM_New);
int GetFF_RInterPQ(double Q2_New, double *FC_New, double *FQ_New, double *FM_New);

/*double Deut_Elastic( double E0, double Theta){{{*/
double Deut_Elastic( double E0,double Theta_Degree){
	double TanThSQ,SinThSQ,TAO,Q2_FM,Ep, Q2;
	double A,B,FF,G_C,G_M,G_Q,GC2,GM2,GQ2;
	double Sig_Mott,REC;
	bool bCoulomb_Distort = kTRUE;//Add Coulomb Distortion when simualte "Exp" data
	//bool bCoulomb_Distort = kFALSE; //When comparing with model or extract, we need to correct it back.
	
	//	int model = 0; cerr<<" Which model? (1->InterP, 2->SoG)"; cin >> model;
	int model = 2; 

	//---------------------------------------------------------------------
	//    Calculate some kinematical quantities
	//---------------------------------------------------------------------
	double Theta = Theta_Degree * PI/180.0;
	TanThSQ = tan(Theta/2.0);
	SinThSQ = pow(sin(Theta/2.0),2);
	Ep= E0/(1.+ 2.*E0/Deut_Mass*SinThSQ); //MeV
	Q2 =  4.0 * E0 * Ep * SinThSQ; //MeV2
	TAO  = Q2/(4.0 * Deut_Mass*Deut_Mass);

	//---------------------------------------------------------------------
	//    Get deuteron form factors
	//---------------------------------------------------------------------
	Q2_FM  = Q2/pow(HBARC,2);
	double Q_FM  = sqrt(Q2_FM);
	int err = 0;
	if(model==2){
		err = GetFF_SOG(Q_FM, &G_C, &G_Q, &G_M);
	//	cerr<<Form("--   SOG: GC=%f, GQ=%f, GM=%f", G_C, G_Q, G_M)<<endl;
	}
		
	if(model==1){
		if(Q2_FM>=0.01&&Q2_FM<=100.0)
			err = GetFF(Q2_FM, &G_C, &G_Q, &G_M);
		else
			err = GetFF_SOG(Q_FM, &G_C, &G_Q, &G_M);

	//		cerr<<Form("--InterP: GC=%f, GQ=%f, GM=%f", G_C, G_Q, G_M)<<endl;
	}

	if(err<0){
		cerr<<"**** Something wrong with getting the Deuteron FF"<<endl;
		return -1;
	}

	GC2 = pow(G_C,2);
	GM2 = pow(G_M,2);
	GQ2 = pow(G_Q,2);
	A   = GC2 + (2.0*TAO/3.0) * GM2 + (8.0*TAO*TAO/9.0) * GQ2;
	B   = (4.0*TAO/3.0) * (1.0+TAO) * GM2;
	FF  = A + B*TanThSQ*TanThSQ;

	//---------------------------------------------------------------------
	//    Get Mott cross section and calculate unpolarized Deuteron
	//    elastic cross section.
	//---------------------------------------------------------------------
	Sig_Mott = 4.0*pow(HBARC,2)*pow((ALPHA*cos(Theta/2.0)*Ep)/Q2,2);//in fm^2/sr
	REC = Ep/E0;                   //recoil factor
	double Sig_Deut = Sig_Mott * FF * REC;     // deuteron cross section in fm^2/sr
	////Debug
//	cerr<<Form("-- Sig_Mot=%f, GC=%f, GM=%f, GQ= %f, A=%f, B=%f, FF=%f, TAO=%f",
//			Sig_Mott, GC2, GM2, GQ2, A, B, FF,TAO)<<endl;

	if(bCoulomb_Distort){
		double Q_FM_Max = 0.88183;//PRAD Max Q_FM, at E0=2200MeV,Theta=4 Deg
		double R = 0.01*(1- Q_FM/Q_FM_Max);
		if(R<0.0) R = 0.0;
		Sig_Deut *= (1+R);
	}

	return Sig_Deut;
}
/*}}}*/

/*SOG{{{*/
int GetFF_SOG(const double aQ_FM, double* aGE, double *aGQ, double *aGM){

  /*SOG Parameterization by Ingo Sick 2012 result{{{*/
  const int N = 24;
  const double Gamma = sqrt(2./3.)*0.600; //fm
  const double R[N] = { 0.000000, 0.400000, 0.800000, 1.200000, 1.600000,
			2.000000, 2.400000, 2.800000, 3.200000, 3.600000,
			4.000000, 4.400000, 4.800000, 5.200000, 5.600000,
			6.000000, 6.500000, 7.000000, 7.500000, 8.000000,
			8.500000, 9.000000, 9.500000, 10.00000};
  const double Q_E[N] = { 0.001586, 0.153416, 0.177804, 0.198969, 0.145431,
			  0.097231, 0.067972, 0.050103, 0.025602, 0.025764,
			  0.017142, 0.012376, 0.008021, 0.006033, 0.003642,
			  0.003163, 0.002188, 0.001342, 0.000851, 0.000538,
			  0.000333, 0.000217, 0.000126, 0.000100};
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

/*Quadratic Interpolation Routine{{{*/

//-----------------------------------------------------------------------------
//    Quadratic Interpolation Routine:
//
//      Calculates Y(Q2_New) given array Y(X) assuming a quadratic
//      dependence:  Y = AX^2 + BX + C
//
//      This routine assumes the X values are arranged in
//      ascending order but does not assume a uniform X spacing
//      of the array.
//-----------------------------------------------------------------------------
int GetFF(double Q2_New, double *FC_New, double *FQ_New, double *FM_New){

	double FC[100],FQ[100],FM[100],Q2_FC[100],Q2_FM[100],Q2_FQ[100];
	//cerr<<"--- What Deuteron FF model? (1->IA, 2->IAMEC, 3->RSC, 4->RSCMEC)  ";
	//int Model; cin >> Model;
	int Model = 4;

	/* Obtain Form Factors{{{*/
	for(int i=0;i<100;i++){
		FC[i] = 0.0;  FQ[i] = 0.0;  FM[i] = 0.0;
		Q2_FC[i] = 0.0; Q2_FQ[i] = 0.0; Q2_FM[i] = 0.0;
	}
	ifstream FC_File; ifstream FQ_File;	ifstream FM_File;
	/* IA{{{*/
	if(Model==1){
		FC_File.open("./DAT/ia.fc");
		FQ_File.open("./DAT/ia.fq");
		FM_File.open("./DAT/ia.fm");
	}
	/*}}}*/
	/* IAMEC{{{*/
	if(Model==2){
		FC_File.open("./DAT/iamec.fc");
		FQ_File.open("./DAT/iamec.fq");
		FM_File.open("./DAT/iamec.fm");
	}
	/*}}}*/
	/* RSC{{{*/
	if(Model==3){
		FC_File.open("./DAT/rsc.fc");
		FQ_File.open("./DAT/rsc.fq");
		FM_File.open("./DAT/rsc.fm");
	}
	/*}}}*/
	/* RSCMEC{{{*/
	if(Model==4){
		FC_File.open("./DAT/rscmec.fc");
		FQ_File.open("./DAT/rscmec.fq");
		FM_File.open("./DAT/rscmec.fm");
	}
	/*}}}*/

	int N_FC = 0;	
	while(!(FC_File.eof())){
		FC_File >> Q2_FC[N_FC] >> FC[N_FC]; 
		//cerr<<Form("--- FC:  Q2 = %f, FC = %f", Q2_FC[N_FC], FC[N_FC])<<endl;	
		N_FC++;
	} 
	N_FC--;
	int N_FQ = 0;	
	while(!(FQ_File.eof())){
		FQ_File >> Q2_FQ[N_FQ] >> FQ[N_FQ]; 
		//cerr<<Form("--- FQ:  Q2 = %f, FQ = %f", Q2_FQ[N_FQ], FQ[N_FQ])<<endl;	
		N_FQ++;
	} 
	N_FQ--;
	int N_FM = 0;	
	while(!(FM_File.eof())){
		FM_File >> Q2_FM[N_FM] >> FM[N_FM];
		//cerr<<Form("--- FM:  Q2 = %f, FM = %f", Q2_FM[N_FM], FM[N_FM])<<endl;	
		N_FM++;
	} 
	N_FM--;
	FC_File.close(); FQ_File.close(); FM_File.close();
	/*}}}*/
	
	FC_New[0] = RINTERPQ(Q2_New, Q2_FC, FC, N_FC);
	FQ_New[0] = RINTERPQ(Q2_New, Q2_FQ, FQ, N_FQ);
	FQ_New[0] *= pow(Deut_Mass/HBARC,2);
	FM_New[0] = RINTERPQ(Q2_New, Q2_FM, FM, N_FM);

	return 0;
}

/*RINTERPQ{{{*/
double RINTERPQ(const double x0, double* x, double* y, const int N){

	int i,j,k;

	if(x0 < x[0] || x0 > x[N-1]){
		cerr<<"*** Extrapolating outside range for X = "<< x0<<Form("out (%f ~ %f)",x[0],x[N-1])<<endl;
		return -111111; 
	}
	else{
		if(x0<=x[0]) {
			i = 1; j = 2; k = 3;
		}
		else if(x0>=x[N-2]){
			i = N-2; j = N-1; k = N-0;
		}
		else{
			for(int l=1;l<=N-2;l++){
				if(x0<=x[l]){
					i=l-1;
					j=l;
					k=l+1;
					break;
				}
			}
		}

		double Det = (x[j]-x[k])*(pow(x[i],2) - x[i]*(x[j]+x[k]) + x[j]*x[k]);
		double  A = ( y[i]*(x[j]-x[k]) - x[i]*(y[j]-y[k]) + y[j]*x[k]
				- x[j]*y[k] )/Det;
		double  B = ( pow(x[i],2)*(y[j]-y[k]) - y[i]*(pow(x[j],2)-pow(x[k],2))
				+ pow(x[j],2)*y[k] - pow(x[k],2)*y[j] )/Det;
		double  C = ( pow(x[i],2)*(x[j]*y[k]-x[k]*y[j])
				- x[i]*(pow(x[j],2)*y[k]-pow(x[k],2)*y[j])
				+ y[i]*(pow(x[j],2)*x[k]-pow(x[k],2)*x[j]) )/Det;

		double RINTERPQ = A*x0*x0 + B*x0 + C;

		return RINTERPQ;
	}
}
/*}}}*/
/*}}}*/

/*Old--Quadratic Interpolation Routine{{{*/

//-----------------------------------------------------------------------------
//    Quadratic Interpolation Routine:
//
//      Calculates Y(Q2_New) given array Y(X) assuming a quadratic
//      dependence:  Y = AX^2 + BX + C
//
//      This routine assumes the X values are arranged in
//      ascending order but does not assume a uniform X spacing
//      of the array.
//-----------------------------------------------------------------------------
int GetFF_RInterPQ(double Q2_New, double *FC_New, double *FQ_New, double *FM_New){

	double FC[100],FQ[100],FM[100],Q2_FC[100],Q2_FM[100],Q2_FQ[100];
	//cerr<<"--- What Deuteron FF model? (1->IA, 2->IAMEC, 3->RSC, 4->RSCMEC)  ";
	//int Model; cin >> Model;
	int Model = 4;

	/* Obtain Form Factors{{{*/
	for(int i=0;i<100;i++){
		FC[i] = 0.0;  FQ[i] = 0.0;  FM[i] = 0.0;
		Q2_FC[i] = 0.0; Q2_FQ[i] = 0.0; Q2_FM[i] = 0.0;
	}
	ifstream FC_File; ifstream FQ_File;	ifstream FM_File;
	/* IA{{{*/
	if(Model==1){
		FC_File.open("./DAT/ia.fc");
		FQ_File.open("./DAT/ia.fq");
		FM_File.open("./DAT/ia.fm");
	}
	/*}}}*/
	/* IAMEC{{{*/
	if(Model==2){
		FC_File.open("./DAT/iamec.fc");
		FQ_File.open("./DAT/iamec.fq");
		FM_File.open("./DAT/iamec.fm");
	}
	/*}}}*/
	/* RSC{{{*/
	if(Model==3){
		FC_File.open("./DAT/rsc.fc");
		FQ_File.open("./DAT/rsc.fq");
		FM_File.open("./DAT/rsc.fm");
	}
	/*}}}*/
	/* RSCMEC{{{*/
	if(Model==4){
		FC_File.open("./DAT/rscmec.fc");
		FQ_File.open("./DAT/rscmec.fq");
		FM_File.open("./DAT/rscmec.fm");
	}
	/*}}}*/

	int N_FC = 0;	
	while(!(FC_File.eof())){
		FC_File >> Q2_FC[N_FC] >> FC[N_FC]; 
		//cerr<<Form("--- FC:  Q2 = %f, FC = %f", Q2_FC[N_FC], FC[N_FC])<<endl;	
		N_FC++;
	} 
	N_FC--;
	int N_FQ = 0;	
	while(!(FQ_File.eof())){
		FQ_File >> Q2_FQ[N_FQ] >> FQ[N_FQ]; 
		//cerr<<Form("--- FQ:  Q2 = %f, FQ = %f", Q2_FQ[N_FQ], FQ[N_FQ])<<endl;	
		N_FQ++;
	} 
	N_FQ--;
	int N_FM = 0;	
	while(!(FM_File.eof())){
		FM_File >> Q2_FM[N_FM] >> FM[N_FM];
		//cerr<<Form("--- FM:  Q2 = %f, FM = %f", Q2_FM[N_FM], FM[N_FM])<<endl;	
		N_FM++;
	} 
	N_FM--;
	FC_File.close(); FQ_File.close(); FM_File.close();
	/*}}}*/

	int I1 = 0, I2 = 0, I3 = 0;
	double DET = 0., A = 0., B = 0., C = 0.;
	/*Get FC {{{*/
	if(Q2_New<Q2_FC[0] || Q2_New > Q2_FC[N_FC-1]){
		cerr<<"*** Extrapolating outside range for FC:  Q2 = "<< Q2_New<<endl;
		cerr<<"    Q2_Min = "<<Q2_FC[0]<<", Q2_Max = "<<Q2_FC[N_FC-1]<<endl;
		return -1;
	}

	if(Q2_New <= Q2_FC[0]){
		I1 = 1; I2 = 2; I3 = 3;
	}
	else if(Q2_New >= Q2_FC[N_FC-2]){
		I1 = N_FC-2; I2 = N_FC-1; I3 = N_FC;
	}
	else{
		for(int I=1;I<N_FC-2;I++){
			if(Q2_New <= Q2_FC[I]){ 
				I1 = I-1; I2 = I; I3 = I+1;
				break;
			}
		}
	}

	DET = (Q2_FC[I2]-Q2_FC[I3])*(Q2_FC[I1]*Q2_FC[I1] - Q2_FC[I1]*(Q2_FC[I2]+Q2_FC[I3]) + Q2_FC[I2]*Q2_FC[I3]);
	A = ( FC[I1]*(Q2_FC[I2]-Q2_FC[I3]) - Q2_FC[I1]*(FC[I2]-FC[I3]) + FC[I2]*Q2_FC[I3]
			- Q2_FC[I2]*FC[I3] )/DET;
	B = ( Q2_FC[I1]*Q2_FC[I1]*(FC[I2]-FC[I3]) - FC[I1]*(Q2_FC[I2]*Q2_FC[I2]-Q2_FC[I3]*Q2_FC[I3])
			+ Q2_FC[I2]*Q2_FC[I2]*FC[I3] - Q2_FC[I3]*Q2_FC[I3]*FC[I2] )/DET;
	C = ( Q2_FC[I1]*Q2_FC[I1]*(Q2_FC[I2]*FC[I3]-Q2_FC[I3]*FC[I2])
			- Q2_FC[I1]*(Q2_FC[I2]*Q2_FC[I2]*FC[I3]-Q2_FC[I3]*Q2_FC[I3]*FC[I2])
			+ FC[I1]*(Q2_FC[I2]*Q2_FC[I2]*Q2_FC[I3]-Q2_FC[I3]*Q2_FC[I3]*Q2_FC[I2]) )/DET;

	FC_New[0] = A*Q2_New*Q2_New + B*Q2_New + C;
	/*}}}*/ 

	/*Get FQ {{{*/
	if(Q2_New<Q2_FQ[0] || Q2_New > Q2_FQ[N_FQ-1]){
		cerr<<"*** Extrapolating outside range for FQ:  Q2 = "<< Q2_New<<endl;
		cerr<<"    Q2_Min = "<<Q2_FQ[0]<<", Q2_Max = "<<Q2_FQ[N_FQ-1]<<endl;
		return -1;
	}

	if(Q2_New <= Q2_FQ[0]){
		I1 = 1; I2 = 2; I3 = 3;
	}
	else if(Q2_New >= Q2_FQ[N_FQ-2]){
		I1 = N_FQ-2; I2 = N_FQ-1; I3 = N_FQ;
	}
	else{
		for(int I=1;I<N_FQ-2;I++){
			if(Q2_New <= Q2_FQ[I]){ 
				I1 = I-1; I2 = I; I3 = I+1;
				break;
			}
		}
	}

	DET = (Q2_FQ[I2]-Q2_FQ[I3])*(Q2_FQ[I1]*Q2_FQ[I1] - Q2_FQ[I1]*(Q2_FQ[I2]+Q2_FQ[I3]) + Q2_FQ[I2]*Q2_FQ[I3]);
	A = ( FQ[I1]*(Q2_FQ[I2]-Q2_FQ[I3]) - Q2_FQ[I1]*(FQ[I2]-FQ[I3]) + FQ[I2]*Q2_FQ[I3]
			- Q2_FQ[I2]*FQ[I3] )/DET;
	B = ( Q2_FQ[I1]*Q2_FQ[I1]*(FQ[I2]-FQ[I3]) - FQ[I1]*(Q2_FQ[I2]*Q2_FQ[I2]-Q2_FQ[I3]*Q2_FQ[I3])
			+ Q2_FQ[I2]*Q2_FQ[I2]*FQ[I3] - Q2_FQ[I3]*Q2_FQ[I3]*FQ[I2] )/DET;
	C = ( Q2_FQ[I1]*Q2_FQ[I1]*(Q2_FQ[I2]*FQ[I3]-Q2_FQ[I3]*FQ[I2])
			- Q2_FQ[I1]*(Q2_FQ[I2]*Q2_FQ[I2]*FQ[I3]-Q2_FQ[I3]*Q2_FQ[I3]*FQ[I2])
			+ FQ[I1]*(Q2_FQ[I2]*Q2_FQ[I2]*Q2_FQ[I3]-Q2_FQ[I3]*Q2_FQ[I3]*Q2_FQ[I2]) )/DET;

	FQ_New[0] = A*Q2_New*Q2_New + B*Q2_New + C;
	FQ_New[0] *= pow(Deut_Mass/HBARC,2);
	/*}}}*/ 

	/*Get FM {{{*/
	if(Q2_New<Q2_FM[0] || Q2_New > Q2_FM[N_FM-1]){
		cerr<<"*** Extrapolating outside range for FM:  Q2 = "<< Q2_New<<endl;
		cerr<<"    Q2_Min = "<<Q2_FM[0]<<", Q2_Max = "<<Q2_FM[N_FM-1]<<endl;
		return -1;
	}

	if(Q2_New <= Q2_FM[0]){
		I1 = 1; I2 = 2; I3 = 3;
	}
	else if(Q2_New >= Q2_FM[N_FM-2]){
		I1 = N_FM-2; I2 = N_FM-1; I3 = N_FM;
	}
	else{
		for(int I=1;I<N_FM-2;I++){
			if(Q2_New <= Q2_FM[I]){ 
				I1 = I-1; I2 = I; I3 = I+1;
				break;
			}
		}
	}

	DET = (Q2_FM[I2]-Q2_FM[I3])*(Q2_FM[I1]*Q2_FM[I1] - Q2_FM[I1]*(Q2_FM[I2]+Q2_FM[I3]) + Q2_FM[I2]*Q2_FM[I3]);
	A = ( FM[I1]*(Q2_FM[I2]-Q2_FM[I3]) - Q2_FM[I1]*(FM[I2]-FM[I3]) + FM[I2]*Q2_FM[I3]
			- Q2_FM[I2]*FM[I3] )/DET;
	B = ( Q2_FM[I1]*Q2_FM[I1]*(FM[I2]-FM[I3]) - FM[I1]*(Q2_FM[I2]*Q2_FM[I2]-Q2_FM[I3]*Q2_FM[I3])
			+ Q2_FM[I2]*Q2_FM[I2]*FM[I3] - Q2_FM[I3]*Q2_FM[I3]*FM[I2] )/DET;
	C = ( Q2_FM[I1]*Q2_FM[I1]*(Q2_FM[I2]*FM[I3]-Q2_FM[I3]*FM[I2])
			- Q2_FM[I1]*(Q2_FM[I2]*Q2_FM[I2]*FM[I3]-Q2_FM[I3]*Q2_FM[I3]*FM[I2])
			+ FM[I1]*(Q2_FM[I2]*Q2_FM[I2]*Q2_FM[I3]-Q2_FM[I3]*Q2_FM[I3]*Q2_FM[I2]) )/DET;

	FM_New[0] = A*Q2_New*Q2_New + B*Q2_New + C;
	/*}}}*/ 

	return 0;
}
/*}}}*/
