double Moller(const double E0, const double Theta_Degree){

	double P0 = sqrt( E0*E0 - Electron_Mass*Electron_Mass);
	double Theta_CM = 2.0 * atan( tan(Theta_Degree*PI/180) * sqrt( (P0+Electron_Mass)/(2.0*Electron_Mass) ) );
    double SinTh_CM = sin(Theta_CM);
    double CosTh_CM = cos(Theta_CM);

	double S = 2*Electron_Mass*(Electron_Mass + E0);//MeV2

	double XS = pow(ALPHA,2)/S * pow( 3+pow(CosTh_CM,2) ,2) / pow(SinTh_CM,4); //MeV^-2
    XS *= pow(HBARC,2);//in fm^2 (MeV2FM from MeV-1 to fm^1)

	return XS;
}

