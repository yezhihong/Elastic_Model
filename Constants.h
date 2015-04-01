//---------------------------------------------------------------------
//    Constants Used In The Calculation
//---------------------------------------------------------------------
const double PI=3.1415926;    
const double Deg2Rad = 3.1415926/180.0;
const double GeV2MeV = 1000.0;

const double HBARC=197.327;          // Conversion constant from MeV to fm^-1, 1*fm^-1=197.327 MeV, 1*MeV = 1./197.327*fm-1, 1*MeV^-1 = 197.327*fm
const double H_C=1239.841;   
const double ALPHA=1./137.035989;    //Fine Struction Constant 
const double Electron_Mass=0.511;    //MeV
const double He4_Mass=3727.40841;         //4He Nuclear Mass of 3He           
const double Deut_Mass = 1875.630;       //Deuteron mass in MeV
const double Proton_Mass = 938.272;       //Proton mass in MeV
const double G_F   = 1.16637E-11;    //Fermi constant

const double MeV2FM = 1/HBARC; //from MeV to fm-1

const double He4_Mag_Moment = 0.0;
const double nu_N = 0.0;
const double Deut_Mag_Moment = 0.857;//,in unit of magneton (nu_N), for s=1,l=0 state; For s=1, l=1, Mag_moment=0.310;
const double Deut_Quad_Moment = 2.286 / HBARC / HBARC ; //fm2 into MeV^-2

const double SINTHW = 0.230;          //sin^2(theta_W)
const double GV = 0.080;          //1-4.0*sinthw
const double K = 0.020;         //GV/4
const double FM2NB = 1e+7; //from ub to fm^2, 1fm^2 = 10mb, 1b=1e-28 m^2, 1fm = 1e-15m
const double NB2FM = 1e-7; //from ub to fm^2, 1fm^2 = 10mb, 1b=1e-28 m^2, 1fm = 1e-15m
const double FM2MB = 1e+1; //from ub to fm^2, 1fm^2 = 10mb, 1b=1e-28 m^2, 1fm = 1e-15m
const double MB2NB = 1e+6; 
const double UB2NB = 1e+3; 
const double CM2NB = 1e+33; 
const double NB2CM = 1e-33; 
const double FM2CM = 1e-26; //from fm^2 to cm^2, 1fm^2 = 10mb, 1b=1e-28 m^2, 1fm = 1e-15m
