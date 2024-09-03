#ifndef CONSTANTS_H
#define CONSTANTS_H
// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <cmath>
#include <math.h>
#include <string>    // String streams                  //
// -------- Processor variables ----------------------- //
#define cst const std::string                           //
// ------------------EW Constants ----------------------//
//Constants found in the MadGraph configuration
double MZ = 91.188; //Z boson mass in GeV
double GF = 1.16639e-5; //Fermi constant
//double alpha_b = 7.29735257e-3;  //alpha
double alpha_s = 0.118; //alpha_s(MZ)
//double alpha_s = 0.05; //alpha_s(MZ)
double alphaM1 = 1.32507e+2; //1/alpha
double alpha = 1/alphaM1;//7.5467711139788835e-3;  //alpha(MZ)
//double tw = (1.+sqrt(1.-4.*M_PI*alpha/(sqrt(2.)*pow(MZ,2.)*GF)))/2.; //tw=cos(theta_w)^2=mW^2/mZ^2
double mw2 = MZ*MZ/2. + sqrt(pow(MZ,4)/4. - alpha*M_PI*MZ*MZ/(GF*sqrt(2)));
double tw = mw2/(MZ*MZ);
double tw2 = pow(0.88190334743339216,2.); // cos(theta_w)^2=MW^2/MZ^2 at scale MZ

double gu=(8./9.*pow(1-tw,2.)+1./4.-2./3.*(1-tw))/(tw*(1-tw)); //g_L^2+g_R^2 for u type quark
double gd=(pow(1-tw,2.)*2./9.+1./4.-(1-tw)/3.)/(tw*(1-tw)); //g_L^2+g_R^2 for d type quark

double muF=MZ;
double muR=MZ;

int Nflav=5;
#endif
