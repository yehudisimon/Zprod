// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <cmath>      	// Mathematical functions       //
#include <iostream>
#include <cstddef>
#include <math.h>
//#include <numbers>
#include <complex>
//#include <boost/math/special_functions/digamma.hpp>
// -------- Classes ----------------------------------- //
#include "process.h"	// Process//
//#include "Constants.h"  //Constants//
// -------- Functions --------------------------------- //
extern "C" {
  double ct18pdf_(int&, double&, double&);              //
  double ctq6pdf_(int&, double&, double&);              //
};

std::complex<double> Psi(std::complex<double>);
//void SetPDF(const int&, double&, double&, double*, double*, double&);
//void SetCouplings(const int&, Process*, double&, double&, double&, double&);
void SetCouplings(const int&, Process*, double&, double&, double&);
void Kinematics(double&, double&, double&, double&, double*, Process*);
void Kinematics2(double&, double&, double&, double&, double*, Process*);
void Kinematics2Plus(double&, double&, double&, double&, double*, Process*);
void SetPDFN(std::complex<double>&,std::complex<double>*,std::complex<double>*,std::complex<double>&,double (*)[8]);
void EvolvePDF(std::complex<double>&, std::complex<double>*, std::complex<double>*, std::complex<double>&);
//void EvolvePDFExp(std::complex<double>&, std::complex<double>*, std::complex<double>*, std::complex<double>&);


// ---------------------------------------------------- //
extern double MZ, alpha, tw, gu, gd, alpha_s, muF, muR;
extern double A [8][8];
extern int Nflav;
std::complex<double> one=1.;
//const double M_EULER = std::numbers::egamma_v;
// ************************************************************************* //
//  TOT: Born integrand .                        //
// ************************************************************************* //
double Born(double *x, size_t dim, void *prm)
{
   // conversion of type void* into Param*
   Process *p = (Process *)prm;

   // Computation of the kinematical variables
   double xa, xb, s, jac;
   Kinematics(xa,xb,jac,s,x,p);

   // PDFs
   double fAB=0;
   
   SetCouplings(0,p,fAB,xa,xb); //Sets PDFs with their appropriate u/d coupling, summed and symmetrized
   
    // result
   return pow(M_PI,2.)*2./3.*alpha*jac*fAB/xa/xb/p->sh*0.38937966e9; //Integration of sigma_0=2*Pi^2*alpha/(3*S), the Born cross section with no partonic (xa or xb) dependance, conversion factor from GeV^-2 to pb
}


// *********************************************************************** //
//  Real Correction Integrand //
// *********************************************************************** //

double Real(double *x, size_t dim, void *prm)
{
   // conversion of type void* into Param*
   Process *p = (Process *)prm;

   // Computation of the kinematical variables
   double xa, xb, s, jac;
   Kinematics2(xa,xb,jac,s,x,p);

   // Couplings and scale setting
   double fAB=0;
   double res=0;
   double z=pow(MZ,2.)/p->sh/xa/xb; //z variable to express functions more easily
   double prefqg=0; //prefactors for real (qg->zq) NLO subprocesses
   double integqg=0; //Function to integrate for real NLO subprocesses needing 2 variables

   prefqg=alpha*alpha_s*M_PI/3./p->sh;

   integqg=(1+7*z)*(1-z)/4.+(pow(z,2.)+pow(1-z,2.))/2.*log((1-z)*MZ/muF/sqrt(z))*2;
   

   SetCouplings(3,p,fAB,xa,xb); //Sets PDFs with their appropriate u/d coupling, summed and symmetrized

   res=prefqg*integqg*fAB*jac/xa/xb*0.38937966e9;
   // result
   return res;
}

// *************************************************** //
// Virtual Correction Integrand //
// *************************************************** //


double Virtual1D(double *x, size_t dim, void *prm) //1D "delta(1-z)" contribution of virtual correction
{
   // conversion of type void* into Param*
   Process *p = (Process *)prm;

   // Computation of the kinematical variables
   double xa, xb, s, jac;
   Kinematics(xa,xb,jac,s,x,p);

   // Couplings and scale setting
   double fAB=0;
   double res=0;
   double prefqq=0; //Prefactors for virtual (qqbar->z) NLO subprocesses
   double integqq=0; //Function to integrate for virtual NLO subprocesses needing 2 variables

   prefqq=alpha*alpha_s*M_PI*8./9./p->sh;
   integqq=3.*log(MZ/muF)+pow(M_PI,2.)/3.-4;

   SetCouplings(0,p,fAB,xa,xb); //Sets PDFs with their appropriate u/d coupling, summed and symmetrized

   res=prefqq*integqq*fAB*jac/xa/xb*0.38937966e9;
   // result
   return res;
}


double Virtual2D(double *x, size_t dim, void *prm) //2D contribution of virtual correction
{
   // conversion of type void* into Param*
   Process *p = (Process *)prm;

   // Computation of the kinematical variables
   double xa, xb, s, jac;
   Kinematics2(xa,xb,jac,s,x,p);

   // Couplings and scale setting
   double fAB=0;
   double res=0;
   double z=pow(MZ,2.)/p->sh/xa/xb; //z variable to express functions more easily
   double prefqq=0; //Prefactors for virtual (qqbar->z) NLO subprocesses
   double integqq=0; //Function to integrate for virtual NLO subprocesses needing 2 variables

   prefqq=alpha*alpha_s*M_PI*8./9./p->sh;
   integqq=-(1+z*z)*log(z)/(1-z);

   SetCouplings(0,p,fAB,xa,xb); //Sets PDFs with their appropriate u/d coupling, summed and symmetrized

   res=prefqq*integqq*fAB*jac/xa/xb*0.38937966e9;
   // result
   return res;
}

double VirtualPlus(double *x, size_t dim, void *prm) //Plus distribution contribution
{
   // conversion of type void* into Param*
   Process *p = (Process *)prm;

   // Computation of the kinematical variables
   double s, jac1, y, z,jac2;
   Kinematics2Plus(y,z,jac1,s,x,p);
   jac2=-log(pow(MZ,2.)/p->sh); // jacobian for 1D +distribution leftovers

   // Couplings and scale setting
   double fAB1=0, fAB2=0;
   double res=0;
   double prefqq=alpha*alpha_s*M_PI*8./9./p->sh; // Prefactor for q qbar > z NLO subprocess

   SetCouplings(1,p,fAB1,y,z); //Sets PDFs with their appropriate u/d coupling, summed and symmetrized
   SetCouplings(2,p,fAB2,y,z); //Sets PDFs with their appropriate u/d coupling, summed and symmetrized

   res=prefqq*0.38937966e9*(jac1*fAB1+jac2*fAB2);
   // result
   return res;
}


//------------------------------------------------------------//
//----------------Resummed computation------------------------//
//------------------------------------------------------------//

double Resummed(double *x, size_t dim, void *prm) // Resummed formula, the numericall integration does the inversed Mellin transform
{

  // Resummation and kinematics variables
  // !! Nb=N*exp(gamma) !!
  double H=0, beta0=23./6., beta1=29./3., C=2.015, Phi=1.7*M_PI/2.;
  double res=0;
  std::complex<double> i(0.0,1.0), N=0, Nm=0, Nb=0, fAB=0, G=0, g1=0, g2=0, jac=0;


   // conversion of type void* into Param*
   Process *p = (Process *)prm;

   //Tau computation
   double tau = pow(MZ,2)/p->sh;

   // Change of variable
   N= C+ (cos(Phi)+i*sin(Phi))*tan(M_PI*x[0]/2.);
   jac = 1./2.*(cos(Phi)+i*sin(Phi))*(1.+pow(tan(M_PI*x[0]/2.),2));
   Nb=N*exp(-Psi(1.));
   //Nm=(N-1.)*exp(-Psi(1.));
   //Np=N+1.;
   
   //Initializing the couplings and PDFs
   std::complex<double> q[5], g, qbar[5];
   SetPDFN(N,q,qbar,g,A); // PDFs with N
   EvolvePDF(N,q,qbar,g);

   // for(int k0=0; k0<5; k0++){
   //   std::cout << "k0 = " << k0 << " qE = " << q[k0] << " qbarE = " << qbar[k0] << std::endl;
   // }
   // std::cout << " g = " << g << std::endl;

   
   for(int i=0; i<Nflav;i++){
     if(i%2==0) fAB+=gd*q[i]*qbar[i];
     else fAB+=gu*q[i]*qbar[i];
   }
   
   
   // fAB+=gd*q[0]*qbar[0]; //d quark
   // fAB+=gd*q[2]*qbar[2]; //s quark
   // fAB+=gd*q[4]*qbar[4]; //b quark
   // fAB+=gu*q[1]*qbar[1]; //u quark
   // fAB+=gu*q[3]*qbar[3]; //c quark

   
   //Global factors
   fAB*=2.*pow(tau,-N);

   if(pow(tau,-N) != pow(tau,-N)) std::cout << "N = " << N << "tau⁻N = " << pow(tau,-N) << std::endl;
   
   //Global H factor, independant of which correction we consider
   H=pow(M_PI,2.)*2./3.*alpha/p->sh*0.38937966e9*(1.+alpha_s/M_PI*4./3.*(2*pow(M_PI,2.)/3.-4.)); //H = sigma_B*(1+alpha_s/(2*Pi)(A_0+Pi^2*A^1_q/6))

   // g1, g2 and corresponding G functions of the N variable, independant of the correction considered (global)
   //g1=16./3./beta0*(1.+log(1.-alpha_s*beta0*log(Nb)/M_PI)*M_PI/(alpha_s*beta0*log(Nb)));
   g1=16./3./beta0*(1.+log(1.-alpha_s*beta0*log(Nb)/M_PI)*M_PI/(alpha_s*beta0*log(Nb)));
   
   g2=-4./beta0*log(1.-alpha_s*beta0/M_PI*log(Nb))+beta1*pow(beta0,-3)*4./3.*(pow(log(1.-alpha_s*beta0*log(Nb)/M_PI),2)+2.*alpha_s*beta0*log(Nb)/M_PI+2.*log(1.-alpha_s*beta0*log(Nb)/M_PI))+(alpha_s/M_PI*beta0*log(Nb)+log(1.-alpha_s/M_PI*beta0*log(Nb)))*(16./3./beta0*log(MZ/muR)-8./3.*pow(beta0,-2)*(67./6.-pow(M_PI,2)/2.-25./9.));

   //g2=-4./beta0*log(1.-alpha_s*beta0/M_PI*log(Nm))+beta1*pow(beta0,-3)*4./3.*(pow(log(1.-alpha_s*beta0*log(Nm)/M_PI),2)+2.*alpha_s*beta0*log(Nm)/M_PI+2.*log(1.-alpha_s*beta0*log(Nm)/M_PI))+(alpha_s/M_PI*beta0*log(Nm)+log(1.-alpha_s/M_PI*beta0*log(Nm)))*(16./3./beta0*log(MZ/muR)-8./3.*pow(beta0,-2)*(67./6.-pow(M_PI,2)/2.-25./9.));
   
   G=std::exp(g1*log(Nb)+g2);
   //G=std::exp(g1*log(Nm)+g2);

   
   //std::cout << " fAB = " << fAB << " G = " << G << " H = " << H << " tau^-N = " << pow(tau,-N) << " N= " << N << " fqq = " << fqq << " fgq = " << fgq << "Eq = " << Eq << " Eg = " << Eg << " s = " << p->sh << " res = " << std::imag(G*fAB*jac)*H << std::endl;
    
   // Result
   res=std::imag(G*fAB*jac)*H;
   //if (res != 0) std::cout << "res=" << G*jac*fAB*H << std::endl;
   
   return res;

}
  

double ResummedHS(double *x, size_t dim, void *prm) // Resummed formula, the numericall integration does the inversed Mellin transform
{

  // Resummation and kinematics variables
  // !! Nb=N*exp(gamma) !!
  double H=0, beta0=23./6., beta1=29./3., C=2.015, Phi=1.7*M_PI/2.;
  double res=0;
  std::complex<double> i(0.0,1.0), N=0, Np=0, Nb=0, fAB=0, G=0, g1=0, g2=0, jac=0;


   // conversion of type void* into Param*
   Process *p = (Process *)prm;

   //Tau computation
   double tau = pow(MZ,2)/p->sh;

   // Change of variable
   N= C+ (cos(Phi)+i*sin(Phi))*tan(M_PI*x[0]/2.);
   jac = 1./2.*(cos(Phi)+i*sin(Phi))*(1.+pow(tan(M_PI*x[0]/2.),2));
   Nb=N*exp(-Psi(1.));
   //Nm=(N-1.)*exp(-Psi(1.));
   Np=N+1.;
   
   //Initializing the couplings and PDFs
   std::complex<double> q[5], g, qbar[5];
   SetPDFN(Np,q,qbar,g,A); // PDFs with N
   EvolvePDF(Np,q,qbar,g);

   for(int i=0; i<Nflav;i++){
     if(i%2==0) fAB+=gd*q[i]*qbar[i];
     else fAB+=gu*q[i]*qbar[i];
   }
     
   //Global factors
   fAB*=2.*pow(tau,-N);

   if(pow(tau,-N) != pow(tau,-N)) std::cout << "N = " << N << "tau⁻N = " << pow(tau,-N) << std::endl;
   
   //Global H factor, independant of which correction we consider
   H=pow(M_PI,2.)*2./3.*alpha/MZ/MZ*0.38937966e9*(1.+alpha_s/M_PI*4./3.*(2*pow(M_PI,2.)/3.-4.)); //H = sigma_B*(1+alpha_s/(2*Pi)(A_0+Pi^2*A^1_q/6))

   g1=16./3./beta0*(1.+log(1.-alpha_s*beta0*log(Nb)/M_PI)*M_PI/(alpha_s*beta0*log(Nb)));
     
   g2=-4./beta0*log(1.-alpha_s*beta0/M_PI*log(Nb))+beta1*pow(beta0,-3)*4./3.*(pow(log(1.-alpha_s*beta0*log(Nb)/M_PI),2)+2.*alpha_s*beta0*log(Nb)/M_PI+2.*log(1.-alpha_s*beta0*log(Nb)/M_PI))+(alpha_s/M_PI*beta0*log(Nb)+log(1.-alpha_s/M_PI*beta0*log(Nb)))*(16./3./beta0*log(MZ/muR)-8./3.*pow(beta0,-2)*(67./6.-pow(M_PI,2)/2.-25./9.));
   
   //G=std::exp(g1*log(Nb)+g2);
   G=std::exp(g1*log(Nb)+g2);

   // Result
   res=std::imag(G*fAB*jac)*H;
   return res;

}


double ResummedHSConvert(double *x, size_t dim, void *prm) // Resummed formula, the numericall integration does the inversed Mellin transform
{

  // Resummation and kinematics variables
  // !! Nb=N*exp(gamma) !!
  double H=0, beta0=23./6., beta1=29./3., C=1.015, Phi=1.7*M_PI/2.;
  double res=0;
  std::complex<double> i(0.0,1.0), N=0, Npb=0, Np=0, Nb=0, fAB=0, G=0, g1=0, g2=0, jac=0;
  
  // conversion of type void* into Param*
  Process *p = (Process *)prm;
  
  //Tau computation
  double tau = pow(MZ,2)/p->sh;
  
  // Change of variable
  N= C+ (cos(Phi)+i*sin(Phi))*tan(M_PI*x[0]/2.);
  jac = 1./2.*(cos(Phi)+i*sin(Phi))*(1.+pow(tan(M_PI*x[0]/2.),2));
  Nb=N*exp(-Psi(1.));
  Npb=(N+1.)*exp(-Psi(1.));
  Np=N+1.;
  
  //Initializing the couplings and PDFs
  std::complex<double> q[5], g, qbar[5];
  SetPDFN(Np,q,qbar,g,A); // PDFs with N
  EvolvePDF(Np,q,qbar,g);
  
  for(int i=0; i<Nflav;i++){
    if(i%2==0) fAB+=gd*q[i]*qbar[i];
    else fAB+=gu*q[i]*qbar[i];
  }
  
  //Global factors
  fAB*=2.*pow(tau,-N);

  if(pow(tau,-N) != pow(tau,-N)) std::cout << "N = " << N << "tau⁻N = " << pow(tau,-N) << std::endl;
  
  //Global H factor, independant of which correction we consider
  H=pow(M_PI,2.)*2./3.*alpha/MZ/MZ*0.38937966e9*(1.+alpha_s/M_PI*4./3.*(2*pow(M_PI,2.)/3.-4.)); //H = sigma_B*(1+alpha_s/(2*Pi)(A_0+Pi^2*A^1_q/6))

  g1=16./3./beta0*(1.+log(1.-alpha_s*beta0*log(Npb)/M_PI)*M_PI/(alpha_s*beta0*log(Npb)));
     
  g2=-4./beta0*log(1.-alpha_s*beta0/M_PI*log(Npb))+beta1*pow(beta0,-3)*4./3.*(pow(log(1.-alpha_s*beta0*log(Npb)/M_PI),2)+2.*alpha_s*beta0*log(Npb)/M_PI+2.*log(1.-alpha_s*beta0*log(Npb)/M_PI))+(alpha_s/M_PI*beta0*log(Npb)+log(1.-alpha_s/M_PI*beta0*log(Npb)))*(16./3./beta0*log(MZ/muR)-8./3.*pow(beta0,-2)*(67./6.-pow(M_PI,2)/2.-25./9.));
   
  //G=std::exp(g1*log(Nb)+g2);
  G=std::exp(g1*log(Npb)+g2);
  
  // Result
  res=std::imag(G*fAB*jac)*H;
  return res;
}


double Expanded(double *x, size_t dim, void *prm) // Expansion from resummed formula, the numericall integration does the inversed Mellin transform
{
  
  // Resummation and kinematics variables
  // !! Nb=N*exp(gamma_E) !!
  double C=2.015; // 2 < C < C_Landau =exp(Pi/2*beta0*alpha_s) ~ 32.2 
  double res=0, Phi=1.7*M_PI/2.;
  std::complex<double> i(0.0,1.0), N=0, Nm=0, Nb=0, fAB=0, fqq=0, fqg=0, jac=0;

  // conversion of type void* into Param*
  Process *p = (Process *)prm;

  // Tau computation
  double tau = pow(MZ,2)/p->sh;   	

  // Change of variable
  N= C+ (cos(Phi)+i*sin(Phi))*tan(M_PI*x[0]/2.);
  jac = 1./2.*(cos(Phi)+i*sin(Phi))*(1.+pow(tan(M_PI*x[0]/2.),2));
  //Np=N+1.;
  Nb=N*exp(-Psi(1.));
  //Nm=(N-1.)*exp(-Psi(1.));
  
  //Initializing the couplings and PDFs
  std::complex<double> q[5], qbar[5], g, Gexp=0, Hexp=0, Counter=0;
  SetPDFN(N,q,qbar,g,A); // PDFs with N
  

  for(int i=0; i<Nflav;i++){
    if(i%2==0) fqq+=gd*q[i]*qbar[i];
    else fqq+=gu*q[i]*qbar[i];
  }

  // Gexp+=gd*q[0]*qbar[0]; //d quark
  // Gexp+=gd*q[2]*qbar[2]; //s quark
  // Gexp+=gd*q[4]*qbar[4]; //b quark
  // Gexp+=gu*q[1]*qbar[1]; //u quark
  // Gexp+=gu*q[3]*qbar[3]; //c quark

  

  // Counter=Gexp*pow(4./3.*(3./2.+1./N/(N+1.)-2.*(Psi(N+1.)-Psi(1.))),2);
  // Counter+=pow(g*(2.+N+N*N)/(2.*N*(N+1.)*(N+2.)),2)*(3.*gd+2.*gu);

  // Counter*=pow(alpha_s/M_PI*log(MZ/Nb/muF),2);
  
  // for(int i=0; i<5;i++){
  //   if(i%2==0) Counter+=g*(q[i]+qbar[i])*gd*(2.+N+N*N)/(2.*N*(N+1.)*(N+2.))*4./3.*(2.+N+N*N)/(N*(N*N-1.));
  //   else Counter+=g*(q[i]+qbar[i])*gu*(2.+N+N*N)/(2.*N*(N+1.)*(N+2.))*4./3.*(2.+N+N*N)/(N*(N*N-1.));
  // }
  
  Gexp=alpha_s*4./3./M_PI*(3.*log(Nb)-2.*pow(log(Nb),2));
  //Gexp=alpha_s*4./3./M_PI*(3.*log(Nm)-2.*pow(log(Nm),2));
  Hexp=alpha_s*4./3./M_PI*(2.*M_PI*M_PI/3.-4.);

  
    
  //EvolvePDFExp(N,q,qbar,g);
  
  // // qqbar contribution symmetrized
  // fqq+=gd*q[0]*qbar[0]; //d quark
  // fqq+=gd*q[2]*qbar[2]; //s quark
  // fqq+=gd*q[4]*qbar[4]; //b quark
  // fqq+=gu*q[1]*qbar[1]; //u quark
  // fqq+=gu*q[3]*qbar[3]; //c quark
   
  // qqbar Cross section normalized by Born
  fqq*=(1.+alpha_s/M_PI*(3.-4.*(Psi(N+1.)-Psi(1.))+2./N/(N+1.))*4./3.*log(MZ/Nb/muF)+Gexp+Hexp);
  //fqq=fqq*(1.+alpha_s/M_PI*(8./3.*pow(log(Nb),2)-16./3.*log(Nb)*log(MZ/muF)+8./3.*(log(Nb)-log(MZ/muF))/N+8.*pow(M_PI,2)/9.-16./3.+4*log(MZ/muF)));
  //fqq*=(1.+alpha_s/M_PI*(-8./3.*pow(log(Nb),2)+16./3.*log(Nb)*(Psi(N+1.)-Psi(1.))-8./3.*log(Nb)/N/(N+1.)+8.*pow(M_PI,2)/9.-16./3.));
  //fqq*=(1.+alpha_s/M_PI*(+8./3.*pow(log(Nb),2)+8./3.*log(Nb)/N+8.*pow(M_PI,2)/9.-16./3.));
  //std::cout << "fqq = " << fqq << " f_u = " << q[1] << " f_bar = " << qbar[1] << " f_c = " << q[3] << "f_cbar = " << qbar[3] << " f_b = " << q[4] << "f_bbar = " << qbar[4] << " f_s = " << q[2] << " f_sbar= " << qbar[2] << " f_d = " << q[0] << " f_dbar = " << qbar[0] << " prefactor = " << 1.+alpha_s/M_PI*(8./3.*pow(log(Nb),2)-16./3.*log(Nb)*log(MZ/muF)+8./3.*(log(Nb)-log(MZ/muF))/N+8.*pow(M_PI,2)/9.-16./3.+4*log(MZ/muF))  <<  std::endl;  

  
  for(int i=0; i<Nflav; i++){
    if(i%2==0) fqg+=gd*g*(q[i]+qbar[i]);
    else fqg+=g*(q[i]+qbar[i])*gu;
  }

  // // qg and qbarg contributions symetrized
  // fqg+=gd*g*(q[0]+qbar[0]); //d and dbar
  // fqg+=gd*g*(q[2]+qbar[2]); //s and dbar
  // fqg+=gd*g*(q[4]+qbar[4]); //b and dbar
  // fqg+=g*gu*(q[1]+qbar[1]); //u and dbar
  // fqg+=g*gu*(q[3]+qbar[3]); //c and dbar
   
  // qg cross section normalized by Born
  fqg*=alpha_s/2./M_PI/N*log(MZ/muF/Nb)*(2.+N+N*N)/(N+1.)/(N+2.);
  
  //Global factors
  //fAB=2.*(fqq+Gexp+Hexp-Counter)*pow(tau,-N); //*2.;?
  fAB=2.*(fqq+fqg)*pow(tau,-N); //*2.;?
  
  
  // Result
  res=std::imag(fAB*jac)*pow(M_PI,2.)*2./3.*alpha/p->sh*0.38937966e9;
  
  return res;
  
}

double ExpandedHS(double *x, size_t dim, void *prm) // Expansion from resummed formula, the numericall integration does the inversed Mellin transform
{
  
  // Resummation and kinematics variables
  // !! Nb=N*exp(gamma_E) !!
  double C=2.015; // 2 < C < C_Landau =exp(Pi/2*beta0*alpha_s) ~ 32.2 
  double res=0, Phi=1.7*M_PI/2.;
  std::complex<double> i(0.0,1.0), N=0, Np=0, Npb=0, Nb=0, fAB=0, fqq=0, fqg=0, jac=0;

  // conversion of type void* into Param*
  Process *p = (Process *)prm;

  // Tau computation
  double tau = pow(MZ,2)/p->sh;   	

  // Change of variable
  N= C+ (cos(Phi)+i*sin(Phi))*tan(M_PI*x[0]/2.);
  jac = 1./2.*(cos(Phi)+i*sin(Phi))*(1.+pow(tan(M_PI*x[0]/2.),2));
  Np=N+1.;
  Nb=N*exp(-Psi(1.));
  Npb=(N+1.)*exp(-Psi(1.));
  
  //Initializing the couplings and PDFs
  std::complex<double> q[5], qbar[5], g, Gexp=0, Hexp=0, Counter=0;
  SetPDFN(Np,q,qbar,g,A); // PDFs with N
  

  for(int i=0; i<Nflav;i++){
    if(i%2==0) fqq+=gd*q[i]*qbar[i];
    else fqq+=gu*q[i]*qbar[i];
  }

  Gexp=alpha_s*4./3./M_PI*(3.*log(Nb)-2.*pow(log(Nb),2));
  Hexp=alpha_s*4./3./M_PI*(2.*M_PI*M_PI/3.-4.);
   
  // qqbar Cross section normalized by Born
  fqq*=(1.+alpha_s/M_PI*(3.-4.*(Psi(Np+1.)-Psi(1.))+2./Np/(Np+1.))*4./3.*log(MZ/Npb/muF)+Gexp+Hexp);

  for(int i=0; i<Nflav; i++){
    if(i%2==0) fqg+=gd*g*(q[i]+qbar[i]);
    else fqg+=g*(q[i]+qbar[i])*gu;
  }
   
  // qg cross section normalized by Born
  fqg*=alpha_s/2./M_PI/Np*log(MZ/muF/Npb)*(2.+Np+Np*Np)/(Np+1.)/(Np+2.);
  
  //Global factors
  //fAB=2.*(fqq+Gexp+Hexp-Counter)*pow(tau,-N); //*2.;?
  fAB=2.*(fqq+fqg)*pow(tau,-N); //*2.;?
  
  
  // Result
  res=std::imag(fAB*jac)*pow(M_PI,2.)*2./3.*alpha/MZ/MZ*0.38937966e9;
  
  return res;
  
}


double ExpandedHSConvert(double *x, size_t dim, void *prm) // Expansion from resummed formula, the numericall integration does the inversed Mellin transform
{
  
  // Resummation and kinematics variables
  // !! Nb=N*exp(gamma_E) !!
  double C=1.015; // 2 < C < C_Landau =exp(Pi/2*beta0*alpha_s) ~ 32.2 
  double res=0, Phi=1.7*M_PI/2.;
  std::complex<double> i(0.0,1.0), N=0, Np=0, Npb=0, Nb=0, fAB=0, fqq=0, fqg=0, jac=0;

  // conversion of type void* into Param*
  Process *p = (Process *)prm;

  // Tau computation
  double tau = pow(MZ,2)/p->sh;   	

  // Change of variable
  N= C+ (cos(Phi)+i*sin(Phi))*tan(M_PI*x[0]/2.);
  jac = 1./2.*(cos(Phi)+i*sin(Phi))*(1.+pow(tan(M_PI*x[0]/2.),2));
  Np=N+1.;
  Nb=N*exp(-Psi(1.));
  Npb=(N+1.)*exp(-Psi(1.));
  
  //Initializing the couplings and PDFs
  std::complex<double> q[5], qbar[5], g, Gexp=0, Hexp=0, Counter=0;
  SetPDFN(Np,q,qbar,g,A); // PDFs with N
  

  for(int i=0; i<Nflav;i++){
    if(i%2==0) fqq+=gd*q[i]*qbar[i];
    else fqq+=gu*q[i]*qbar[i];
  }

  Gexp=alpha_s*4./3./M_PI*(3.*log(Npb)-2.*pow(log(Npb),2));
  Hexp=alpha_s*4./3./M_PI*(2.*M_PI*M_PI/3.-4.);
   
  // qqbar Cross section normalized by Born
  fqq*=(1.+alpha_s/M_PI*(3.-4.*(Psi(Np+1.)-Psi(1.))+2./Np/(Np+1.))*4./3.*log(MZ/Npb/muF)+Gexp+Hexp);

  for(int i=0; i<Nflav; i++){
    if(i%2==0) fqg+=gd*g*(q[i]+qbar[i]);
    else fqg+=g*(q[i]+qbar[i])*gu;
  }
   
  // qg cross section normalized by Born
  fqg*=alpha_s/2./M_PI/Np*log(MZ/muF/Npb)*(2.+Np+Np*Np)/(Np+1.)/(Np+2.);
  
  //Global factors
  //fAB=2.*(fqq+Gexp+Hexp-Counter)*pow(tau,-N); //*2.;?
  fAB=2.*(fqq+fqg)*pow(tau,-N); //*2.;?
  
  
  // Result
  res=std::imag(fAB*jac)*pow(M_PI,2.)*2./3.*alpha/MZ/MZ*0.38937966e9;
  
  return res;
  
}


//------------------------------------------------------------//
//--------------------Total Integrand-------------------------//
//------------------------------------------------------------//


double Tot(double *x, size_t dim, void *prm)
{
  double res=0;
   	
	Process *p = (Process *)prm;
	
	if(p->NLO==0) res = Born(x,dim,prm); 
	else if(p->NLO==1) res= Born(x,dim,prm) + Virtual1D(x,dim,prm) + Virtual2D(x,dim,prm) + VirtualPlus(x,dim,prm) +  Real(x,dim,prm);
	else if(p->NLO==2) res = Resummed(x,dim,prm);
	else if(p->NLO==3) res = Expanded(x,dim,prm);
	else if(p->NLO==4) res = ResummedHS(x,dim,prm);
	else if(p->NLO==5) res = ExpandedHS(x,dim,prm);
	else if(p->NLO==6) res = ResummedHSConvert(x,dim,prm);
	else if(p->NLO==7) res = ExpandedHSConvert(x,dim,prm);
	else std::cout << "Process ID ERROR ID = " << p->NLO << std::endl;
	
	return res;
}
