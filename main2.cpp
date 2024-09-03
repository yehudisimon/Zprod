// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <iostream>   // In/Out streams                 //
#include <cmath>
#include <fstream>
// -------- Classes ----------------------------------- //
#include "process.h"      // Process definition         //
#include "messages.h"     // Message services           //
#include "PDFfitAN.h"     // Fit coefficients

#include "LHAPDF/LHAPDF.h"

using namespace std;
using namespace LHAPDF;


//#include "Constants.h"



extern "C" {
  double ct18pdf_(int&, double&, double&);              //                                                           
  double ctq6pdf_(int&, double&, double&);              //                                                            
};

// -------- Functions --------------------------------- //
void Integrate(double&,double&,double&,Process*);  // Integration
void pdfFit(double&, double (*)[8], double&, double, double, double, Process*); //


// ---------------------------------------------------- //
//
extern double MZ, muF, muR, A1min;
extern double A[8][8];
// ************************************************************************* //
//  Main code                                                                //
// ************************************************************************* //
int main(int argc, char* argv[])
{
  // Checking the arguments
  //if(argc!=5)
  if(argc!=2)
  {
    //std::cout << "This function requires 1 argument\n"
    std::cout << "This function requires 2 argument\n"
              << "  - PDF ID (1=CT18 NLO, 2=CT14LN, 3=CT14LL, 4=CTEQ6M)\n";
    exit(0);
  }	
  

  for (const string& p : paths())
    cout << p << endl;
  cout << "@" << findFile("lhapdf.conf") << "@" << endl;
  cout << "List of available PDFs:" << endl;
  for (const string& s : availablePDFSets())
    cout << " " << s << endl;
      
 
  // std::cout << "making PDF" << std::endl;
  // const PDF* pdf = mkPDF("CT18NLO", 0);
    
  //Opening writing file
  std::ofstream LO;
  std::ofstream NLO;
  std::ofstream Expand;
  std::ofstream Resummed;
  
  // Initializing message services
  InitMessages();

  // Setup of the collision setup, pdfs, scales, etc...
  const int nsteps=100, j0=0;//, pdf=int(argv[1]);
  double s0=300., sf=1e6, sc=s0;
  //const char *J;
  //J = std::to_string(j0).c_str();
  
  
  Process *Proc = new Process(argv[1], sc, j0);

  Proc->Init(); //PDF initialisation
  
  //Fitting PDFs to get Ai coefficients to include in the F_j(N) Mellin transform of PDFs
  //if(Proc->NLO==2 || Proc->NLO==3) pdfFit(A1min,A,muF,-1.6,-1.6,-1.6,Proc); 
 
  
  //Fitting gluon PDF...#chisq/dof = 5.08633e-05
  int i0=0;
  A[i0][0] =  3.07320 ;
  A[i0][1] =  -1.40521 ;
  A[i0][2] =  13.19817 ;
  A[i0][3] =  -11.39889;
  A[i0][4] =  85.40350 ;
  A[i0][5] =  -324.82879;
  A[i0][6] =  595.24922 ;
  A[i0][7] =  -374.83165 ;
 
  //Fitting valence down quark PDF...#chisq/dof = 1.03978e-07
  i0=1;
  A[i0][0] =  0.10616;
  A[i0][1] =  -1.37718;
  A[i0][2] =  5.59590 ;
  A[i0][3] =  0.66416 ;
  A[i0][4] =  11.26630 ;
  A[i0][5] =  -8.78590 ;
  A[i0][6] =  48.02226 ;
  A[i0][7] =  -41.87545 ;


  //Fitting valence up quark PDF...#chisq/dof = 1.29063e-07
  i0=2;
  A[i0][0] =  0.14032;
  A[i0][1] =  -1.35731;
  A[i0][2] =  4.01971 ;
  A[i0][3] =  -2.72381 ;
  A[i0][4] =  44.79915 ;
  A[i0][5] =  -112.29431;
  A[i0][6] =  198.33459 ;
  A[i0][7] =  -115.50177 ;
  
  //Fitting sea down quark PDF...#chisq/dof = 7.17778e-08
  i0=3;
  A[i0][0] =  0.11747;
  A[i0][1] =  -1.36943;
  A[i0][2] =  3.04292 ;
  A[i0][3] =  -2.30502;
  A[i0][4] =  12.30102 ;
  A[i0][5] =  -44.87680;
  A[i0][6] =  61.78924 ;
  A[i0][7] =  -28.50200;

  //strange quark
  i0=4;
  A[i0][0] =  0.14421;
  A[i0][1] =  -1.35119;
  A[i0][2] =  10.34872 ;
  A[i0][3] =  -6.60674 ;
  A[i0][4] =  36.68691 ;
  A[i0][5] =  -114.60239;
  A[i0][6] =  184.16840 ;
  A[i0][7] =  -105.63378;

  //Fitting bottom quark PDF...#chisq/dof = 1.7599e-08
  i0=5;
  A[i0][0] =  0.09036;
  A[i0][1] =  -1.37489;
  A[i0][2] =  14.52458 ;
  A[i0][3] =  -11.21468 ;
  A[i0][4] =  81.30858 ;
  A[i0][5] =  -305.25437;
  A[i0][6] =  562.65945 ;
  A[i0][7] =  -361.30952;


  //Fitting sea up quark PDF...#chisq/dof = 8.35594e-09
  i0=6;
  A[i0][0] =  0.17969;
  A[i0][1] =  -1.33666;
  A[i0][2] =  10.04120 ;
  A[i0][3] =  -7.24955 ;
  A[i0][4] =  53.29604 ;
  A[i0][5] =  -185.10252;
  A[i0][6] =  303.01405 ;
  A[i0][7] =  -166.09590 ;

  //Fitting charm quark PDF...#chisq/dof = 1.61266e-08
  i0=7;
  A[i0][0] =  0.12213;
  A[i0][1]=  -1.36013 ;
  A[i0][2] =  10.60837 ;
  A[i0][3] =  -8.77261 ;
  A[i0][4] =  51.19531 ;
  A[i0][5] =  -160.11007;
  A[i0][6] =  257.14663 ;
  A[i0][7] =  -148.30440 ;
 

  LO.open("LOXsec.txt");
  NLO.open("NLOXsec.txt");
  Resummed.open("ResumXsec.txt");
  Expand.open("ExpandedXsec.txt");
  

  for(int ic=0; ic <= nsteps; ic++) 
    {

      for(int kID=0; kID < 4; kID++)
	{
	  //J = std::to_string(kID).c_str();
	  Process *Proc = new Process(argv[1], sc, kID);
	  
	  // LO/NLO/Resummed/Expanded computation
	  std::cout << "k= " << kID << " and ID = " << Proc->NLO << std::endl;
	  double res0=0., err0=0., chi0=0.;
	  Integrate(res0,err0,chi0,Proc);

	  // Writing in output file
	  if(kID==0) LO << res0 << "," << sc << std::endl;
	  else if(kID==1) NLO << res0 << "," << sc << std::endl;
	  else if(kID==2) Resummed << res0 << "," << sc << std::endl;
	  else if(kID==3) Expand << res0 << "," << sc << std::endl;
	  else std::cout << "ERROR in Process ID, k =  " << kID << std::endl;
	}
      
      //Increment of center of mass energy
      sc=sc*pow(sf/s0,1./nsteps);
      
    }

  Resummed.close();
  Expand.close();
  LO.close();
  NLO.close();
  // End of the program
  return 0;
}
