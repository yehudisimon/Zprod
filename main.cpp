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

  const std::string & pathstr="/home/yehudi/LHAPDF6/share/LHAPDF";

  setPaths(pathstr);
  for (const string& p : paths())
    cout << p << endl;
  
  cout << "@" << findFile("lhapdf.conf") << "@" << endl;
  cout << "List of available PDFs:" << endl;
  for (const string& s : availablePDFSets())
    cout << " " << s << endl;
      
 
  // std::cout << "making PDF" << std::endl;
  //const PDF* pdf = mkPDF("CT18NLO", 0);
    
  //Opening writing file
  std::ofstream LO;
  std::ofstream NLO;
  std::ofstream Expand;
  std::ofstream Resummed;
  std::ofstream ExpandHS;
  std::ofstream ResummedHS;
  std::ofstream ExpandHSConvert;
  std::ofstream ResummedHSConvert;
  
  // Initializing message services
  InitMessages();

  // Setup of the collision setup, pdfs, scales, etc...
  const int nsteps=10, j0=0;//, pdf=int(argv[1]);
  double s0=300., sf=1e6, sc=s0;
  //const char *J;
  //J = std::to_string(j0).c_str();
  
  
  Process *Proc = new Process(argv[1], sc, j0);
  

  Proc->Init(); //PDF initialisation
  
  //Fitting PDFs to get Ai coefficients to include in the F_j(N) Mellin transform of PDFs
  //if(Proc->NLO==2 || Proc->NLO==3) pdfFit(A1min,A,muF,-1.6,-1.6,-1.6,Proc); 

  // Only for CT18NLO
  //Fitting gluon PDF...#chisq/dof = 5.08633e-05
  int i0=0;
  A[i0][0] =  3.36311 ;
  A[i0][1] =  -1.35686 ;
  A[i0][2] =  12.01851 ;
  A[i0][3] =  -9.06557;
  A[i0][4] =  59.82888 ;
  A[i0][5] =  -209.61905;
  A[i0][6] =  358.78810 ;
  A[i0][7] =  -214.30364;

  //Fitting valence down quark PDF...#chisq/dof = 1.03978e-07
  i0=1;
  A[i0][0] =  0.17088;
  A[i0][1] =  -1.32599;
  A[i0][2] =  4.86054;
  A[i0][3] =  -2.74299;
  A[i0][4] =  27.12173 ;
  A[i0][5] =  -51.56600;
  A[i0][6] =  58.71166;
  A[i0][7] =  -29.80953;

  //Fitting valence up quark PDF...#chisq/dof = 1.29063e-07
  i0=2;
  A[i0][0] = 0.15891;
  A[i0][1] = -1.33284;
  A[i0][2] = 4.04498 ;
  A[i0][3] =  -1.28483 ;
  A[i0][4] = 25.16123;
  A[i0][5] =  -21.03617;
  A[i0][6] =  37.34074 ;
  A[i0][7] =  -32.12262 ;

  //Fitting sea down quark PDF...#chisq/dof = 7.17778e-08
  i0=3;
  A[i0][0] = 0.18060;
  A[i0][1] =  -1.32147;
  A[i0][2] =  4.80454 ;
  A[i0][3] =  -4.01851;
  A[i0][4] =  19.31558 ;
  A[i0][5] =   -55.90824;
  A[i0][6] =  71.69672;
  A[i0][7] =  -33.26465;

  //strange quark
  i0=4;
  A[i0][0] = 0.19715;
  A[i0][1] = -1.31408;
  A[i0][2] = 12.42845 ;
  A[i0][3] =  -7.36475 ;
  A[i0][4] =  44.88323;
  A[i0][5] =  -146.11711;
  A[i0][6] = 243.91655;
  A[i0][7] = -145.10893;

  //Fitting bottom quark PDF...#chisq/dof = 1.7599e-08
  i0=5;
  A[i0][0] = 0.10914;
  A[i0][1] =  -1.34332;
  A[i0][2] =  16.06492;
  A[i0][3] =   -9.85651;
  A[i0][4] =   68.81877;
  A[i0][5] = -260.94287;
  A[i0][6] = 489.90364;
  A[i0][7] = -321.74473;

  //Fitting sea up quark PDF...#chisq/dof = 8.35594e-09
  i0=6;
  A[i0][0] =  0.19950;
  A[i0][1] =  -1.31377;
  A[i0][2] =  10.66072;
  A[i0][3] =  -5.55560;
  A[i0][4] =  37.16200;
  A[i0][5] = -128.66354;
  A[i0][6] = 217.17546;
  A[i0][7] = -124.50014;

  //Fitting charm quark PDF...#chisq/dof = 1.61266e-08
  i0=7;
  A[i0][0] =  0.14115;
  A[i0][1]=  -1.33602;
  A[i0][2] = -0.27091;
  A[i0][3] =  -6.04545;
  A[i0][4] = 15.21004;
  A[i0][5] = -19.76661;
  A[i0][6] =  13.13392;
  A[i0][7] =  -3.53289;

  LO.open("LOXsec.txt");
  NLO.open("NLOXsec.txt");
  Resummed.open("ResumXsec.txt");
  Expand.open("ExpandedXsec.txt");

  ResummedHS.open("ResumXsecHS.txt");
  ExpandHS.open("ExpandedXsecHS.txt");
  
  ResummedHSConvert.open("ResumXsecHSConvert.txt");
  ExpandHSConvert.open("ExpandedXsecHSConvert.txt");
  

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
	  // else if(kID==4) ResummedHS << res0 << "," << sc << std::endl;
	  // else if(kID==5) ExpandHS << res0 << "," << sc << std::endl;
	  // else if(kID==6) ResummedHSConvert << res0 << "," << sc << std::endl;
	  // else if(kID==7) ExpandHSConvert << res0 << "," << sc << std::endl;
	  else std::cout << "ERROR in Process ID, k =  " << kID << std::endl;
	}
      
      //Increment of center of mass energy
      sc=sc*pow(sf/s0,1./nsteps);
      
    }

  Resummed.close();
  Expand.close();
  ResummedHS.close();
  ExpandHS.close();
  ResummedHSConvert.close();
  ExpandHSConvert.close();
  LO.close();
  NLO.close();
  // End of the program
  return 0;
}
