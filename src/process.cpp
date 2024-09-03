// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <cmath>      // Mathematical functions         //
#include <cstring>   // In/Out streams                 //
#include <sstream>    // String streams                 //
// -------- Classes ----------------------------------- //
//#include "gsl/gsl_matrix.h"  // GSL                     //
#include "messages.h"        // Message services        //
#include "process.h"         // Process                 //
#include "LHAPDF/LHAPDF.h"

using namespace LHAPDF;
// ---------------------------------------------------- //

// ************************************************************************* //
//  Constructor                                                              //
// ************************************************************************* //
//Process::Process(cst &s1, cst &s2, cst&s3, cst&s4)
Process::Process(cst &s1, const double&s2, const int s3)
{
   // Printing information
   info("Initializing the collision setup");

   // PDF set
   //pdf=s1;
   std::istringstream is1(s1); is1 >> pdf;

   // Center of mass energy
   sh=s2; sh*=sh;
   
   // Cross section ID=0 LO, 1 NLO, 2 Resummed, 3 Expanded
   NLO=s3;
   //std::istringstream is3(s3); is3 >> NLO;

   //Loading the PDF grids
   std::string name;
   if(pdf==1) name = "CT18NLO";
   else if(pdf==2) name = "CT14lo";
   else if (pdf==3) name = "CT14llo";
   else if (pdf==4) name = "cteq6";

   F=mkPDF(name,0);


}


// ************************************************************************* //
//  Initialization of the pdfs, couplings,...                                //
// ************************************************************************* //
void Process::Init()
{
     // Loading the PDF grids
   // std::string name;
   // if(pdf==1) name = "CT18NLO";
   // else if(pdf==2) name = "CT14lo";
   // else if (pdf==3) name = "CT14llo";
   // else if (pdf==4) name = "cteq6";
   // const PDF* F=mkPDF(name,0);

  // Loading the PDF grids
   std::string name;
   if(pdf==1) name = "pdf/i2Tn2.00.pds";
   else if(pdf==2) name = "pdf/CT14LN.pds";
   else if (pdf==3) name = "pdf/CT14LL.pds";
   else if (pdf==4) name = "pdf/cteq6m.tbl";
   char myname[40]; strcpy(myname,name.c_str());
   if(pdf<4) setct18_(myname);
   else {int dum = 1; setctq6_(dum);}
}

