#ifndef PROCESS_H
#define PROCESS_H
// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <string>    // String streams                  //
#include "LHAPDF/LHAPDF.h"
// -------- Processor variables ----------------------- //
using namespace LHAPDF;
#define cst const std::string                           //
// -------- Functions --------------------------------- //
extern "C" { void setct18_(char[40]); void setctq6_(int&); };                    //
// ---------------------------------------------------- //

// ************************************************************************* //
// Definition of the class Process                                           //
// ************************************************************************* //
class Process
{
public:
  Process()  { };
  //Process(cst&,cst&,cst&);
  Process(cst &, const double&, const int);
  ~Process() { };
  
  // Kinematical quantities and scales
  double sh;               // Hadronic center of mass energy
  int NLO;              // Cross section ID 0=LO, 1=NLO, 2=Resummed, 3=Expanded 
  int pdf;              // PDF identifier

  const PDF* F;
  
  // Functions
  void Init();          // Init of PDFs, etc...
};
#endif
