// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <cmath>      	// Mathematical functions       //
#include <sstream>    	// String streams               //
#include <string>	// Strings                      //
#include <iostream>
#include <fstream>
#include <complex>
// -------- Classes ----------------------------------- //
#include "process.h"      	// Process definition   //
#include <gsl/gsl_monte_vegas.h>// GSL                  //
#include "LHAPDF/LHAPDF.h"
// -------- Functions --------------------------------- //

using namespace LHAPDF;
//double Born(double*,size_t,void*);
//double Real(double*,size_t,void*);
//double Virtual(double*,size_t,void*);
double Tot(double*,size_t,void*);

void DisplayXsec(const double&,const double&,const std::string&); //
// ---------------------------------------------------- //

// ************************************************************************* //
//  Main integration routine, steering vegas.                                //
// ************************************************************************* //
void Integrate(double& res, double& err, double& chi, Process* proc)
{
// Initializing the random number generator 
   gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
   time_t seed = time(NULL);
   gsl_rng_set(r,seed);

    // Integrand
    gsl_monte_function I; I.dim=0;
    size_t calls=0;

    // Number of integrations and calls
    calls=10000;
    I.dim=2;


    // Selecting the good integrand + includes in the variable factor the
    // constant pieces of the squared matrix element

I.f= &Tot;
I.params=proc;

    // Integration bounds
    double xmin[I.dim], xmax[I.dim];
    for(size_t i=0; i<I.dim; i++) { xmin[i]=0.; xmax[i]=1.; } // should be called once

    // Initialization
    DisplayXsec(res,err,"-init");

    // Integration
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(I.dim);

    // Warm-up
    s->stage=0;  s->iterations=5;
    s->verbose=0;
    gsl_monte_vegas_integrate(&I, xmin, xmax, I.dim, calls/50, r, s, &res, &err);
    chi = s->chisq;
    DisplayXsec(res,err,"Warm-up"); 
    double prec=1e9;

    // Real things: stops if the precision reaches 1% or if one oscillates
    int counter=1; s->iterations=5;
    while(prec>1e-2 && counter<=10)
    {
       std::ostringstream ocnt; ocnt << counter; std::string cntstr= ocnt.str();
       s->stage=1;  
       gsl_monte_vegas_integrate(&I, xmin, xmax, I.dim, calls, r, s, &res, &err);
       prec=std::abs(err/res);
       DisplayXsec(res,err,"Refine-"+cntstr);
       counter++;
       calls*=5; 
    }
    chi = s->chisq;
    DisplayXsec(res,err,"final");

    // Cleaning the memory and closing the file
    gsl_monte_vegas_free(s);
    gsl_rng_free(r);
}



