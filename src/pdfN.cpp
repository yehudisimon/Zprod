// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2016 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Implements the PDF module and the PDF fit (needed for Mellin space PDFs).

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <complex>
//#include "gsl_all.h"
//#include "LHAPDF/LHAPDF.h"
//#include "utils.h"
//#include "pdf.h"
//#include "maths.h"
#include "process.h"

#include "gsl/gsl_blas.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_monte_vegas.h"
#include "gsl/gsl_multifit_nlin.h"
#include "gsl/gsl_sf.h"
#include "gsl/gsl_types.h"

using namespace std;

#define NFLVR 5   // Number of active flavours.
#define NDATA 1000 // Number of data points for the fit

void SetPDF(const int&, double&, double&, double*, double*, double&);


// function used to sample the PDFs.
double sampling(double xmin, double xmax,const size_t n, int i) {
    return pow(xmin, 1.0 - (1.0 + (double)i) / (xmax + (double)n));
    //return xmin + (double)i / (double)n * (xmax - xmin) - xmin;
    //double step_size = (xmax - xmin)/(double)n;
    //return(xmin + (double)i * step_size);

}


// Struct for GSL fit.
struct data {
    size_t n;
    double *y;
    double *sigma;
    double xr;
};

int expb_f(const gsl_vector *x, void *data, gsl_vector *f) {
    size_t n = ((struct data *)data)->n;
    double *y = ((struct data *)data)->y;
    double *sigma = ((struct data *)data)->sigma;
    double xr = ((struct data *)data)->xr;

    double A[8];
    for (size_t i0 = 0; i0 < 8; i0++) {
        A[i0] = gsl_vector_get(x, i0);
    }

    for (size_t i0 = 0; i0 < n; i0++) {
      //double t = pow(xr, 1.0 - (1.0 + (double)i0) / (1.0 + (double)n));
      double t = sampling(xr,1.0,n,i0);
        double Yi = A[0] * pow(t, A[1]) * pow(1.0 - t, A[2])
                    * (1.0 + A[3] * sqrt(t) + A[4] * t + A[5] * pow(t, 1.5)
                       + A[6] * pow(t, 2.0) + A[7] * pow(t, 2.5));
        gsl_vector_set(f, i0, (Yi - y[i0]) / sigma[i0]);
    }

    return GSL_SUCCESS;
}

int expb_df(const gsl_vector *x, void *data, gsl_matrix *J) {
    size_t n     = ((struct data *)data)->n;
    double *sigma = ((struct data *)data)->sigma;
    double xr    = ((struct data *)data)->xr;

    double A[8];
    for (size_t i0 = 0; i0 < 8; i0++) {
        A[i0] = gsl_vector_get(x, i0);
    }

    for (size_t i0 = 0; i0 < n; i0++) {
      //double t = pow(xr, 1.0 - (1.0 + (double)i0) / (1.0 + (double)n));
        double t = sampling(xr,1.0,n,i0);
        double s = sigma[i0];
        double e = A[0] * pow(t, A[1]) * pow(1.0 - t, A[2])
                   * (1.0 + A[3] * sqrt(t) + A[4] * t + A[5] * pow(t, 1.5)
                      + A[6] * pow(t, 2.0) + A[7] * pow(t, 2.5)) / s;
        gsl_matrix_set(J, i0, 0, e / A[0]);
        gsl_matrix_set(J, i0, 1, e * log(t));
        gsl_matrix_set(J, i0, 2, e * log(1.0 - t));

        e = A[0] * pow(t, A[1]) * pow(1.0 - t, A[2]) / s;
        gsl_matrix_set(J, i0, 3, e * sqrt(t));
        gsl_matrix_set(J, i0, 4, e * t);
        gsl_matrix_set(J, i0, 5, e * pow(t, 1.5));
        gsl_matrix_set(J, i0, 6, e * pow(t, 2.0));
        gsl_matrix_set(J, i0, 7, e * pow(t, 2.5));
    }

    return GSL_SUCCESS;
}

int expb_fdf(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J) {
    expb_f(x, data, f);
    expb_df(x, data, J);

    return GSL_SUCCESS;
}

void Fit(double A[8], double E[8], int flag, double xr, double muf, double weight_valence, double weight_sea, double weight_gluon, Process *proc) {
    // Defines the function to minimize.
    const size_t n = NDATA;
    const size_t p = 8;
    
    double y[NDATA], sigma[NDATA];
    struct data d = { n, y, sigma, xr };

    //std::cout << "Entering Fit" << std::endl;
    gsl_multifit_function_fdf f;
    f.f = &expb_f;
    f.df = &expb_df;
    f.fdf = &expb_fdf;
    f.n = n;
    f.p = p;
    f.params = &d;

    // PDFs to fit.
    double gamma_weight = -1.6;
    for (size_t i0 = 0; i0 < n; i0++)
      {
	//double t = pow(xr, 1.0 - (1.0 + (double)i0) / (1.0 + (double)n));
	double t = sampling(xr,1.0,n,i0);
	double q[2][5];
	double g;
	double qq;
	//std::cout << "Setting the PDF in x space" << std::endl;
	SetPDF(proc->pdf,t,muf,q[0],q[1],g);
	//std::cout << "PDF set" << std::endl;
	//pdfX(g, q, t, Q2);
	switch (flag) {
	case 0:
	  qq = g;
	  gamma_weight = weight_gluon;
	  break;
	case 1:
	  qq = q[0][0]; //d
	  gamma_weight = weight_valence;
	  break;
	case 2:
	  qq = q[0][1]; //u
	  gamma_weight = weight_valence;
	  break;
	case 3:
	  qq = q[1][0]; //dbar
	  gamma_weight = weight_sea;
	  break;
	case 4:
	  qq = q[1][2]; //sbar or s
	  gamma_weight = weight_sea;
	  break;
	case 5:
	  qq = q[1][4]; //bbar or b
	  gamma_weight = weight_sea;
	  break;
	case 6:
	  qq = q[1][1]; //ubar
	  gamma_weight = weight_sea;
	  break;
	case 7:
	  qq = q[1][3]; //cbar or c
	  gamma_weight = weight_sea;
	  break;
	default:
	  fprintf(stderr,
		  "error: while retrieving PDF: unkown flag %d\n", flag);
	  exit(1);
	}

	//std::cout << "Allocating the current PDF" << std::endl;
        y[i0] = qq;
        sigma[i0] = pow(t,gamma_weight); // DeFlorian-like pdf-weights t^(-1.6); DeJonathan t^(-1)
    }

    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc(T, n, p);


    // First guess for the parameters.
    double x_init[8] = {1.0, -1.4, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    if (flag == 1 || flag == 2) {
        x_init[1] = -0.6;
    }

    gsl_vector_view x = gsl_vector_view_array(x_init, p);
    gsl_multifit_fdfsolver_set(s, &f, &x.vector);

    int status = GSL_CONTINUE;
    unsigned int iter = 0;

    //std::cout << "Doing the actual fit" << std::endl;
    // Does the fit.
    while (status == GSL_CONTINUE && iter < 500000) {
        iter++;
        status = gsl_multifit_fdfsolver_iterate(s);
        if (status) {
            break;
        }
        status = gsl_multifit_test_delta(s->dx, s->x, 1.0e-12, 1.0e-5);
    }
    //std::cout << "Exiting fit" << std::endl;

    gsl_matrix *covar = gsl_matrix_alloc(p, p);
    gsl_multifit_covar(s->J, 0.0, covar); // J ?
    //double c = pow(gsl_blas_dnrm2(s->f), 2) / ((double)(n - p));
    //c = (c > 1.0 ? c : 1.0);
    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof));
    
    // Gets parameters.
    for (size_t i0 = 0; i0 < 8; i0++) {
        A[i0] = 0.0;
    }
    for (size_t i0 = 0; i0 < p; i0++) {
        A[i0] = gsl_vector_get(s->x, i0);
    }
    for (size_t i0 = 0; i0 < p; i0++) {
        E[i0] = c * gsl_matrix_get(covar, i0, i0);
    }

    
    printf("#chisq/dof = %g\n",  pow(chi, 2.0) / dof);
    printf ("#status = %s\n", gsl_strerror (status));

    gsl_multifit_fdfsolver_free(s);
    gsl_matrix_free(covar);
}

void pdfFit(double &A1MIN, double A[8][8], double &muf, double weight_valence, double weight_sea, double weight_gluon, Process *p) {
    const int nf = NFLVR;
    fprintf(stderr,"Performing PDF fit with %d flavors, Q^2 = mu_F^2 = %g\n and weights: valence: x^%g, sea: x^%g, gluon: x^%g and xmin = 1e-7 \n Fit function: f = A0 * x^A1 * (1 - x)^A2 * ( 1 + A3 * x^(1/2) + A4 * x + A5 * x^(3/2) + A6 * x^2 + A7 * x^(5/2) )\n", nf, muf, weight_valence, weight_sea, weight_gluon);

    for (int i0 = 0; i0 < 8; i0++) {
        double err[8];
        switch (i0) {
        case 0:
            fprintf(stderr, "Fitting gluon PDF...");
            break;
        case 1:
            fprintf(stderr, "Fitting valence down quark PDF...");
            break;
        case 2:
            fprintf(stderr, "Fitting valence up quark PDF...");
            break;
        case 3:
            fprintf(stderr, "Fitting sea down quark PDF...");
            break;
        case 4:
            fprintf(stderr, "Fitting strange quark PDF...");
            break;
        case 5:
            fprintf(stderr, "Fitting bottom quark PDF...");
            break;
        case 6:
            fprintf(stderr, "Fitting sea up quark PDF...");
            break;
        case 7:
            fprintf(stderr, "Fitting charm quark PDF...");
            break;
        }
        fflush(stderr);
	//std::cout << "Arriving to Fit" << std::endl;
        //Fit(A[i0], err, i0, xmin, muf); // xmin = M^2/Smax >= 1e-7
        Fit(A[i0], err, i0, 1e-7, muf, weight_valence,  weight_sea,  weight_gluon,p); 
        fprintf(stderr, " done.\n");
        fprintf(stderr, "Fit result:\n");

        for (int i1 = 0; i1 < 8; i1++) {
          fprintf(stderr, "A%i =  %.5f # +-%.5f\n",i1, A[i0][i1], err[i1]);
        }
        if (A[i0][1] < A1MIN) {
            A1MIN = A[i0][1];
        }
    }
}
