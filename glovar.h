#include <stdio.h>
#include "math.h"
#include "stdlib.h"
#include "stdint.h"
#include "string.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf.h>

// Isotope setup
#define A_NUC 76 // Atomic Mass
#define A_FACTOR 9.155 // [MeV] Average nuclear excitation energy
#define B_OSC 2.927
#define Z_ATOM 32 // Atomic Number

// Technical parameters
#define COR_FAC 1

// Physical constants
#define ALPHA_FS 0.007297352
#define R_NUC 1.2 // [fm]
#define M_ELECTRON 0.511 // [MeV]
#define M_NEUTRON 939.57 // [MeV]
#define E_BETA 2.55 // [MeV]
#define G_AXIAL 1.275 // Axial coupling
#define G_MAGNETIC 4.7
#define G_TENSOR 0.99
#define G_VECTOR 1.0

#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) > (b) ? (a):(b))

#define PION_MASS 139.570 //[MeV]
#define LAMBDA_V 850.0 //[MeV]
#define LAMBDA_A 1040.0 // [MeV]
#define KAPPA_1 3.7
#define MU 1
