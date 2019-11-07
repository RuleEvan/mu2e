#ifndef WFN_H
#define WFN_H
#include "charge.h"

void mu2e();
double wfn_sq(gsl_spline* g_spline, gsl_spline* f_spline, gsl_interp_accel *acc, double r);
double g1lep(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r, double q);
double g1lep2(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r, double q1, double q2);
double M_op_sd(int ijp, int ij, int j_op, double q, double b);
double b_osc(int a_nuc);

double newlep1(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r);

double newlep2(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r);

double newlep3(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r);

double core1(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r);

double core2(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r);

double core3(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r);

double MJ0_s1_s1(double x);
double MJ0_d3_d3(double x);
double MJ0_d5_d5(double x);


double SigmaJ0_1p1_1s1(double x);
double SigmaJ0_2s1_1p1(double x);
double SigmaJ0_1d3_1p3(double x);

double SigmaJ1_2s1_1s1(double x);
double SigmaJ1_2s1_2s1(double x);
double SigmaJ1_1d3_2s1(double x);
double SigmaJ1_1d3_1s1(double x);
double SigmaJ1_1d3_1d3(double x);
double SigmaJ1_1d5_1d3(double x);
double SigmaJ1_1d5_1d5(double x);
double SigmaJ1_1s1_1s1(double x);
double SigmaJ1_1p1_1p1(double x);
double SigmaJ1_1p3_1p1(double x);
double SigmaJ1_1p3_1p3(double x);

double SigmaJ2_1p3_1s1(double x);
double SigmaJ2_2s1_1p3(double x);
double SigmaJ2_1d3_1p1(double x);
double SigmaJ2_1d3_1p3(double x);
double SigmaJ2_1d5_1p1(double x);
double SigmaJ2_1d5_1p3(double x);
double SigmaJ2_1d5_1d5(double x);

double SigmaJ3_1p3_1p3(double x);
double SigmaJ3_1d5_1s1(double x);
double SigmaJ3_1d3_1d3(double x);
double SigmaJ3_1d5_2s1(double x);
double SigmaJ3_1d5_1d3(double x);
double SigmaJ3_1d5_1d5(double x);

double SigmaJ4_1d5_1p3(double x);

double SigmaJ5_1d5_1d5(double x);

double lep_int_w1(double (*f) (double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r);


#endif
