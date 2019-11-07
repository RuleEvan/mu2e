#ifndef ROMBERG_H
#define ROMBERG_H
#include "glovar.h"
double RombergIntegrator(double (*f)(double), double a, double b, double tol); 
double Romberg2Vars(double (*f)(double, double), double a, double b, double r, double tol);
double Romberg3Vars(double (*f)(double, double, int), double a, double b, double p, int iv, double tol);
double RombergSpline2Integrator(double (*f)(gsl_spline*, gsl_spline*, gsl_interp_accel* , double), gsl_spline* g_spline, gsl_spline* f_spline, gsl_interp_accel* acc, double a, double b, double tol);

double RombergSplineFunIntegrator(double (*f)(double (*lep) (double), gsl_spline*, gsl_spline*, gsl_spline*, gsl_spline*, gsl_interp_accel* , gsl_interp_accel* , gsl_interp_accel* , gsl_interp_accel* , double), double (*lep) (double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel* acc1, gsl_interp_accel* acc2, gsl_interp_accel* acc3, gsl_interp_accel* acc4, double a, double b, double tol);


double RombergSpline4Integrator(double (*f)(gsl_spline*, gsl_spline*, gsl_spline*, gsl_spline*, gsl_interp_accel*, gsl_interp_accel*, gsl_interp_accel*, gsl_interp_accel*, double, double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel* acc1, gsl_interp_accel* acc2, gsl_interp_accel* acc3, gsl_interp_accel* acc4, double q, double a, double b, double tol);

double RombergSpline42Integrator(double (*f)(gsl_spline*, gsl_spline*, gsl_spline*, gsl_spline*, gsl_interp_accel*, gsl_interp_accel*, gsl_interp_accel*, gsl_interp_accel*, double, double, double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel* acc1, gsl_interp_accel* acc2, gsl_interp_accel* acc3, gsl_interp_accel* acc4, double q1, double q2, double a, double b, double tol);

double RombergSpline41Integrator(double (*f)(gsl_spline*, gsl_spline*, gsl_spline*, gsl_spline*, gsl_interp_accel*, gsl_interp_accel*, gsl_interp_accel*, gsl_interp_accel*, double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel* acc1, gsl_interp_accel* acc2, gsl_interp_accel* acc3, gsl_interp_accel* acc4, double a, double b, double tol);


#endif
