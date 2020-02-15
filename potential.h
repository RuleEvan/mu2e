#ifndef POTENTIAL_H
#define POTENTIAL_H
#include "angular.h"

double finite_q_alpha_pot_1(double r, int l, double q, double m_pi);
double finite_q_alpha_int_1(double r, int l, double q, double m_pi, double alpha);
double compute_rel_potential(double np, double lp, double n, double l, int J, double qt, int iv);
double talmi_rel(double p, int iv, int J, double qt);
double talmi_integrand_rel(double q, int iv, double p, double J, double qt);
double v_cm_finite_q(double r, int l, double q);

double b_coeff(double n, double l, double np, double lp, double p);
double compute_potential(double n, double np, double l, double lp, int iv);
double a_coeff(double n, double l, double k);
double v_light_limit(double q);
double v_pion_f1(double r);
double v_pion_f2(double r);
double v_light_limit_d(double r);
double v_pion_g1(double r);
double v_pion_g2(double r);
double v_pion_NN_g1(double r);
double v_pion_NN_g2(double r);
double talmi(double p, int iv);
double talmi_integrand(double q, double p, int iv);
double h_MM_GT(double r);
double h_MM_T(double r);
double h_AA_GT_sd(double r);
double h_AA_T_sd(double r);
double h_AP_GT_sd(double r);
double h_AP_T_sd(double r);
double h_PP_GT_sd(double r);
double h_PP_T_sd(double r);
double h_MM_GT_sd(double r);
double h_MM_T_sd(double r);
double h_F_sd(double r);
double h_AA_GT(double r);
double h_AA_T(double r);
double h_AP_GT(double r);
double h_AP_T(double r);
double h_PP_GT(double r);
double h_PP_T(double r);
double h_F(double r);
#endif
