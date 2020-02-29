#ifndef MATRIX_ELEMENT_H
#define MATRIX_ELEMENT_H
#include "wfn.h"
double compute_matrix_element_sigma_tau_plus(char *density_file, int iv);
double compute_matrix_element_tau_plus(char* density_file, int iv);
double compute_matrix_element_M_GT(char *density_file);
double compute_matrix_element_M_F(char *density_file);
double compute_matrix_element_M_F1(char *density_file);
double compute_matrix_element_M_I2(char *density_file);
double compute_matrix_element_M_I1(char *density_file);
double compute_matrix_element_M_J2(char *density_file);
double compute_matrix_element_M_J1(char *density_file);
double compute_total_matrix_element_TT(char * density_file, int iv);

double compute_matrix_element_sigma_0_finite_q_op5(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, double qt, int J1, int J2, int t12);

double compute_matrix_element_sigma_0_finite_q(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, double qt, int J, int t12);
double compute_total_matrix_element_sigma_0_finite_q(char *density_file, double qt, int J);
double compute_total_matrix_element_sigma_0_finite_q_op5(char *density_file, double qt, int J1, int J2);

double compute_matrix_element_sigma_0(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, int iv);
double compute_matrix_element_TT(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, int iv);
double cme_1_sigma();
int get_l(int n, int j);
double compute_total_matrix_element_sigma_0(char *density_file);
double compute_matrix_element_y2_finite_q(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J1, int J2);
double compute_total_matrix_element_y2_finite_q(char* density_file, double q, int J1, int J2);
double compute_matrix_element_y2_finite_q_alt(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J1, int J2);
double compute_total_matrix_element_y2_finite_q_alt(char* density_file, double q, int J1, int J2);

#endif
