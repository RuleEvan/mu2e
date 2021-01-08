#ifndef MATRIX_ELEMENT_H
#define MATRIX_ELEMENT_H
#include "wfn.h"

// Finite q Spline functions
double compute_total_matrix_element_finite_q_op1(char* density_file, double q, int J);
double compute_total_matrix_element_finite_q_op1_alt(char* density_file, double q, int J);

double compute_total_matrix_element_finite_q_op3(char* density_file, double q, int J1, int J2);
double compute_total_matrix_element_finite_q_op4(char* density_file, double q, int J1, int J2);
double compute_total_matrix_element_finite_q_op5(char *density_file, double qt, int J1, int J2);
double compute_total_matrix_element_finite_q_op6(char *density_file, double qt, int J1, int J2, int J);
double compute_total_matrix_element_finite_q_op7(char *density_file, double qt, int J1, int J2, int J);


double compute_matrix_element_finite_q_op1(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J, gsl_spline *f_spline, gsl_interp_accel *acc);
double compute_matrix_element_finite_q_op3(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J1, int J2, gsl_spline *f_spline, gsl_interp_accel *acc);
double compute_matrix_element_finite_q_op4(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J1, int J2, gsl_spline *f_spline, gsl_interp_accel *acc);
double compute_matrix_element_finite_q_op5(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double qt, int J1, int J2, gsl_spline *f_spline, gsl_interp_accel *acc);
double compute_matrix_element_finite_q_op6(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double qt, int J1, int J2, int J, gsl_spline* f_spline, gsl_interp_accel *acc);
double compute_matrix_element_finite_q_op7(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double qt, int J1, int J2, int J, gsl_spline* f_spline, gsl_interp_accel *acc);

// q = 0 functions
double compute_matrix_element_sigma_0(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, int iv);
double compute_matrix_element_TT(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, int iv);
double compute_matrix_element_I2(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12);
double compute_relative_density(double r, int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12);

double compute_total_matrix_element_sigma_0(char *density_file);
double compute_total_relative_density(char *density_file, double r);
double compute_total_matrix_element_TT(char * density_file, int iv);
double compute_total_matrix_element_I2(char * density_file);
double compute_total_matrix_element_I2_alt(char * density_file);
double compute_total_matrix_element_I2_cj(char * density_file);
void translate_cj_to_er(char* density_file);
void generate_bigstick_int_file();

int get_l(int n, int j);
int get_n(int index);
int get_j(int index);

void convert_density_file(char *density_file);

#endif
