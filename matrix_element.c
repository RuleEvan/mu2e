#include "matrix_element.h"

/* Two-body matrix elements */

double compute_matrix_element_finite_q_op1(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J, gsl_spline *f_spline, gsl_interp_accel *acc) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double t12 = it12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  double m4 = 0.0;
  // Lambda = Lambdap and S = SP
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    int s_max = MIN(lambda + j12, 1);
    int s_min = abs(lambda - j12);
    for (int s = s_min; s <= s_max; s++) {
      // JJ -> LS coupling factors
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      // Un-reduce matrix element wrt total J
      fact *= sqrt(2*j12 + 1.0);
      // Un-reduced matrix element of sig(1) dot sig(2)
      double m1 = pow(-1.0, 1.0 + s)*six_j(s,0.5,0.5,1.0,0.5,0.5)*6.0;
      // Reduced matrix element of tau(1) dot tau(2)
      m1 *= pow(-1.0, 1.0 + t12)*sqrt(2*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*6.0;

      m1 *= fact;
      m1 *= compute_radial_matrix_element_finite_q_op1(J, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12, q, f_spline, acc);
      m4 += m1;
    }
  }

  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= 1/sqrt(2);}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= 1/sqrt(2);}

  return m4;
}

double compute_matrix_element_finite_q_op3(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J1, int J2, gsl_spline *f_spline, gsl_interp_accel *acc) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;
  double t12 = it12/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
    
  double m4 = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  int lambda_min = MAX(abs(l1 - l2), abs(j12 - 1));
  int lambda_max = MIN(l1 + l2, j12 + 1);
  int lambdap_min = MAX(abs(l1p - l2p), abs(j12p - 1));
  int lambdap_max = MIN(l1p + l2p, j12p + 1);
  for (int lambda = lambda_min; lambda <= lambda_max; lambda++) {
    lambdap_min = MAX(lambdap_min, lambda - 2);
    lambdap_max = MIN(lambdap_max, lambda + 2);
    if (lambdap_min > lambdap_max) {continue;}
    for (int lambdap = lambdap_min; lambdap <= lambdap_max; lambdap++) {
      // JJ -> LS coupling factors
      // S = SP = 1 required by spin tensor operator
      double fact = sqrt(3.0*(2*lambda + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt(3.0*(2*lambdap + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, 1, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, 1, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      // Un-reduce wrt total J
      fact *= sqrt(2*j12 + 1);
      // Unreduced decoupling of Y2 dot S2 matrix element
      fact *= pow(-1.0, lambda + 1 + j12)*six_j(j12, 1, lambdap, 2, lambda, 1);
      // Reduced matrix element of tau(1) dot tau(2)
      fact *= pow(-1.0, 1.0 + t12)*sqrt(2.0*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*6.0;
      // Reduced matrix element of [sig(1) x sig(2)]_2
      fact *= 2.0*sqrt(5);
     
      fact *= compute_radial_matrix_element_finite_q_op3(J1, J2, n1p, l1p, n2p, l2p, lambdap, n1, l1, n2, l2, lambda, 1, t12, q, f_spline, acc);
      m4 += fact;
    }
  }
  if ((in1 == in2) && (j1 == j2)) {m4 *= 1/sqrt(2);}
  if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1/sqrt(2);}  

  return m4;
}

double compute_matrix_element_finite_q_op4(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J1, int J2, gsl_spline *f_spline, gsl_interp_accel *acc) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;
  double t12 = it12/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
    
  double m4 = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  int lambda_min = MAX(abs(l1 - l2), abs(j12 - 1));
  int lambda_max = MIN(l1 + l2, j12 + 1);
  int lambdap_min = MAX(abs(l1p - l2p), abs(j12p - 1));
  int lambdap_max = MIN(l1p + l2p, j12p + 1);
  for (int lambda = lambda_min; lambda <= lambda_max; lambda++) {
    lambdap_min = MAX(lambdap_min, lambda - 2);
    lambdap_max = MIN(lambdap_max, lambda + 2);
    if (lambdap_min > lambdap_max) {continue;}
    for (int lambdap = lambdap_min; lambdap <= lambdap_max; lambdap++) {
      //  JJ -> LS coupling factors
      double fact = sqrt(3.0*(2*lambda + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt(3.0*(2*lambdap + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, 1, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, 1, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      // Unreduce wrt total J
      fact *= sqrt(2*j12 + 1);
      // Unreduced decoupling of Y2 dot S2
      fact *= pow(-1.0, lambda + j12 + 1)*six_j(j12, 1, lambdap, 2, lambda, 1);
      // Reduced matrix element of S2
      fact *= 2.0*sqrt(5);
      // Reduced matrix element of tau(1) dot tau(2)
      fact *= pow(-1.0, 1.0 + t12)*sqrt(2.0*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*6.0;

      fact *= compute_radial_matrix_element_finite_q_op4(J1, J2, n1p, l1p, n2p, l2p, lambdap, n1, l1, n2, l2, lambda, 1, t12, q, f_spline, acc);
      m4 += fact;
    }
  }
  if ((in1 == in2) && (j1 == j2)) {m4 *= 1/sqrt(2);}
  if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1/sqrt(2);}  

  return m4;
}

double compute_matrix_element_finite_q_op5(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J1, int J2, gsl_spline *f_spline, gsl_interp_accel *acc) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double t12 = it12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;
//  double t12p = it12p/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  double m4 = 0.0;
  // Lambda = Lambdap and S = Sp
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    int s_max = MIN(lambda + j12, 1);
    int s_min = abs(lambda - j12);
    for (int s = s_min; s <= s_max; s++) {
      // JJ -> LS coupling factors
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      // Uncouple wrt total J
      fact *= sqrt(2*j12 + 1.0);
      // Un-reduced matrix element of sigma(1) dot sigma(2)
      double m1 = pow(-1.0, 1.0 + s)*six_j(s,0.5,0.5,1.0,0.5,0.5)*6.0;
      // Reduced matrix element of tau(1) dot tau(2)
      m1 *= pow(-1.0, 1.0 + t12)*sqrt(2*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*6.0;
      m1 *= fact;
      m1 *= compute_radial_matrix_element_finite_q_op5(J1, J2, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12, q, f_spline, acc);
      m4 += m1;
    }
  }
  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= 1/sqrt(2);}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= 1/sqrt(2);}

  return m4;
}


double compute_matrix_element_finite_q_op6(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J1, int J2, int J, gsl_spline *f_spline, gsl_interp_accel *acc) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double t12 = it12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;
//  double t12p = it12p/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  double m4 = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  int lambda_min = MAX(abs(l1 - l2), abs(j12 - 1));
  int lambda_max = MIN(l1 + l2, j12 + 1);
  int lambdap_min = MAX(abs(l1p - l2p), abs(j12p - 1));
  int lambdap_max = MIN(l1p + l2p, j12p + 1);
  for (int lambda = lambda_min; lambda <= lambda_max; lambda++) {
    lambdap_min = MAX(lambdap_min, lambda - 2);
    lambdap_max = MIN(lambdap_max, lambda + 2);
    if (lambdap_min > lambdap_max) {continue;}
    for (int lambdap = lambdap_min; lambdap <= lambdap_max; lambdap++) {
      int s_max = MIN(lambda + j12, 1);
      int s_min = abs(lambda - j12);
      for (int s = s_min; s <= s_max; s++) {
        int sp_max = MIN(lambdap + j12p, 1);
        int sp_min = abs(lambdap - j12p);
        for (int sp = sp_min; sp <= sp_max; sp++) {
	  if (s == sp) {continue;}
	  // JJ -> LS coupling factors
          double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
          fact *= sqrt((2*lambdap + 1)*(2*sp + 1)*(2*j1p + 1)*(2*j2p + 1));
          fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
          fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, sp, j1p, j2p, j12p);
	  // Unreduce wrt total J
          fact *= sqrt(2*j12 + 1.0);
	  // Un-reduced decoupling of Y1 dot S1
	  fact *= pow(-1.0, lambda + sp + j12)*six_j(j12, sp, lambdap, 1, lambda, s);
          if (fact == 0.0) {continue;}
	  // Reduced matrix element of S1
          double m1 = sqrt(6.0);
	  // Reduced matrix element of tau(1) dot tau(2)
          m1 *= pow(-1.0, 1.0 + t12)*sqrt(2*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*6.0;
          m1 *= fact;
          m1 *= compute_radial_matrix_element_finite_q_op6(J1, J2, J, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12, q, f_spline, acc);
          m4 += m1;
	}
      }
    }
  }
  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= 1/sqrt(2);}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= 1/sqrt(2);}

  return m4;
}

double compute_matrix_element_finite_q_op7(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J1, int J2, int J, gsl_spline *f_spline, gsl_interp_accel *acc) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;
  double t12 = it12/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
    
  double m4 = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  int lambda_min = MAX(abs(l1 - l2), abs(j12 - 1));
  int lambda_max = MIN(l1 + l2, j12 + 1);
  int lambdap_min = MAX(abs(l1p - l2p), abs(j12p - 1));
  int lambdap_max = MIN(l1p + l2p, j12p + 1);
  for (int lambda = lambda_min; lambda <= lambda_max; lambda++) {
    lambdap_min = MAX(lambdap_min, lambda - 2);
    lambdap_max = MIN(lambdap_max, lambda + 2);
    if (lambdap_min > lambdap_max) {continue;}
    for (int lambdap = lambdap_min; lambdap <= lambdap_max; lambdap++) {
      // JJ -> LS coupling factors
      double fact = sqrt(3.0*(2*lambda + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt(3.0*(2*lambdap + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, 1, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, 1, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      // Unreduce wrt total J
      fact *= sqrt(2*j12 + 1);
      // Unreduced decoupling of Y2 dot S2
      fact *= pow(-1.0, lambda + j12 + 1)*six_j(j12, 1, lambdap, 2, lambda, 1);
      // Reduced matrix element of tau(1) dot tau(2)
      fact *= pow(-1.0, 1.0 + t12)*sqrt(2.0*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*6.0;
      // Reduced matrix element of S2
      fact *= 2.0*sqrt(5);

      fact *= compute_radial_matrix_element_finite_q_op7(J1, J2, J, n1p, l1p, n2p, l2p, lambdap, n1, l1, n2, l2, lambda, 1, t12, q, f_spline, acc);
      m4 += fact;
    }
  }
  if ((in1 == in2) && (j1 == j2)) {m4 *= 1/sqrt(2);}
  if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1/sqrt(2);}  

  return m4;
}

double compute_total_matrix_element_finite_q_op1(char* density_file, double q, int J) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;
 
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *f_spline = gsl_spline_alloc(gsl_interp_linear, NSPLINE + 1);
  double delta_r = (RMAX - RMIN)/(1.0*NSPLINE);
  double *f_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *r_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));

  for (int i = 0; i <= NSPLINE; i++) {
    double r = RMIN + i*delta_r;
    r_array[i] = r;
    f_array[i] = finite_q_alpha_pot_1(r, J, q, M_PION);
//    printf("%g, %g, %g\n", q, r_array[i], f_array[i]);
  } 

  gsl_spline_init(f_spline, r_array, f_array, NSPLINE + 1);

  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {
    double m4 = compute_matrix_element_finite_q_op1(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, q, J, f_spline, acc); 
    mat += m4*density;
  }

  mat *= -1.0/(12.0*sqrt(M_PI));
 
  gsl_spline_free(f_spline);
  gsl_interp_accel_free(acc);  

  return mat;
}

double compute_total_matrix_element_finite_q_op3(char* density_file, double q, int J1, int J2) {
  if ((J1 % 2) || (J2 % 2)) {return 0.0;}
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *f_spline = gsl_spline_alloc(gsl_interp_linear, NSPLINE + 1);
  double delta_r = (RMAX - RMIN)/(1.0*NSPLINE);
  double *f_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *r_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));

  for (int i = 0; i <= NSPLINE; i++) {
    double r = RMIN + i*delta_r;
    r_array[i] = r;
    f_array[i] = finite_q_alpha_pot_2(r, J1, q, M_PION);
  } 

  gsl_spline_init(f_spline, r_array, f_array, NSPLINE + 1);

  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    double m4 = compute_matrix_element_finite_q_op3(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, q, J1, J2, f_spline, acc); 
    mat += m4*density;
  }
  
  int phase = pow(-1, (J1 - J2)/2.0); 
  mat *= phase*sqrt((2.0*J1 + 1.0)*(2.0*J2 + 1.0))*clebsch_gordan(J1, J2, 2.0, 0.0, 0.0, 0.0)/(5.0*sqrt(24.0*M_PI));
  
  gsl_spline_free(f_spline);
  gsl_interp_accel_free(acc);  
   
  return mat;
}

double compute_total_matrix_element_finite_q_op4(char* density_file, double q, int J1, int J2) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  if ((J1 % 2) || (J2 % 2)) {return 0.0;}
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;
 
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *f_spline = gsl_spline_alloc(gsl_interp_linear, NSPLINE + 1);
  double delta_r = (RMAX - RMIN)/(1.0*NSPLINE);
  double *f_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *r_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));

  for (int i = 0; i <= NSPLINE; i++) {
    double r = RMIN + i*delta_r;
    r_array[i] = r;
    f_array[i] = finite_q_alpha_pot_3(r, J2, q, M_PION);
  } 

  gsl_spline_init(f_spline, r_array, f_array, NSPLINE + 1);

  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    double m4 = compute_matrix_element_finite_q_op4(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, q, J1, J2, f_spline, acc); 
    mat += m4*density;
  }
  
  mat *= sqrt(2.0*J2 + 1.0)*clebsch_gordan(J2, 2, J1, 0.0, 0.0, 0.0)*sqrt(8.0*M_PI/15.0)/(8.0*M_PI);
  
  gsl_spline_free(f_spline);
  gsl_interp_accel_free(acc);  
                       
  return mat;
}

double compute_total_matrix_element_finite_q_op5(char* density_file, double q, int J1, int J2) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *f_spline = gsl_spline_alloc(gsl_interp_linear, NSPLINE + 1);
  double delta_r = (RMAX - RMIN)/(1.0*NSPLINE);
  double *f_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *r_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));

  for (int i = 0; i <= NSPLINE; i++) {
    double r = RMIN + i*delta_r;
    r_array[i] = r;
    f_array[i] = finite_q_alpha_pot_4(r, J1, q, M_PION);
  } 

  gsl_spline_init(f_spline, r_array, f_array, NSPLINE + 1);


  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {
    double m4 = compute_matrix_element_finite_q_op5(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, q, J1, J2, f_spline, acc); 
    mat += m4*density;
  }
  
  int phase = pow(-1.0, (1.0 + J1 - J2)/2.0);
  mat *= phase/(12.0*sqrt(M_PI))*(2.0*J1 + 1.0)/(2.0*J2 + 1.0)*pow(clebsch_gordan(J1, 1.0, J2, 0.0, 0.0, 0.0), 2.0);   

  gsl_spline_free(f_spline);
  gsl_interp_accel_free(acc);  
                      
  return mat;
}


double compute_total_matrix_element_finite_q_op6(char* density_file, double q, int J1, int J2, int J) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *f_spline = gsl_spline_alloc(gsl_interp_linear, NSPLINE + 1);
  double delta_r = (RMAX - RMIN)/(1.0*NSPLINE);
  double *f_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *r_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));

  for (int i = 0; i <= NSPLINE; i++) {
    double r = RMIN + i*delta_r;
    r_array[i] = r;
    f_array[i] = finite_q_alpha_pot_5(r, J1, q, M_PION);
  } 

  gsl_spline_init(f_spline, r_array, f_array, NSPLINE + 1);

  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {
    double m4 = compute_matrix_element_finite_q_op6(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, q, J1, J2, J, f_spline, acc); 
    mat += m4*density;
  }
  
  int phase = pow(-1.0, (1.0 + J1 - J2)/2.0);
  mat *= -phase/(4.0*sqrt(M_PI))*(2.0*J1 + 1.0)*clebsch_gordan(J1, 1.0, J2, 0.0, 0.0, 0.0)*clebsch_gordan(J1, 1.0, J, 0.0, 0.0, 0.0)*six_j(J1, 1, J2, 1, J, 1);    
 
  gsl_spline_free(f_spline);
  gsl_interp_accel_free(acc);  
                     
  return mat;
}

double compute_total_matrix_element_finite_q_op7(char* density_file, double q, int J1, int J2, int J) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;
 
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *f_spline = gsl_spline_alloc(gsl_interp_linear, NSPLINE + 1);
  double delta_r = (RMAX - RMIN)/(1.0*NSPLINE);
  double *f_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *r_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));

  for (int i = 0; i <= NSPLINE; i++) {
    double r = RMIN + i*delta_r;
    r_array[i] = r;
    f_array[i] = finite_q_alpha_pot_4(r, J1, q, M_PION);
  } 

  gsl_spline_init(f_spline, r_array, f_array, NSPLINE + 1);

  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {
    double m4 = compute_matrix_element_finite_q_op7(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, q, J1, J2, J, f_spline, acc); 
    mat += m4*density;
  }
                       
  int phase = pow(-1.0, (1.0 + J1 - J2)/2.0);
  mat *= -phase/(4.0*sqrt(M_PI))*(2.0*J1 + 1.0)*clebsch_gordan(J1, 1.0, J2, 0.0, 0.0, 0.0)*clebsch_gordan(J1, 1.0, J, 0.0, 0.0, 0.0)*six_j(J1, 1, J2, 2, J, 1);    

  return mat;
}


// Simplified expressions at q = 0 for cross-check

double compute_total_matrix_element_sigma_0(char* density_file) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;
 
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    double m4 = compute_matrix_element_sigma_0(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, 2); 
    mat += m4*density;
  }
  
  return mat;
}

double compute_matrix_element_sigma_0(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, int iv) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
    
    // The angular momentum are doubled in the file
  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double t12 = it12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;
  //  double t12p = it12p/2.0;

  int l1, l2, l1p, l2p;
     
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;

  double m4 = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    int s_max = MIN(lambda + j12, 1);
    int s_min = abs(lambda - j12);
    for (int s = s_min; s <= s_max; s++) {
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      fact *= sqrt(2*j12 + 1.0);
      double m1 = pow(-1.0, 1.0 + s)*six_j(s,0.5,0.5,1.0,0.5,0.5)*6.0;
      m1 *= pow(-1.0, 1.0 + t12)*sqrt(2.0*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*6.0;
      m1 *= fact;
      m1 *= compute_radial_matrix_element_scalar(iv, n1p, l1p, n2p, l2p, n1, l1, n2, l2, lambda, s, t12);
      m4 += m1;
    }
  }
  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= 1/sqrt(2);}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= 1/sqrt(2);}
                        
  return m4;
}

double compute_total_matrix_element_TT(char* density_file, int iv) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;
 
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    double m4 = compute_matrix_element_TT(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, iv); 
    mat += m4*density;
  }
                        
  return mat;
}


double compute_matrix_element_TT(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, int iv) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double t12 = it12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;

  double m4 = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  int lambda_min = MAX(abs(l1 - l2), abs(j12 - 1));
  int lambda_max = MIN(l1 + l2, j12 + 1);
  int lambdap_min = MAX(abs(l1p - l2p), abs(j12p - 1));
  int lambdap_max = MIN(l1p + l2p, j12p + 1);
  for (int lambda = lambda_min; lambda <= lambda_max; lambda++) {
    lambdap_min = MAX(lambdap_min, lambda - 2);
    lambdap_max = MIN(lambdap_max, lambda + 2);
    if (lambdap_min > lambdap_max) {continue;}
    for (int lambdap = lambdap_min; lambdap <= lambdap_max; lambdap++) {
      double fact = sqrt((2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, 1, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, 1, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      fact *= sqrt(2*lambda + 1)*sqrt(2*lambdap + 1);
      fact *= pow(-1.0, j12 + 1)*3.0;
      fact *= pow(-1.0, 1.0 + t12)*sqrt(2.0*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*6.0;
      fact *= six_j(j12, 1, lambdap, 2, lambda, 1);
      fact *= 2.0*sqrt(5.0)*sqrt(2*j12p + 1);
      fact *= compute_radial_matrix_element_y2(iv, n1p, l1p, n2p, l2p, lambdap, n1, l1, n2, l2, lambda, 1, t12);
      m4 += fact;
    }
  }
  if ((in1 == in2) && (j1 == j2)) {m4 *= 1/sqrt(2);}
  if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1/sqrt(2);}  
   m4 *= 1.0/(8.0*M_PI)*sqrt(8.0*M_PI/15.0);   

  return m4;
}

double compute_total_matrix_element_I2(char* density_file) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;
 
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    double m4 = compute_matrix_element_I2(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12); 
    mat += m4*density;
  }
                        
  return mat;
}


double compute_matrix_element_I2(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12) {
  // Computes matrix elements of two-body identity operator

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double t12 = it12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;
  if (j12 != j12p) {return 0.0;}

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;

  double m4 = 0.0;
  if ((in1 == in1p) && (j1 == j1p) && (in2 == in2p) && (j2 == j2p)) {
    m4 += 1.0;
  }
  if ((in1 == in2p) && (j1 == j2p) && (in2 == in1p) && (j2 == j1p)) {
    m4 += pow(-1.0, j1 + j2 + j12 + t12);
  }

  if ((in1 == in2) && (j1 == j2)) {m4 *= 1/sqrt(2);}
  if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1/sqrt(2);}  

  m4 *= sqrt(2.0*j12p + 1.0)*sqrt(2*t12 + 1.0);   

  return m4;
}


int get_l(int n, int j) {
  int l;
  if (n == 0) {
    l = 0;
   } else if (n == 1) {
     l = 1;
   } else if (n == 2) {
     if (j == 3 || j == 5) {
     l = 2;
     } else {
     l = 0;
     }
   }
  return l;
}

