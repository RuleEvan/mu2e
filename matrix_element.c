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
    if ((lambda < abs(l1p - l2p)) || (lambda > l1p + l2p)) {continue;} 
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
      double rm = compute_radial_matrix_element_finite_q_op1(J, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12, q, f_spline, acc);
      m4 += m1*rm;
    }
  }

  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= 1.0/sqrt(2.0);}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= 1.0/sqrt(2.0);}

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
      // Un-reduce wrt total J
      fact *= sqrt(2*j12 + 1);
      // Unreduced decoupling of Y2 dot S2 matrix element
      fact *= pow(-1.0, lambda + 1 + j12)*six_j(j12, 1, lambdap, 2, lambda, 1);
      // Reduced matrix element of tau(1) dot tau(2)
      fact *= pow(-1.0, 1.0 + t12)*sqrt(2.0*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*6.0;
      // Reduced matrix element of [sig(1) x sig(2)]_2
      fact *= 2.0*sqrt(5);
     
      if (fact == 0.0) {continue;}
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
      // Unreduce wrt total J
      fact *= sqrt(2*j12 + 1);
      // Unreduced decoupling of Y2 dot S2
      fact *= pow(-1.0, lambda + j12 + 1)*six_j(j12, 1, lambdap, 2, lambda, 1);
      // Reduced matrix element of S2
      fact *= 2.0*sqrt(5);
      // Reduced matrix element of tau(1) dot tau(2)
      fact *= pow(-1.0, 1.0 + t12)*sqrt(2.0*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*6.0;

      if (fact == 0.0) {continue;}
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
    if ((lambda < abs(l1p - l2p)) || (lambda > l1p + l2p)) {continue;} 
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
          m1 *= compute_radial_matrix_element_finite_q_op6(J1, J2, J, n1p, l1p, n2p, l2p, lambdap, n1, l1, n2, l2, lambda, s, t12, q, f_spline, acc);
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
  gsl_spline *f_spline = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  double delta_r = (RMAX - RMIN)/(1.0*NSPLINE);
  double *f_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *r_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));
//  printf("Beginning splines\n");
  for (int i = 0; i <= NSPLINE; i++) {
    double r = RMIN + i*delta_r;
    r_array[i] = r;
    f_array[i] = finite_q_alpha_pot_1(r, J, q, M_PION);
 //   printf("%g, %g\n", r_array[i], f_array[i]);
  }
//  exit(0);
  printf("Done splines\n");

  gsl_spline_init(f_spline, r_array, f_array, NSPLINE + 1);

  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;

  printf("Starting\n");
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
    double m4 = compute_matrix_element_finite_q_op1(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, q, J, f_spline, acc); 
    mat += m4*density;

  }

  mat *= -1.0/(12.0*sqrt(M_PI));
  gsl_spline_free(f_spline);
  gsl_interp_accel_free(acc);  
  
  return mat;
}

void convert_density_file(char *density_file) {

  FILE *in_file;
  FILE *out_file;
  in_file = fopen(density_file, "r");
  out_file = fopen("Al27_exsp_converted_J0_T0.dens", "w");
  double mat = 0.0;
 
  int ind1, ind2, ind1p, ind2p;
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%f\n", &ind1p, &ind2p, &ind1, &ind2, &ij12p, &it12p, &it12, &density) == 8) {
    in1 = get_n(ind1);
    ij1 = get_j(ind1);
    in2 = get_n(ind2);
    ij2 = get_j(ind2);
    in1p = get_n(ind1p);
    ij1p = get_j(ind1p);
    in2p = get_n(ind2p);
    ij2p = get_j(ind2p);
    ij12 = ij12p;
    fprintf(out_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", in1p, ij1p, in2p, ij2p, 2*ij12p, 2*it12p, in1, ij1, in2, ij2, 2*ij12, 2*it12, density); 

  }
  fclose(in_file);
  fclose(out_file);
  return;
}

double compute_total_matrix_element_finite_q_op1_alt(char* density_file, double q, int J) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;
 
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *f_spline = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
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

  int ind1, ind2, ind1p, ind2p;
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%f\n", &ind1p, &ind2p, &ind1, &ind2, &ij12p, &it12p, &it12, &density) == 8) {
    in1 = get_n(ind1);
    ij1 = get_j(ind1);
    in2 = get_n(ind2);
    ij2 = get_j(ind2);
    in1p = get_n(ind1p);
    ij1p = get_j(ind1p);
    in2p = get_n(ind2p);
    ij2p = get_j(ind2p);
    ij12 = ij12p;


    double m4 = compute_matrix_element_finite_q_op1(in1p, ij1p, in2p, ij2p, 2*ij12p, in1, ij1, in2, ij2, 2*ij12, 2*it12, q, J, f_spline, acc); 
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
  gsl_spline *f_spline = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
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
    if (fabs(density) < pow(10, -8)) {continue;}

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
  gsl_spline *f_spline = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
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
    if (fabs(density) < pow(10, -8)) {continue;}

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
  gsl_spline *f_spline = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
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
    if (fabs(density) < pow(10, -8)) {continue;}

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
  gsl_spline *f_spline = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
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
    if (fabs(density) < pow(10, -8)) {continue;}
    
    double m4 = compute_matrix_element_finite_q_op6(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, q, J1, J2, J, f_spline, acc);
    if (fabs(m4) < pow(10, -8)) {continue;}
    mat += m4*density;
  }
  
  double phase = pow(-1.0, (1.0 + J1 - J2)/2.0);
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
  gsl_spline *f_spline = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
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
    if (fabs(density) < pow(10, -8)) {continue;}

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
  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= sqrt(1.0 - pow(-1.0, j12 + t12))/2.0;}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= sqrt(1.0 - pow(-1.0, j12 + t12))/2.0;}
                        
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

double compute_total_matrix_element_I2_cj(char* density_file) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;
 
  int ind1, ind2, ind1p, ind2p;
  int j12, j12p, jmin, jmax;
  float density;
 
  while(fscanf(in_file, "    %d   %d   %d   %d     %d     %f\n", &ind1p, &ind2p, &ind1, &ind2, &j12, &density) == 6) {
     double m = 0.0;
     double j1p, j1, j2p, j2;
     double mt1, mt2;
     if (ind1 == 1 || ind1 == 4) {j1 = 1.5;}
     if (ind1 == 2 || ind1 == 5) {j1 = 2.5;}
     if (ind1 == 3 || ind1 == 6) {j1 = 0.5;}
     if (ind2 == 1 || ind2 == 4) {j2 = 1.5;}
     if (ind2 == 2 || ind2 == 5) {j2 = 2.5;}
     if (ind2 == 3 || ind2 == 6) {j2 = 0.5;}
     if (ind1 <= 3) {mt1 = 0.5;}
     else {mt1 = -0.5;}
     if (ind2 <= 3) {mt2 = 0.5;}
     else {mt2 = -0.5;}


     if ((ind1p == ind1) && (ind2p == ind2)) {
       m = 1.0;
     }
  //   if ((ind1p == ind2) && (ind2p == ind1)) {
  //     m += pow(-1.0, 1 + j1 + j2 - j12);
  //   }
     if (ind1 == ind2) {m *= 1.0/sqrt(2.0);}
     if (ind1p == ind2p) {m *= 1.0/sqrt(2.0);}
     m *= sqrt(2*j12 + 1); 
     mat += m*density;
//    double m4 = compute_matrix_element_I2(in1p, ij1p, in2p, ij2p, 2*ij12p, in1, ij1, in2, ij2, 2*ij12, 2*it12); 
   // mat += m4*density;
  }
                        
  return mat;
}


void translate_cj_to_er(char* density_file) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;
 
  int ind1, ind2, ind1p, ind2p;
  int j12, j12p, jmin, jmax;
  float density;
  double *translated = (double*) calloc(3*3*3*3*6*2, sizeof(double));
  double *pn_mat = (double*) calloc(6*6*6*6*6, sizeof(double));

  while(fscanf(in_file, "    %d   %d   %d   %d     %d     %f\n", &ind1p, &ind2p, &ind1, &ind2, &j12, &density) == 6) {
     int ij1p, ij1, ij2p, ij2;
     double j1, j2, j1p, j2p;
     
     if (ind1 == 1 || ind1 == 4) {ij1 = 3; j1 = 1.5;}
     else if (ind1 == 2 || ind1 == 5) {ij1 = 5; j1 = 2.5;}
     else if (ind1 == 3 || ind1 == 6) {ij1 = 1; j1 = 0.5;}
    
     if (ind2 == 1 || ind2 == 4) {ij2 = 3; j2 = 1.5;}
     else if (ind2 == 2 || ind2 == 5) {ij2 = 5; j2 = 2.5;}
     else if (ind2 == 3 || ind2 == 6) {ij2 = 1; j2 = 0.5;}
    
     if (ind1p == 1 || ind1p == 4) {ij1p = 3; j1p = 1.5;}
     else if (ind1p == 2 || ind1p == 5) {ij1p = 5; j1p = 2.5;}
     else if (ind1p == 3 || ind1p == 6) {ij1p = 1; j1p = 0.5;}
     
     if (ind2p == 1 || ind2p == 4) {ij2p = 3; j2p = 1.5;}
     else if (ind2p == 2 || ind2p == 5) {ij2p = 5; j2p = 2.5;}
     else if (ind2p == 3 || ind2p == 6) {ij2p = 1; j2p = 0.5;}
     
     // Inital matrix element
     int index = (ind1p - 1) + 6*((ind2p - 1) + 6*((ind1 - 1) + 6*((ind2 - 1) + 6*j12)));
     pn_mat[index] = density;
     // swap a <-> b
     index = (ind2p - 1) + 6*((ind1p - 1) + 6*((ind1 - 1) + 6*((ind2 - 1) + 6*j12)));
     pn_mat[index] = density*pow(-1.0, 1 + j1p + j2p - j12);
     // swap c <-> d
     index = (ind1p - 1) + 6*((ind2p - 1) + 6*((ind2 - 1) + 6*((ind1 - 1) + 6*j12)));
     pn_mat[index] = density*pow(-1.0, 1 + j1 + j2 - j12);
     // swap a <-> b and c <-> d
     index = (ind2p - 1) + 6*((ind1p - 1) + 6*((ind2 - 1) + 6*((ind1 - 1) + 6*j12)));
     pn_mat[index] = density*pow(-1.0, j1 + j2 + j1p + j2p);
     // cd ab
     index = (ind1 - 1) + 6*((ind2 - 1) + 6*((ind1p - 1) + 6*((ind2p - 1) + 6*j12)));
     pn_mat[index] = density;
     // cd ba
     index = (ind1 - 1) + 6*((ind2 - 1) + 6*((ind2p - 1) + 6*((ind1p - 1) + 6*j12)));
     pn_mat[index] = density*pow(-1.0, j1p + j2p + 1 + j12);
     // dc ab    
     index = (ind2 - 1) + 6*((ind1 - 1) + 6*((ind1p - 1) + 6*((ind2p - 1) + 6*j12)));
     pn_mat[index] = density*pow(-1.0, 1 + j12 + j1 + j2);
     // dc ba
     index = (ind2 - 1) + 6*((ind1 - 1) + 6*((ind2p - 1) + 6*((ind1p - 1) + 6*j12)));
     pn_mat[index] = density*pow(-1.0, j1 + j2 + j1p + j2p);
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_no_core_2body_J0_T0_0_0.dens", "r");
  
  int in1p, in2p, in1, in2, ij1, ij2, ij12, ij1p, ij2p, ij12p, it12, it12p;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {
    if (ij1p == 3) {ind1p = 1;}
    else if( ij1p == 5) {ind1p = 2;}
    else if (ij1p == 1) {ind1p = 3;}
    if (ij2p == 3) {ind2p = 1;}
    else if( ij2p == 5) {ind2p = 2;}
    else if (ij2p == 1) {ind2p = 3;}
    if (ij1 == 3) {ind1 = 1;}
    else if( ij1 == 5) {ind1 = 2;}
    else if (ij1 == 1) {ind1 = 3;}
    if (ij2 == 3) {ind2 = 1;}
    else if( ij2 == 5) {ind2 = 2;}
    else if (ij2 == 1) {ind2 = 3;}
    
    j12 = ij12/2;
    int t12 = it12/2;

    double result = 0.0;
    double norm = 1.0;
    double normp = 1.0;
    if (ind1p == ind2p) {normp *= 1.0/sqrt(2.0);}
    if (ind1 == ind2) {norm *= 1.0/sqrt(2.0);}

    if (t12 == 0) {
      int index  = (ind1p - 1) + 6*((ind2p + 3 - 1) + 6*((ind1 - 1) + 6*((ind2 + 3 - 1) + 6*j12)));
      result += pn_mat[index];
      index  = (ind1p  + 3 - 1) + 6*((ind2p  - 1) + 6*((ind1 - 1) + 6*((ind2 + 3 - 1) + 6*j12)));
      result -= pn_mat[index];
      index  = (ind1p - 1) + 6*((ind2p + 3 - 1) + 6*((ind1 + 3 - 1) + 6*((ind2 - 1) + 6*j12)));
      result -= pn_mat[index];
      index  = (ind1p + 3 - 1) + 6*((ind2p - 1) + 6*((ind1 + 3 - 1) + 6*((ind2 - 1) + 6*j12)));
      result += pn_mat[index];
      result *= sqrt(2.0)/2.0*norm*normp;
    } else if (t12 == 1) {
      int index  = (ind1p - 1) + 6*((ind2p + 3 - 1) + 6*((ind1 - 1) + 6*((ind2 + 3 - 1) + 6*j12)));
      result += 0.5*pn_mat[index]*norm*normp;
      index  = (ind1p  + 3 - 1) + 6*((ind2p  - 1) + 6*((ind1 - 1) + 6*((ind2 + 3 - 1) + 6*j12)));
      result += 0.5*pn_mat[index]*norm*normp;
      index  = (ind1p - 1) + 6*((ind2p + 3 - 1) + 6*((ind1 + 3 - 1) + 6*((ind2 - 1) + 6*j12)));
      result += 0.5*pn_mat[index]*norm*normp;
      index  = (ind1p + 3 - 1) + 6*((ind2p - 1) + 6*((ind1 + 3 - 1) + 6*((ind2 - 1) + 6*j12)));
      result += 0.5*pn_mat[index]*norm*normp;
      index  = (ind1p - 1) + 6*((ind2p - 1) + 6*((ind1 - 1) + 6*((ind2 - 1) + 6*j12)));
      result += pn_mat[index];
      index  = (ind1p + 3 - 1) + 6*((ind2p + 3 - 1) + 6*((ind1 + 3 - 1) + 6*((ind2 + 3 - 1) + 6*j12)));
      result += pn_mat[index];
      result *= sqrt(2.0)/sqrt(3.0);
    }

   //  if (ind1p == ind2p) {result *= 1.0/sqrt(2.0);}
     //if (ind1 == ind2) {result *= 1.0/sqrt(2.0);}

    if (fabs(result/density - 1.0) < pow(10, -4)) {printf("    %d   %d   %d   %d     %d      %d      % 1.7f      % 1.7f\n", ind1p, ind2p, ind1, ind2, j12, t12, result, density);}
    else {printf("Failure! %d %d %d %d %d %d %g %g\n", ind1p, ind2p, ind1, ind2, j12, t12, result, density);}
  }

  fclose(in_file);
/*
     double m = 0.0;
          double mt1, mt1p, mt2, mt2p;
     
     if (ind1 == 1 || ind1 == 4) {ij1 = 3; j1 = 1.5;}
     else if (ind1 == 2 || ind1 == 5) {ij1 = 5; j1 = 2.5;}
     else if (ind1 == 3 || ind1 == 6) {ij1 = 1; j1 = 0.5;}
    
     if (ind2 == 1 || ind2 == 4) {ij2 = 3; j2 = 1.5;}
     else if (ind2 == 2 || ind2 == 5) {ij2 = 5; j2 = 2.5;}
     else if (ind2 == 3 || ind2 == 6) {ij2 = 1; j2 = 0.5;}
    
     if (ind1p == 1 || ind1p == 4) {ij1p = 3; j1p = 1.5;}
     else if (ind1p == 2 || ind1p == 5) {ij1p = 5; j1p = 2.5;}
     else if (ind1p == 3 || ind1p == 6) {ij1p = 1; j1p = 0.5;}
     
     if (ind2p == 1 || ind2p == 4) {ij2p = 3; j2p = 1.5;}
     else if (ind2p == 2 || ind2p == 5) {ij2p = 5; j2p = 2.5;}
     else if (ind2p == 3 || ind2p == 6) {ij2p = 1; j2p = 0.5;}
     
     if (ind1p > 3 && ind2p > 3 && ind1 > 3 && ind2 > 3) {
       int t12 = 1;
       double fact = 1.0;
    //   if (ind1p != ind1 || ind2p != ind2) {fact = 2.0;}
       int index = (ind1p - 4) + 3*((ind2p - 4) + 3*((ind1 - 4) + 3*((ind2 - 4) + 3*(j12 + 6*t12))));
       translated[index] += sqrt(2.0)/sqrt(3.0)*density*fact;

     } else if (ind1p <=3 && ind2p <=3 && ind1 <= 3 && ind2 <= 3) {
       int t12 = 1;
       double fact = 1.0;
     //  if (ind1p != ind1 || ind2p != ind2) {fact = 2.0;}
       int index = (ind1p - 1) + 3*((ind2p - 1) + 3*((ind1 - 1) + 3*((ind2 - 1) + 3*(j12 + 6*t12))));
       translated[index] += sqrt(2.0)/sqrt(3.0)*density*fact;

     } else if (ind1p <= 3 && ind2p > 3 && ind1 <= 3 && ind2 > 3) {
       for (int t12 = 0; t12 <= 1; t12++) {
         double cg1 = 0.5;
         double cg3 = clebsch_gordan(t12, t12, 0, 0, 0, 0);
         if (ind1p != ind1 || ind2p != ind2) {cg3 *= 2.0;}
         if (ind1 == ind2 - 3) {cg3 *= 1.0/sqrt(2.0);}
         if (ind1p == ind2p - 3) {cg3 *= 1.0/sqrt(2.0);}
	 if ((ind1p == ind2p - 3) && (ind1 == ind2 - 3)) {
           int index = (ind1p - 1) + 3*((ind2p - 4) + 3*((ind1 - 1) + 3*((ind2 - 4) + 3*(j12 + 6*t12))));
	   double sym = (1.0 + pow(-1.0, j1 + j2 + j12 + t12))*(1.0 + pow(-1.0, j1p + j2p + j12 + t12));
	   translated[index] += sqrt(2.0)*cg1*cg3*sym*density; 

	 } else if ((ind1p > ind2p - 3) && (ind1 == ind2 - 3)) {
           int index = (ind1p - 1) + 3*((ind2p - 4) + 3*((ind1 - 1) + 3*((ind2 - 4) + 3*(j12 + 6*t12))));
           double sym = 1.0 + pow(-1.0, j1 + j2 + j12 + t12);
	   translated[index] += sqrt(2.0)*cg1*cg3*sym*density; 

         } else if ((ind1p < ind2p - 3) && (ind1 == ind2 - 3)) {
           int index = (ind2p - 4) + 3*((ind1p - 1) + 3*((ind1 - 1) + 3*((ind2 - 4) + 3*(j12 + 6*t12))));
           double sym = 1.0 + pow(-1.0, j1 + j2 + j12 + t12);
	   translated[index] += sqrt(2.0)*cg1*cg3*pow(-1.0,  j1p + j2p - j12 - t12)*sym*density; 

         } else if ((ind1p == ind2p - 3) && (ind1 > ind2 - 3)) {
           int index = (ind1p - 1) + 3*((ind2p - 4) + 3*((ind1 - 1) + 3*((ind2 - 4) + 3*(j12 + 6*t12))));
           double sym = 1.0 + pow(-1.0,  j1p + j2p + j12 + t12);
	   translated[index] += sqrt(2.0)*cg1*cg3*sym*density; 

         } else if ((ind1p == ind2p - 3) && (ind1 < ind2 - 3)) {
           int index = (ind1p - 1) + 3*((ind2p - 4) + 3*((ind2 - 4) + 3*((ind1 - 1) + 3*(j12 + 6*t12))));
           double sym = 1.0 + pow(-1.0,  j1p + j2p + j12 + t12);
	   translated[index] += sqrt(2.0)*cg1*cg3*sym*pow(-1.0,  j1 + j2 - j12 - t12)*density; 

         } else if ((ind1p > ind2p - 3) && (ind1 > ind2 - 3)) {
           int index = (ind1p - 1) + 3*((ind2p - 4) + 3*((ind1 - 1) + 3*((ind2 - 4) + 3*(j12 + 6*t12))));
	   translated[index] += sqrt(2.0)*cg1*cg3*density; 

         } else if ((ind1p > ind2p - 3) && (ind1 < ind2 - 3)) {
           int index = (ind1p - 1) + 3*((ind2p - 4) + 3*((ind2 - 4) + 3*((ind1 - 1) + 3*(j12 + 6*t12))));
	   translated[index] += sqrt(2.0)*cg1*cg3*pow(-1.0,  j1 + j2 - j12 - t12)*density; 

         } else if ((ind1p < ind2p - 3) && (ind1 > ind2 - 3)) {
           int index = (ind2p - 4) + 3*((ind1p - 1) + 3*((ind1 - 1) + 3*((ind2 - 4) + 3*(j12 + 6*t12))));
	   translated[index] += sqrt(2.0)*cg1*cg3*pow(-1.0,  j1p + j2p - j12 - t12)*density; 

         } else if ((ind1p < ind2p - 3) && (ind1 < ind2 - 3)) {
           int index = (ind2p - 4) + 3*((ind1p - 1) + 3*((ind2 - 4) + 3*((ind1 - 1) + 3*(j12 + 6*t12))));
	   translated[index] += sqrt(2.0)*cg1*cg3*pow(-1.0, j1 + j2 + j1p + j2p)*density; 

         }
      }
    }

  }
  fclose(in_file);  
  int *j_shell = (int*) malloc(sizeof(int)*3);
  j_shell[0] = 3;
  j_shell[1] = 5;
  j_shell[2] = 1;

  in_file = fopen("isotope_data/al27/density/al27-al27_no_core_2body_J0_T0_0_0.dens", "r");
  
  double* my_dens = (double*) calloc(3*3*3*3*6*2, sizeof(double));


  int in1p, in2p, in1, in2, ij1, ij2, ij12, ij1p, ij2p, ij12p, it12, it12p;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {
    j12 = ij12/2.0;
    int t12 = it12/2;
    if (ij1p == 3) {ind1p = 1;}
    else if( ij1p == 5) {ind1p = 2;}
    else if (ij1p == 1) {ind1p = 3;}
    if (ij2p == 3) {ind2p = 1;}
    else if( ij2p == 5) {ind2p = 2;}
    else if (ij2p == 1) {ind2p = 3;}
    if (ij1 == 3) {ind1 = 1;}
    else if( ij1 == 5) {ind1 = 2;}
    else if (ij1 == 1) {ind1 = 3;}
    if (ij2 == 3) {ind2 = 1;}
    else if( ij2 == 5) {ind2 = 2;}
    else if (ij2 == 1) {ind2 = 3;}
  
    if (ind1p < ind2p && ind1 < ind2) {
      int index = (ind2p - 1) + 3*((ind1p - 1) + 3*((ind2 - 1) + 3*((ind1 - 1) + 3*(j12 + 6*t12))));
      my_dens[index] += density*pow(-1.0, (ij1 + ij2 + ij1p + ij2p)/2.0);
    } else if (ind1p < ind2p && ind1 >= ind2) {
      int index = (ind2p - 1) + 3*((ind1p - 1) + 3*((ind1 - 1) + 3*((ind2 - 1) + 3*(j12 + 6*t12))));
      my_dens[index] += density*pow(-1.0, (ij1p + ij2p)/2.0 + j12 + t12);
    } else if (ind1p >= ind2p && ind1 < ind2) {
      int index = (ind1p - 1) + 3*((ind2p - 1) + 3*((ind2 - 1) + 3*((ind1 - 1) + 3*(j12 + 6*t12))));
      my_dens[index] += density*pow(-1.0, (ij1 + ij2)/2.0 + j12 + t12);
    } else if (ind1p >= ind2p && ind1 >= ind2) {
      int index = (ind1p - 1) + 3*((ind2p - 1) + 3*((ind1 - 1) + 3*((ind2 - 1) + 3*(j12 + 6*t12))));
      my_dens[index] += density;
    }
  }

  int num_right = 0;
  for (ind1p = 1; ind1p <= 3; ind1p++) {
    int ij1p = j_shell[ind1p - 1];
    double j1p = ij1p/2.0;
    for (ind2p = 1; ind2p <= ind1p; ind2p++) {
     int ij2p = j_shell[ind2p - 1];
     double j2p = ij2p/2.0; 
     for (ind1 = 1; ind1 <= 3; ind1++) {
       int ij1 = j_shell[ind1 - 1];
       double j1 = ij1/2.0;
	for (ind2 = 1; ind2 <= ind1; ind2++) {
          int ij2 = j_shell[ind2 - 1];
	  double j2 = ij2/2.0;
	  for (int t12 = 0; t12 <= 1; t12++) {
	    for (int j12 = 0; j12 < 6; j12++) {
	      int index = (ind1p - 1) + 3*((ind2p - 1) + 3*((ind1 - 1) + 3*((ind2 - 1) + 3*(j12 + 6*t12))));
	      if ((fabs(translated[index]) < pow(10, -8)) || fabs(my_dens[index]) < pow(10, -8)) {continue;}
	      //printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%g\n", 2, ij1p, 2, ij2p, 2*j12, 2*t12, 2, ij1, 2, ij2, 2*j12, 2*t12, translated[index]);
	      if (fabs(translated[index] - my_dens[index]) > pow(10, -3)) {printf("error: %d %d %d %d %d %d %g %g\n", ind1p, ind2p, ind1, ind2, j12, t12, translated[index], my_dens[index]);} else {
		     printf("success: %d %d %d %d %d %d %g %g\n", ind1p, ind2p, ind1, ind2, j12, t12, translated[index], my_dens[index]); num_right++;}

     	    }
	  }
	}
      }
    }
  }
  printf("num right: %d\n", num_right);
*/
  return;
}


double compute_total_matrix_element_I2_alt(char* density_file) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;
 
  int ind1, ind2, ind1p, ind2p;
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%f\n", &ind1p, &ind2p, &ind1, &ind2, &ij12p, &it12p, &it12, &density) == 8) {
    in1 = get_n(ind1);
    ij1 = get_j(ind1);
    in2 = get_n(ind2);
    ij2 = get_j(ind2);
    in1p = get_n(ind1p);
    ij1p = get_j(ind1p);
    in2p = get_n(ind2p);
    ij2p = get_j(ind2p);
    ij12 = ij12p;
    
    double m4 = compute_matrix_element_I2(in1p, ij1p, in2p, ij2p, 2*ij12p, in1, ij1, in2, ij2, 2*ij12, 2*it12); 
    mat += m4*density;
  }
                        
  return mat;
}

int get_j(int index) {
  int j;
  if (index == 1 || index == 2 || index == 4) {j = 1;} 
  if (index == 3 || index == 5) {j = 3;}
  if (index == 6) {j = 5;}

   return j;
}

int get_n(int index) {
  int n;
  if (index == 1) {n = 0;}
  if (index == 2 || index == 3) {n = 1;}
  if (index == 4 || index == 5 || index == 6) {n = 2;}

  return n;
}

double compute_total_relative_density(char* density_file, double r) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;
 
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    double m4 = compute_relative_density(r, in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12); 
    mat += m4*density;
  }
                        
  return mat;
}

double compute_relative_density(double r, int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12) {
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
      m1 *= compute_relative_density_scalar(r, n1p, l1p, n2p, l2p, n1, l1, n2, l2, lambda, s, t12);
      m4 += m1;
    }
  }
  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= 1/sqrt(2);}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= 1/sqrt(2);}
                        
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
  int ind1p = get_shell_index(in1p, ij1p); 
  int ind2p = get_shell_index(in2p, ij2p); 
  int ind1 = get_shell_index(in1, ij1); 
  int ind2 = get_shell_index(in2, ij2); 
//  if (ind1p < ind2p) {return 0.0;}
//  if (ind1 < ind2) {return 0.0;}
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

  if ((in1 == in2) && (j1 == j2)) {m4 *= 1.0/sqrt(2.0);}
  if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1.0/sqrt(2.0);}  

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
   } else if (n == 3) {
     if (j == 1 || j == 3) {
       l = 1;
     } else {
       l = 3;
     }
   } else if (n == 4) {
     if (j == 9 || j == 7) {
       l = 4;
     }
   } else {
    printf("get_l error: Model space too large!\n"); exit(0);}
  return l;
}

void generate_bigstick_int_file() {
  int *n_shell = (int*) malloc(sizeof(int)*11);
  int *j_shell = (int*) malloc(sizeof(int)*11);
  int *l_shell = (int*) malloc(sizeof(int)*11);
  int *i_core = (int*) malloc(sizeof(int)*11);
  double q = 103.583;

  n_shell[0] = 3;
  j_shell[0] = 3;
  l_shell[0] = 1;
  i_core[0] = 0;

  n_shell[1] = 3;
  j_shell[1] = 5;
  l_shell[1] = 3;
  i_core[1] = 0;

  n_shell[2] = 3;
  j_shell[2] = 1;
  l_shell[2] = 1;
  i_core[2] = 0;

  n_shell[3] = 4;
  j_shell[3] = 9;
  l_shell[3] = 4;
  i_core[3] = 0;

  n_shell[4] = 3;
  j_shell[4] = 7;
  l_shell[4] = 3;
  i_core[4] = 1;

  n_shell[5] = 2;
  j_shell[5] = 3;
  l_shell[5] = 2;
  i_core[5] = 1;

  n_shell[6] = 2;
  j_shell[6] = 1;
  l_shell[6] = 0;
  i_core[6] = 1;

  n_shell[7] = 2;
  j_shell[7] = 5;
  l_shell[7] = 2;
   i_core[7] = 1;

  n_shell[8] = 1;
  j_shell[8] = 1;
  l_shell[8] = 1;
  i_core[8] = 1;

  n_shell[9] = 1;
  j_shell[9] = 3;
  l_shell[9] = 1;
  i_core[9] = 1;

  n_shell[10] = 0;
  j_shell[10] = 1;
  l_shell[10] = 0;
  i_core[10] = 1;

  /*  
  n_shell[0] = 3;
  j_shell[0] = 7;
  l_shell[0] = 3;
  i_core[0] = 0;

  n_shell[1] = 3;
  j_shell[1] = 3;
  l_shell[1] = 1;
  i_core[1] = 0;

  n_shell[2] = 3;
  j_shell[2] = 5;
  l_shell[2] = 3;
  i_core[2] = 0;

  n_shell[3] = 3;
  j_shell[3] = 1;
  l_shell[3] = 1;
  i_core[3] = 0;

  n_shell[4] = 2;
  j_shell[4] = 3;
  l_shell[4] = 2;
  i_core[4] = 1;

  n_shell[5] = 2;
  j_shell[5] = 1;
  l_shell[5] = 0;
  i_core[5] = 1;

  n_shell[6] = 2;
  j_shell[6] = 5;
  l_shell[6] = 2;
   i_core[6] = 1;

  n_shell[7] = 1;
  j_shell[7] = 1;
  l_shell[7] = 1;
  i_core[7] = 1;

  n_shell[8] = 1;
  j_shell[8] = 3;
  l_shell[8] = 1;
  i_core[8] = 1;

  n_shell[9] = 0;
  j_shell[9] = 1;
  l_shell[9] = 0;
  i_core[9] = 1;
*/
  FILE* out_file;
  out_file = fopen("ni58_tbme.int", "w");

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *f1_spline = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  gsl_spline *f4_spline = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  gsl_spline *f5_spline = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  gsl_spline *f6_spline = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);

  double delta_r = (RMAX - RMIN)/(1.0*NSPLINE);
  double *f1_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *f4_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *f5_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *f6_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));

  double *r_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));

  for (int i = 0; i <= NSPLINE; i++) {
    double r = RMIN + i*delta_r;
    r_array[i] = r;
    f1_array[i] = finite_q_alpha_pot_1(r, 0, q, M_PION);
    f4_array[i] = finite_q_alpha_pot_3(r, 0, q, M_PION);
    f5_array[i] = finite_q_alpha_pot_4(r, 1, q, M_PION);
    f6_array[i] = finite_q_alpha_pot_5(r, 0, q, M_PION);

  }

  gsl_spline_init(f1_spline, r_array, f1_array, NSPLINE + 1);
  gsl_spline_init(f4_spline, r_array, f4_array, NSPLINE + 1);
  gsl_spline_init(f5_spline, r_array, f5_array, NSPLINE + 1);
  gsl_spline_init(f6_spline, r_array, f6_array, NSPLINE + 1);


  for (int ia = 0; ia < 11; ia++) {
    int n1p = n_shell[ia];
    int j1p = j_shell[ia];
    int l1p = l_shell[ia];
    int ic1p = i_core[ia];
    for (int ib = 0; ib <= ia; ib++) {
      int n2p = n_shell[ib];
      int j2p = j_shell[ib];
      int l2p = l_shell[ib];
      int ic2p = i_core[ib];
      int j12p_min = abs(j1p - j2p)/2;
      int j12p_max = (j1p + j2p)/2;
      for (int ig = 0; ig < 11; ig++) {
	int n1 = n_shell[ig];
	int j1 = j_shell[ig];
	int l1 = l_shell[ig];
	int ic1 = i_core[ig];
	if ((ic1p == 1) && (ic2p == 1)) {
          if (ic1 == 0) {continue;}
	  if ((ig != ia) && (ig != ib)) {continue;}
	}
        for (int id = 0; id <= ig; id++) {
	  int n2 = n_shell[id];
	  int j2 = j_shell[id];
	  int l2 = l_shell[id];
	  int ic2 = i_core[id];
	  if (ic1p + ic2p != ic1 + ic2) {continue;}
	  if (ic1p == 1) {
            if ((ia != ig) && (ia != id)) {continue;}
	  }
	  if (ic2p == 1) {
	    if ((ib != ig) && (ib != id)) {continue;}
	  }
	  if (pow(-1, l1 + l2) != pow(-1, l1p + l2p)) {continue;}
	  int j12_min = abs(j1 - j2)/2;
	  int j12_max = (j1 + j2)/2;
	  for (int j12 = MAX(j12_min, j12p_min); j12 <= MIN(j12_max, j12p_max); j12++) {
            for (int t12 = 0; t12 <= 1; t12++) {
	    if ((ia == ib) && ((j12 + t12) % 2 == 0)) {continue;}
	    if ((ig == id) && ((j12 + t12) % 2 == 0)) {continue;}
	    printf("%d %d %d %d\n", ia, ib, ig, id);
	    //  if (pow(-1.0, (j1 + j2)/2 + j12 + 1 + t12) != -1.0) {continue;}
            //  if (pow(-1.0, (j1p + j2p)/2 + j12 + 1 + t12) != -1.0) {continue;}
	      double op1 =  -1.0/(12.0*sqrt(M_PI))*compute_matrix_element_finite_q_op1(n1p, j1p, n2p, j2p, 2*j12, n1, j1, n2, j2, 2*j12, 2*t12, q, 0, f1_spline, acc); 
              double op4 = clebsch_gordan(0, 2, 2, 0.0, 0.0, 0.0)*sqrt(8.0*M_PI/15.0)/(8.0*M_PI)*compute_matrix_element_finite_q_op4(n1p, j1p, n2p, j2p, 2*j12, n1, j1, n2, j2, 2*j12, 2*t12, q, 2, 0, f4_spline, acc);
              double op5 = -1.0/(12.0*sqrt(M_PI))*(2.0 + 1.0)*pow(clebsch_gordan(1.0, 1.0, 0.0, 0.0, 0.0, 0.0), 2.0)*compute_matrix_element_finite_q_op5(n1p, j1p, n2p, j2p, 2*j12, n1, j1, n2, j2, 2*j12, 2*t12, q, 1, 0, f5_spline, acc);
              double op6 = 1.0/(12.0*sqrt(M_PI))*compute_matrix_element_finite_q_op6(n1p, j1p, n2p, j2p, 2*j12, n1, j1, n2, j2, 2*j12, 2*t12, q, 0, 1, 1, f6_spline, acc);

	      double mat = 2.0*(op1 + op4 + op5 + op6);
              if (fabs(mat) > pow(10, -8)) {
		mat /= sqrt((2*j12 + 1)*(2*t12 + 1));
	 	if (ia != ig || ib != id) {mat *= 0.5;}
	    //    if (ib != id) {mat *= 0.5;}
	        fprintf(out_file, "  %d  %d  %d  %d    %d  %d  %g\n", ia + 1, ib + 1, ig + 1, id + 1, j12, t12, mat);
              }
	    }
	  }
	}
      }
    }
  }

  fclose(out_file);

  return;
}
          
