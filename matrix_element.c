#include "matrix_element.h"

/* One-body matrix elements */

/* Two-body matrix elements */

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


double compute_matrix_element_tau_plus(char* density_file, int iv) {
  double mat = 0.0;
  FILE *in_file;
  in_file = fopen(density_file, "r");

  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;

  // Each line of the file corresponds to a nuclear shell
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {
    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j2 = ij2/2.0;
    double j12 = ij12/2.0;
    double t12 = it12/2.0;
    double j1p = ij1p/2.0;
    double j2p = ij2p/2.0;
    double j12p = ij12p/2.0;
    double t12p = it12p/2.0;
    if (j12 != j12p) {continue;}
    int lambda, s;
    int l1, l2, l1p, l2p;
    
    if (j1 == 4.5) {l1 = 4;}
    else if (j1 == 2.5) {l1 = 3;}
    else {l1 = 1;}
    if (j2 == 4.5) {l2 = 4;}
    else if (j2 == 2.5) {l2 = 3;}
    else {l2 = 1;}
    if (j1p == 4.5) {l1p = 4;}
    else if (j1p == 2.5) {l1p = 3;}
    else {l1p = 1;}
    if (j2p == 4.5) {l2p = 4;}
    else if (j2p == 2.5) {l2p = 3;}
    else {l2p = 1;}

    
    // The N's listed in the input file are energy quanta, we want radial quantum numbers
    double n1 = (in1)/2.0;
    double n2 = (in2)/2.0;
    double n1p = (in1p)/2.0;
    double n2p = (in2p)/2.0;
    double m4 = 0.0;
    // Convert from JJ to LS coupling (L is lambda)
    int lambda_min = MAX(abs(l1 - l2), abs(l1p - l2p));
    int lambda_max = MIN(l1 + l2, l1p + l2p);
    for (lambda = lambda_min; lambda <= lambda_max; lambda++) {
      int s_max = MIN(lambda + j12, 1);
      int s_min = abs(lambda - j12);
      for (s = s_min; s <= s_max; s++) {
        double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
        fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
        fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
        fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
        if (fact == 0.0) {continue;}
        fact *= sqrt(2*j12 + 1.0);
        double m1 = sqrt(5.0);
        m1 *= fact;
        if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m1 *= 1/sqrt(2);}
        if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m1 *= 1/sqrt(2);}  

        m1 *= compute_radial_matrix_element_scalar(iv, n1p, l1p, n2p, l2p, n1, l1, n2, l2, lambda, s, t12);
        m4 += m1;
      }
    }
    
    mat += density*m4;
  }
                        
  return mat;
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

double compute_matrix_element_sigma_0_finite_q_op5(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, double q, int J1, int J2, int it12) {
  int phase = pow(-1.0, (1.0 + J1 - J2)/2.0);

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
      m1 *= pow(-1.0, 1.0 + t12)*sqrt(2*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*6.0;
      m1 *= fact;
      m1 *= compute_radial_matrix_element_J_dot_J_alt(J1, J2, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12, q);
      m4 += m1;
    }
  }
  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= 1/sqrt(2);}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= 1/sqrt(2);}
  m4 *= phase/(12.0*sqrt(M_PI))*(2.0*J1 + 1.0)/(2.0*J2 + 1.0)*pow(clebsch_gordan(J1, 1.0, J2, 0.0, 0.0, 0.0), 2.0);    

  return m4;
}


double compute_matrix_element_sigma_0_finite_q(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, double q, int J, int it12) {

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
      m1 *= pow(-1.0, 1.0 + t12)*sqrt(2*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*6.0;
      m1 *= fact;
      m1 *= compute_radial_matrix_element_J_dot_J(J, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12, q);
      m4 += m1;
    }
  }
  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= 1/sqrt(2);}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= 1/sqrt(2);}
                        
  return m4;
}

double compute_total_matrix_element_sigma_0_finite_q_op5(char* density_file, double q, int J1, int J2) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;
 
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {
    double m4 = compute_matrix_element_sigma_0_finite_q_op5(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, q, J1, J2, it12); 
    mat += m4*density;
  }
                        
  return mat;
}


double compute_total_matrix_element_sigma_0_finite_q(char* density_file, double q, int J) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;
 
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {
    double m4 = compute_matrix_element_sigma_0_finite_q(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, q, J, it12); 
    mat += m4*density;
  }
                        
  return mat;
}

double compute_total_matrix_element_y2_finite_q_alt(char* density_file, double q, int J1, int J2) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;
 
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    double m4 = compute_matrix_element_y2_finite_q_alt(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, q, J1, J2); 
    mat += m4*density;
  }
                        
  return mat;
}


double compute_total_matrix_element_y2_finite_q(char* density_file, double q, int J1, int J2) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;
 
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    double m4 = compute_matrix_element_y2_finite_q(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, q, J1, J2); 
    mat += m4*density;
  }
                        
  return mat;
}


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

double compute_matrix_element_y2_finite_q_alt(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J1, int J2) {
  if ((J1 % 2) || (J2 % 2)) {return 0.0;}

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
      double fact = sqrt((2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, 1, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, 1, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      fact *= sqrt(2*lambda + 1)*sqrt(2*lambdap + 1);
      fact *= pow(-1.0, lambda + j12 + 1)*3.0;
      fact *= pow(-1.0, 1.0 + t12)*sqrt(2.0*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*6.0;
      fact *= six_j(j12, 1, lambdap, 2, lambda, 1);
      fact *= 2.0*sqrt(5)*sqrt(2*j12 + 1);
      fact *= compute_radial_matrix_element_Y2_q_alt(J1, J2, n1p, l1p, n2p, l2p, lambdap, n1, l1, n2, l2, lambda, 1, t12, q);
      m4 += fact;
    }
  }
  if ((in1 == in2) && (j1 == j2)) {m4 *= 1/sqrt(2);}
  if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1/sqrt(2);}  
  m4 *= sqrt(2.0*J2 + 1.0)*clebsch_gordan(J2, 2, J1, 0.0, 0.0, 0.0)*sqrt(8.0*M_PI/15.0)/(8.0*M_PI);

  return m4;
}


double compute_matrix_element_y2_finite_q(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J1, int J2) {
  if ((J1 % 2) || (J2 % 2)) {return 0.0;}
  int phase = pow(-1, (J1 - J2)/2.0); 

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
      double fact = sqrt((2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, 1, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, 1, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      fact *= sqrt(2*lambda + 1)*sqrt(2*lambdap + 1);
      fact *= pow(-1.0, j12 + 1)*3.0;
      fact *= six_j(j12, 1, lambdap, 2, lambda, 1);
      fact *= pow(-1.0, 1.0 + t12)*sqrt(2.0*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*6.0;
      fact *= 2.0*sqrt(5)*sqrt(2*j12 + 1);
      fact *= compute_radial_matrix_element_Y2_q(J1, J2, n1p, l1p, n2p, l2p, lambdap, n1, l1, n2, l2, lambda, 1, t12, q);
      m4 += fact;
    }
  }
  if ((in1 == in2) && (j1 == j2)) {m4 *= 1/sqrt(2);}
  if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1/sqrt(2);}  
  m4 *= phase*sqrt((2.0*J1 + 1.0)*(2.0*J2 + 1.0)/5.0)*clebsch_gordan(J1, J2, 2.0, 0.0, 0.0, 0.0)*sqrt(8.0*M_PI/15.0)/(8.0*M_PI);

  return m4;
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


double compute_matrix_element_sigma_tau_plus(char* density_file, int iv) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
  double mat = 0.0;
 
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;

  // Each line of the file corresponds to a nuclear shell
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j2 = ij2/2.0;
    double j12 = ij12/2.0;
    double t12 = it12/2.0;
    double j1p = ij1p/2.0;
    double j2p = ij2p/2.0;
    double j12p = ij12p/2.0;
    double t12p = it12p/2.0;
    if (j12 != j12p) {continue;}

    int l1, l2, l1p, l2p;
     
    if (j1 == 4.5) {l1 = 4;}
    else if (j1 == 2.5) {l1 = 3;}
    else {l1 = 1;}
    if (j2 == 4.5) {l2 = 4;}
    else if (j2 == 2.5) {l2 = 3;}
    else {l2 = 1;}
    if (j1p == 4.5) {l1p = 4;}
    else if (j1p == 2.5) {l1p = 3;}
    else {l1p = 1;}
    if (j2p == 4.5) {l2p = 4;}
    else if (j2p == 2.5) {l2p = 3;}
    else {l2p = 1;}
    
    // The N's listed in the input file are energy quanta, we want radial quantum numbers
    double n1 = (in1)/2.0;
    double n2 = (in2)/2.0;
    double n1p = (in1p)/2.0;
    double n2p = (in2p)/2.0;
 
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
        m1 *= sqrt(5);
        m1 *= fact;
        m1 *= compute_radial_matrix_element_scalar(iv, n1p, l1p, n2p, l2p, n1, l1, n2, l2, lambda, s, t12);
        m4 += m1;
      }
    }
    if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= 1/sqrt(2);}
    if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= 1/sqrt(2);}

    mat += m4*density;
  }
                        
  return mat;
}


double compute_matrix_element_M_GT(char* density_file) {
  // no radial part
  // Computes the total nuclear matrix element for double Gamow-Teller operator sigma(1) dot sigma(2) [tau+(1) tau+(2)]_2
  // Uses the density matrix method to decompose into two-body matrix elements
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  
  // Open the file containing density matrix coefficients
  FILE *in_file;
  in_file = fopen(density_file, "r");
  
  double mat = 0.0;
  int i;
  // Each line of the file corresponds to a nuclear shell
  float density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13){

    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j2 = ij2/2.0;
    double j12 = ij12/2.0;
    double t12 = it12/2.0;
    double j1p = ij1p/2.0;
    double j2p = ij2p/2.0;
    double j12p = ij12p/2.0;
    double t12p = it12p/2.0;
    if (j12 != j12p) {continue;} 
    int lambda, s;
    int l1, l2, l1p, l2p;
  
    if (j1 == 4.5) {l1 = 4;}
    else if (j1 == 2.5) {l1 = 3;}
    else {l1 = 1;}
    if (j2 == 4.5) {l2 = 4;}
    else if (j2 == 2.5) {l2 = 3;}
    else {l2 = 1;}
    if (j1p == 4.5) {l1p = 4;}
    else if (j1p == 2.5) {l1p = 3;}
    else {l1p = 1;}
    if (j2p == 4.5) {l2p = 4;}
    else if (j2p == 2.5) {l2p = 3;}
    else {l2p = 1;}
//    if ((j2 > j1) || (j2p > j1p)) {continue;}
    // The N's listed in the input file are energy quanta, we want radial quantum numbers
    double n1 = (in1 - l1)/2.0;
    double n2 = (in2 - l2)/2.0;
    double n1p = (in1p - l1p)/2.0;
    double n2p = (in2p - l2p)/2.0;
    double m4 = 0.0;
    // Convert from JJ to LS coupling (L is lambda)
    for (lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
      int s_max = MIN(lambda + j12, 1);
      int s_min = abs(lambda - j12);
      for (s = s_min; s <= s_max; s++) {
        double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
        fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
        fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
        fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
        if (fact == 0.0) {continue;}
        fact *= sqrt(2*j12 + 1.0);
        double m1 = pow(-1.0, 1.0 + s)*six_j(s,0.5,0.5,1.0,0.5,0.5)*6.0;
        m1 *= fact;
        m1 *= sqrt(5);
        double anti_symm = 0.0;
        if ((n1 == n1p) && (l1 == l1p) && (n2 == n2p) && (l2 == l2p)) {anti_symm = 1.0;}
        if ((n1 == n2p) && (l1 == l2p) && (n2 == n1p) && (l2 == l1p)) {anti_symm += pow(-1.0, t12 + l1 + l2 + lambda + s + 1);}
        if (anti_symm == 0) {continue;}
        m1 *= anti_symm;
        m4 += m1;
      }
    }
    if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= 1/sqrt(2.0);}
    if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= 1/sqrt(2.0);}

    mat += m4*density;
  }
  fclose(in_file);     
  
  return mat;
}

double compute_matrix_element_M_F(char* density_file) {
  // No radial part
  // Computes the total nuclear matrix element for double Fermi operator [tau+ tau+]_2
  // Uses the density matrix method to decompose into two-body matrix elements
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;

  // Open the file containing density matrix coefficients
  FILE *in_file;
  in_file = fopen(density_file, "r");
  
  double mat = 0.0;
    // Each line of the file corresponds to a nuclear shell
  float density;
  int n_spec;
  while (fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13){
    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j2 = ij2/2.0;
    double j12 = ij12/2.0;
    double t12 = it12/2.0;
    double j1p = ij1p/2.0;
    double j2p = ij2p/2.0;
    double j12p = ij12p/2.0;
    double t12p = it12p/2.0;
    double m4 = 0.0;
    if (t12 != 1 || t12p != 1) {continue;}
    if (j12 != j12p) {continue;}
    if ((in1 == in1p) && (j1 == j1p) && (in2 == in2p) && (j2 == j2p)) {
      m4 = 1.0;
    }
   
    if ((in1 == in2p) && (j1 == j2p) && (in2 == in1p) && (j2 == j1p)) {
      m4 += pow(-1.0, j1 + j2 + j12 + t12);
    }
    
    if (m4 == 0) {continue;}
    if ((in1 == in2) && (j1 == j2)) {m4 *= 1.0/sqrt(2.0);}
    if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1.0/sqrt(2.0);}
    m4 *= sqrt(2.0*j12p + 1.0);
    mat += m4*density*sqrt(5.0);
  }
  fclose(in_file);
                        
  return mat;
}

double compute_matrix_element_M_I1(char *density_file) {
  // No radial part
  // Computes the total nuclear matrix element for one-body identity operator
  int in1, ij1; 
  int in1p, ij1p;

  // Open the file containing density matrix coefficients
  FILE *in_file;
  in_file = fopen(density_file, "r");
  
  double mat = 0.0;
  float density;
  while(fscanf(in_file, "%d %d %d %d %f\n", &in1p, &ij1p, &in1, &ij1, &density) == 5) {
    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j1p = ij1p/2.0;
    double m4 = 0.0;
    if ((in1 == in1p) && (j1 == j1p)) {
      m4 = 1.0;
    } else {
      continue;
    }
    m4 *= sqrt(2.0*j1 + 1.0)*sqrt(2);
    mat += m4*density;
  }
  fclose(in_file);
                        
  return mat;
}

double compute_matrix_element_M_F1(char *density_file) {
  // No radial part
  // Computes the total nuclear matrix element for one-body Fermi operator
  int in1, ij1; 
  int in1p, ij1p;

  // Open the file containing density matrix coefficients
  FILE *in_file;
  in_file = fopen(density_file, "r");
  
  double mat = 0.0;
  float density;
  while(fscanf(in_file, "%d %d %d %d %f\n", &in1p, &ij1p, &in1, &ij1, &density) == 5) {
    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j1p = ij1p/2.0;
    double m4 = 0.0;
    if ((in1 == in1p) && (j1 == j1p)) {
      m4 = 1.0;
    } else {
      continue;
    }
    m4 *= sqrt(2.0*j1 + 1.0)*sqrt(6)/2;
    mat += m4*density;
  }
  fclose(in_file);
                        
  return mat;
}

double compute_matrix_element_M_J1(char *density_file) {
  // No radial part
  // Computes the total nuclear matrix element for one-body Fermi operator
  int in1, ij1; 
  int in1p, ij1p;

  // Open the file containing density matrix coefficients
  FILE *in_file;
  in_file = fopen(density_file, "r");
  
  double mat = 0.0;
  float density;
  while(fscanf(in_file, "%d %d %d %d %f\n", &in1p, &ij1p, &in1, &ij1, &density) == 5) {
    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j1p = ij1p/2.0;
    double m4 = 0.0;
    if ((in1 == in1p) && (j1 == j1p)) {
      m4 = 1.0;
    } else {
      continue;
    }
    m4 *= sqrt(2.0)*sqrt(j1*(j1 + 1)*(2*j1 + 1))/2;
    mat += m4*density;
  }
  fclose(in_file);
                        
  return mat;
}


double compute_matrix_element_M_I2(char *density_file) {
  // No radial part
  // Computes the total nuclear matrix element for two-body identity operator 1x1
  // Uses the density matrix method to decompose into two-body matrix elements
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;

  // Open the file containing density matrix coefficients
  FILE *in_file;
  in_file = fopen(density_file, "r");
  
  double mat = 0.0;
  int i;
  float density;
  int n_spec;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {
    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j2 = ij2/2.0;
    double j12 = ij12/2.0;
    double t12 = it12/2.0;
    double j1p = ij1p/2.0;
    double j2p = ij2p/2.0;
    double j12p = ij12p/2.0;
    double t12p = it12p/2.0;
    double m4 = 0.0;
    if (t12 != t12p) {continue;}
    if (j12 != j12p) {continue;}
    if ((in1 == in1p) && (j1 == j1p) && (in2 == in2p) && (j2 == j2p)) {
      m4 = 1.0;
    }
   
    if ((in1 == in2p) && (j1 == j2p) && (in2 == in1p) && (j2 == j1p)) {
      m4 += pow(-1.0, j1 + j2 + j12 + t12);
    }
    if (m4 == 0) {continue;}
    if ((in1 == in2) && (j1 == j2)) {m4 *= 1.0/sqrt(2.0);}
    if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1.0/sqrt(2.0);}
    m4 *= sqrt(2.0*j12p + 1.0)*sqrt(2*t12 + 1);
    mat += m4*density;
  }
  fclose(in_file);
                        
  return mat;
}


double compute_matrix_element_M_J2(char* density_file) {
  // No radial part
  // Computes the total nuclear matrix element for the operator [J x J]_2
  // Uses the density matrix method to decompose into two-body matrix elements
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;

  // Open the file containing density matrix coefficients
  FILE *in_file;
  in_file = fopen(density_file, "r");
  
  double mat = 0.0;
  int i;
  // Each line of the file corresponds to a nuclear shell
  double density;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j2 = ij2/2.0;
    double j12 = ij12/2.0;
    double t12 = it12/2.0;
    double j1p = ij1p/2.0;
    double j2p = ij2p/2.0;
    double j12p = ij12p/2.0;
    double t12p = it12p/2.0;
    double m4 = 0.0;
    if (t12 != t12p) {continue;}
    if ((in1 == in1p) && (j1 == j1p) && (in2 == in2p) && (j2 == j2p)) {
      m4 = nine_j(j1p, j1, 1, j2p, j2, 1, j12p, j12, 2);
    }
   
    if ((in1 == in2p) && (j1 == j2p) && (in2 == in1p) && (j2 == j1p)) {
      m4 += pow(-1.0,  j1 + j2 + j12 + t12)*nine_j(j1p, j2, 1, j2p, j1, 1, j12p, j12, 2);
    }
    if (m4 == 0) {continue;}
    m4 *= sqrt(j1*(2*j1 + 1)*(j1 + 1))*sqrt(j2*(2*j2 + 1)*(j2 + 1))*sqrt(5);
    m4 *= sqrt(2.0*j12 + 1.0)*sqrt(2*j12p + 1.0);
    
    if ((in1 == in2) && (j1 == j2)) {m4 *= 1.0/sqrt(2.0);}
    if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1.0/sqrt(2.0);}
    mat += m4*density*sqrt(2*t12 + 1);
  }
  fclose(in_file);
                        
  return mat;
}

