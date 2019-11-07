#include "wfn.h"

void mu2e() {

  FILE* in_file;
  in_file = fopen("dirac_bound_mu.dat", "r");
  int np_mu = 359;
  int np_e = 2000;
  double a0mu = 255.9; // [fm]
  double a0e = 52917.7; // [fm]
  double *r_arr_mu = (double*) malloc(sizeof(double)*np_mu);
  double *g_arr_mu = (double*) malloc(sizeof(double)*np_mu);
  double *f_arr_mu = (double*) malloc(sizeof(double)*np_mu);
  
  double *r_arr_e = (double*) malloc(sizeof(double)*np_e);
  double *g_arr_e = (double*) malloc(sizeof(double)*np_e);
  double *f_arr_e = (double*) malloc(sizeof(double)*np_e);

  for (int i = 0; i < np_mu; i++) {
    double r, f, g;
    fscanf(in_file, "%lf %lf %lf\n", &r, &g, &f);
    r_arr_mu[i] = a0mu*r;
    g_arr_mu[i] = g/sqrt(a0mu);
    f_arr_mu[i] = f/sqrt(a0mu);
  }
  fclose(in_file);
 
  in_file = fopen("dirac_free_e.dat", "r");
  
  for (int i = 0; i < np_e; i++) {
    double r, f, g;
    fscanf(in_file, "%lf %lf %lf\n", &r, &g, &f);
    r_arr_e[i] = a0e*r;
    g_arr_e[i] = g/sqrt(a0e)*sqrt(2.0*ALPHA_FS)/sqrt(27.2114*pow(10, -6));
    f_arr_e[i] = f/sqrt(a0e)*sqrt(2.0*ALPHA_FS)/sqrt(27.2114*pow(10, -6));
  }
  fclose(in_file);


  gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
  gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
  gsl_interp_accel *acc3 = gsl_interp_accel_alloc();
  gsl_interp_accel *acc4 = gsl_interp_accel_alloc();

  gsl_spline *g_mu = gsl_spline_alloc(gsl_interp_cspline, np_mu);
  gsl_spline *f_mu = gsl_spline_alloc(gsl_interp_cspline, np_mu);
  gsl_spline *g_e = gsl_spline_alloc(gsl_interp_cspline, np_e);
  gsl_spline *f_e = gsl_spline_alloc(gsl_interp_cspline, np_e);

  gsl_spline_init(g_mu, r_arr_mu, g_arr_mu, np_mu);
  gsl_spline_init(f_mu, r_arr_mu, f_arr_mu, np_mu);
  gsl_spline_init(g_e, r_arr_e, g_arr_e, np_e);
  gsl_spline_init(f_e, r_arr_e, f_arr_e, np_e);

  double mat = sqrt(2.0)*(1.799*RombergSpline41Integrator(&newlep1, g_mu, f_mu, g_e, f_e, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001) + 1.21246*RombergSpline41Integrator(&newlep2, g_mu, f_mu, g_e, f_e, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001) + 8.97138*RombergSpline41Integrator(&newlep3, g_mu, f_mu, g_e, f_e, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001));
  mat *= -32.0*M_PI*sqrt(105.657)/(2.0*sqrt(M_PI))*clebsch_gordan(0.5, 0.5, 0, -0.5, 0.5, 0)/sqrt(6.0);
  double core_mat = 2.0*sqrt(2.0)*(sqrt(2.0)*RombergSpline41Integrator(&core1, g_mu, f_mu, g_e, f_e, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001) + sqrt(2.0)*RombergSpline41Integrator(&core2, g_mu, f_mu, g_e, f_e, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001) + 2*RombergSpline41Integrator(&core3, g_mu, f_mu, g_e, f_e, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001));

  core_mat *= 32.0*M_PI/(2.0*sqrt(M_PI))*sqrt(105.657);

  printf("T = 0: %g\n", mat + core_mat);
  mat = sqrt(6.0)*(0.33169*RombergSpline41Integrator(&newlep1, g_mu, f_mu, g_e, f_e, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001) + 0.17295*RombergSpline41Integrator(&newlep2, g_mu, f_mu, g_e, f_e, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001) + 0.66729*RombergSpline41Integrator(&newlep3, g_mu, f_mu, g_e, f_e, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001));
  mat *= 32.0*M_PI*sqrt(105.657)/(2.0*sqrt(M_PI))*clebsch_gordan(0.5, 0.5, 1, -0.5, 0.5, 0)/sqrt(3.0*(2*5/2+1));
  printf("T = 1: %g\n", mat);

  double *sigma_mat = (double*) calloc(36, sizeof(double));

  sigma_mat[7] = RombergSplineFunIntegrator(&lep_int_w1, &SigmaJ1_2s1_2s1, g_mu, f_mu, g_e, f_e, acc1, acc2, acc3, acc4, 0.0001, 30.0, 0.00001);

   sigma_mat[9] = RombergSplineFunIntegrator(&lep_int_w1, &SigmaJ1_1d3_2s1, g_mu, f_mu, g_e, f_e, acc1, acc2, acc3, acc4, 0.0001, 30.0, 0.00001);
 
   sigma_mat[19] = -RombergSplineFunIntegrator(&lep_int_w1, &SigmaJ1_1d3_2s1, g_mu, f_mu, g_e, f_e, acc1, acc2, acc3, acc4, 0.0001, 30.0, 0.00001);
 
   sigma_mat[21] = RombergSplineFunIntegrator(&lep_int_w1, &SigmaJ1_1d3_1d3, g_mu, f_mu, g_e, f_e, acc1, acc2, acc3, acc4, 0.0001, 30.0, 0.00001);
 
   sigma_mat[23] = RombergSplineFunIntegrator(&lep_int_w1, &SigmaJ1_1d5_1d3, g_mu, f_mu, g_e, f_e, acc1, acc2, acc3, acc4, 0.0001, 30.0, 0.00001);
 
   sigma_mat[33] = -RombergSplineFunIntegrator(&lep_int_w1, &SigmaJ1_1d5_1d3, g_mu, f_mu, g_e, f_e, acc1, acc2, acc3, acc4, 0.0001, 30.0, 0.00001);
 
   sigma_mat[35] = RombergSplineFunIntegrator(&lep_int_w1, &SigmaJ1_1d5_1d5, g_mu, f_mu, g_e, f_e, acc1, acc2, acc3, acc4, 0.0001, 30.0, 0.00001);


 printf("%g %g %g %g %g %g %g\n", sigma_mat[7], sigma_mat[9], sigma_mat[19], sigma_mat[21], sigma_mat[23], sigma_mat[33], sigma_mat[35]);
 
 // for (int i = 0; i < 300; i++) {
 //   printf("%g, %g, %g, %g\n", i/20.0, SigmaJ0_1p1_1s1(i/20.0 + 0.00001), SigmaJ0_2s1_1p1(i/20.0 + 0.00001), SigmaJ0_1d3_1p3(i/20.0 + 0.00001));
 // }
  printf("Test: %g\n", SigmaJ5_1d5_1d5(2.0));
  in_file = fopen("al27-al27_core_J0_T0_0_0.dens", "r");


  float density;
  mat = 0.0;
  int in1p, ij1p, in2p, ij2p, ij12p, it12p, in1, ij1, in2, ij2, ij12, it12;
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
    if (t12 != t12p) {continue;}
    // Compute J = 1 operators
    mat += density*nine_j(j1p, j1, 1, j2p, j2, 1, j12p, j12, 0)*pow(-1.0, 1 + t12)*sqrt(2.0*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*(2*j12 + 1)*6;
      //*sigma_mat[ij1p + 6*ij1]*sigma_mat[ij2p + 6*ij2]
  //  if (mat != 0.0) {printf("n1p: %d j1p: %g n2p: %d j2p: %g n1: %d j1: %g n2: %d j2: %g\n", in1p, j1p, in2p, j2p, in1, j1, in2, j2);}
  }
  mat *= 64.0*113.06/sqrt(4.0*M_PI)*(-1)*clebsch_gordan(0.5, 0.5, 0.0, -0.5, 0.5, 0.0)/sqrt(6.0);
  mat *= 3.0*clebsch_gordan(1, 1, 0, 0, 0, 0);
  printf("Final: %g\n", mat);
 

 
  return;
}

double wfn_sq(gsl_spline* g_spline, gsl_spline* f_spline, gsl_interp_accel *acc, double r) {
  double psi2 = pow(gsl_spline_eval(g_spline, r, acc), 2.0) + pow(gsl_spline_eval(f_spline, r, acc), 2.0);

  return psi2;
}

double newlep1(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

double glep = 1.0/(3.0*sqrt(2.0))/pow(105.658*b_osc(27)/197.3, 3.0)*exp(-pow(r/b_osc(27), 2.0))*pow(3.0 - 2*pow(r/b_osc(27) , 2.0), 2.0)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}

double core1(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

double glep = sqrt(2.0)/pow(105.658*b_osc(27)/197.3, 3.0)*exp(-pow(r/b_osc(27), 2.0))*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}

double core2(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

double glep = 2.0*sqrt(2.0)/3.0/pow(105.658*b_osc(27)/197.3, 3.0)*exp(-pow(r/b_osc(27), 2.0))*pow(r/b_osc(27), 2.0)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}

double core3(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

double glep = 4.0/3.0/pow(105.658*b_osc(27)/197.3, 3.0)*exp(-pow(r/b_osc(27), 2.0))*pow(r/b_osc(27), 2.0)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}


double newlep2(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

  double glep = 8.0/(15.0)/pow(105.658*b_osc(27)/197.3, 3.0)*exp(-pow(r/b_osc(27), 2.0))*pow(r/b_osc(27) , 4.0)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}

double newlep3(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

  double glep = 4.0/(5.0)*sqrt(2.0/3.0)/pow(105.658*b_osc(27)/197.3, 3.0)*exp(-pow(r/b_osc(27), 2.0))*pow(r/b_osc(27) , 4.0)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}



double g1lep(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r, double q) {
  double glep = 1.0/(2.0*sqrt(M_PI))*sqrt(105.658)*gsl_sf_bessel_j0(q*r/197.3)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));
  return glep;
}

double g1lep2(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r, double q1, double q2) {
  double glep = 1.0/(2.0*sqrt(M_PI))*sqrt(105.658)*gsl_sf_bessel_j0(q1*r/197.3)*gsl_sf_bessel_j0(q2*r/197.3)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));
  return glep;
}

double lep_int_w1(double (*f) (double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

  double glep = 1.0/(2.0*sqrt(M_PI))*sqrt(105.658)*f(r)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}

double M_op_sd(int ijp, int ij, int j_op, double q, double b) {
  double y = pow(b*q/(2.0*197.3), 2.0);
  printf("b: %g y: %g\n", b, y);
  double mat = pow(4.0*M_PI, -0.5)*pow(y, (j_op - 2.0)/2.0)*exp(-y);

  if (j_op == 0) {
    if ((ij == 1) && (ijp == 1)) {
      mat *= sqrt(2.0)*(y - 4.0/3.0*y*y + 2.0/3.0*y*y*y);
    } else if ((ij == 3) && (ijp == 3)) {
      mat *= 2.0*(y - 4.0/3.0*y*y + 4.0/15.0*y*y*y);
    } else if ((ij == 5) && (ijp == 5)) {
      mat *= sqrt(6.0)*(y - 4.0/3.0*y*y + 4.0/15.0*y*y*y);
    } else {
      mat = 0.0;
    }
  } else if (j_op == 2) {
    if ((ij == 1) && (ijp == 3)) {
      mat *= sqrt(10.0)*2.0/15.0*(y - 0.5*y*y);
    } else if ((ij == 3) && (ijp == 1)) {
      mat *= -sqrt(10.0)*2.0/15.0*(y - 0.5*y*y);
    } else if ((ij == 3) && (ijp == 3)) {
      mat *= 28.0/15.0*(-y + 2.0/7.0*y*y);
    } else if ((ijp == 5) && (ij == 3)) {
      mat *= sqrt(21.0)*4.0/15.0*(-y + 2.0/7.0*y*y);
    } else if ((ijp == 3) && (ij == 5)) {
      mat *= -sqrt(21.0)*4.0/15.0*(-y + 2.0/7.0*y*y);
    } else if ((ijp == 5) && (ij == 5)) {
      mat *= sqrt(21.0)*8.0/15.0*(-y + 2.0/7.0*y*y);
    } else if ((ijp = 5) && (ij == 1)) {
      mat *= sqrt(15.0)*8.0/15.0*(-y + 0.5*y*y);
    } else {
      mat = 0.0;
    }
  } else if (j_op == 4) {
    if ((ijp == 5) && (ij == 3)) {
      mat *= sqrt(14.0)*8.0/35.0*y;
    } else if ((ijp == 3) && (ij == 5)) {
      mat *= -sqrt(14.0)*8.0/35.0*y;
    } else if ((ijp == 5) && (ij == 5)) {
      mat *= sqrt(7.0)*8.0/35.0*y;
    } else {
      mat *= 0.0;
    }
  } else {
    mat *= 0.0;
  }
  return mat;
}

double MJ0_s1_s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 2.0/pow(w, 3.0)*exp(-pow(z/w, 2.0))*(w*w + 2.0*pow(w, 4) + 2.0*z*z);
  mat += sqrt(M_PI)/z*exp(w*w)*(3.0 + 4.0*w*w + 2.0*pow(w, 4))*(exp(2.0*z)*gsl_sf_erfc(w + z/w) - exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *= -1.0/(24.0*sqrt(2.0));

  return mat;
}

double MJ0_d3_d3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/pow(w, 3.0)*exp(-pow(z/w, 2.0))*(7.0*w*w + 2.0*pow(w, 4) + 2.0*z*z);
  mat += sqrt(M_PI)/z*exp(w*w)*(15.0 + 20.0*w*w + 4.0*pow(w, 4))*(exp(2.0*z)*gsl_sf_erfc(w + z/w) - exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *= -1.0/(120.0);

  return mat;
}

double MJ0_d5_d5(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/pow(w, 3.0)*exp(-pow(z/w, 2.0))*(7.0*w*w + 2.0*pow(w, 4) + 2.0*z*z);
  mat += sqrt(M_PI)/z*exp(w*w)*(15.0 + 20.0*w*w + 4.0*pow(w, 4))*(exp(2.0*z)*gsl_sf_erfc(w + z/w) - exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *= -sqrt(6.0)/(240.0);

  return mat;
}

double SigmaJ1_2s1_2s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*pow(w, 5))*exp(-pow(z/w, 2.0))*(4.0*pow(w, 6) + 2.0*pow(w, 8) - w*w*z*z + 2.0*pow(z, 4) + 3.0*pow(w, 4) + 2.0*pow(w, 4)*z*z);
  mat += sqrt(M_PI)/(z*z)*exp(w*w)*(3.0 + 4.0*w*w + 2.0*pow(w, 4))*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(48.0*sqrt(2.0));

  return mat;
}

double SigmaJ1_1d3_2s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -2.0/(z*pow(w, 5))*exp(-pow(z/w, 2.0))*(4.0*pow(w, 6) + 2.0*pow(w, 8) - w*w*z*z + 2.0*pow(z, 4) + 2.0*pow(w, 4)*z*z);
  mat += sqrt(M_PI)/(z*z)*exp(w*w)*(2.0*w*w + pow(w, 4))*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(12.0*sqrt(10.0));

  return mat;
}

double SigmaJ1_1d3_1d3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*pow(w, 5))*exp(-pow(z/w, 2.0))*(8.0*pow(w, 6) + 4.0*pow(w, 8) - 2.0*w*w*z*z + 4.0*pow(z, 4) - 15.0*pow(w, 4) + 4.0*pow(w, 4)*z*z);
  mat += sqrt(M_PI)/(z*z)*exp(w*w)*(-15.0 + 8.0*w*w + 4.0*pow(w, 4))*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(240.0*sqrt(5.0));

  return mat;
}

double SigmaJ1_1d5_1d3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/(z*pow(w, 5))*exp(-pow(z/w, 2.0))*(9.0*pow(w, 6) + 2.0*pow(w, 8) + 4.0*w*w*z*z + 2.0*pow(z, 4) + 5.0*pow(w, 4) + 2.0*pow(w, 4)*z*z);
  mat -= sqrt(M_PI)/(z*z)*exp(w*w)*(5.0 + 9.0*w*w + 2.0*pow(w, 4))*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(40.0*sqrt(5.0));

  return mat;
}

double SigmaJ1_1d5_1d5(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*pow(w, 5))*exp(-pow(z/w, 2.0))*(28.0*pow(w, 6) + 4.0*pow(w, 8) + 18.0*w*w*z*z + 4.0*pow(z, 4) + 35.0*pow(w, 4) + 4.0*pow(w, 4)*z*z);
  mat += sqrt(M_PI)/(z*z)*exp(w*w)*(35.0 + 28.0*w*w + 4.0*pow(w, 4))*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(80.0*sqrt(70.0));

  return mat;
}

double SigmaJ1_1s1_1s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*w)*exp(-pow(z/w, 2.0));
  mat += sqrt(M_PI)/(z*z)*exp(w*w)*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(16.0*sqrt(2.0));

  return mat;
}

double SigmaJ1_1p1_1p1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*pow(w, 3))*exp(-pow(z/w, 2.0))*(2.0*pow(w, 4) + 2.0*z*z - w*w);
  mat += sqrt(M_PI)/(z*z)*exp(w*w)*(2.0*w*w - 1.0)*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(48.0*sqrt(2.0));

  return mat;
}

double SigmaJ1_1p3_1p1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/(z*pow(w, 3))*exp(-pow(z/w, 2.0))*(pow(w, 4) + z*z + w*w);
  mat -= sqrt(M_PI)/(z*z)*exp(w*w)*(w*w + 1.0)*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(24.0);

  return mat;
}

double SigmaJ1_1p3_1p3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*pow(w, 3))*exp(-pow(z/w, 2.0))*(2.0*pow(w, 4) + 2.0*z*z + 5.0*w*w);
  mat += sqrt(M_PI)/(z*z)*exp(w*w)*(2.0*w*w + 5.0)*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(48.0*sqrt(5.0));

  return mat;
}

double SigmaJ0_1p1_1s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 2.0/(pow(w, 2))*exp(-pow(z/w, 2.0));
  mat += sqrt(M_PI)/(z)*exp(w*w)*w*(exp(2*z)*gsl_sf_erfc(w + z/w) - exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(8.0*sqrt(3.0));

  return mat;
}

double SigmaJ2_1p3_1s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(pow(z*w, 2))*exp(-pow(z/w, 2.0))*(2.0*z*z + 3.0*w*w);
  mat += sqrt(M_PI)/(z*z*z)*exp(w*w)*w*(-exp(2*z)*(3.0 - 6.0*z + 4.0*z*z)*gsl_sf_erfc(w + z/w) + (3.0 + 6.0*z + 4.0*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(16.0*sqrt(6.0));

  return mat;
}

double SigmaJ3_1d3_1d3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/(w*w*pow(z*w, 3))*exp(-pow(z/w, 2.0))*(15.0*pow(w, 6)*(2*w*w - 1.0) + 2*pow(w, 4)*(-5 + 8*w*w + 4*pow(w, 4))*z*z + 4*w*w*(2*w*w - 1.0)*pow(z, 4) + 8*pow(z, 6));
  mat -= sqrt(M_PI)/(z*z*z*z)*exp(w*w)*w*w*(2*w*w - 1.0)*(exp(2*z)*(-15.0 + 30*z - 24*z*z + 8*z*z*z)*gsl_sf_erfc(w + z/w) + (15 + 30*z + 24*z*z + 8*z*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(160.0*sqrt(5.0));

  return mat;
}

double SigmaJ3_1d5_2s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/(w*w*pow(z*w, 3))*exp(-pow(z/w, 2.0))*(15.0*pow(w, 6)*(w*w + 2.0) + 2*pow(w, 4)*(2 + w*w)*(5 + 2*w*w)*z*z + 4*w*w*(2 + w*w)*pow(z, 4) + 4*pow(z, 6));
  mat -= sqrt(M_PI)/(z*z*z*z)*exp(w*w)*w*w*(w*w + 2.0)*(exp(2*z)*(-15.0 + 30*z - 24*z*z + 8*z*z*z)*gsl_sf_erfc(w + z/w) + (15 + 30*z + 24*z*z + 8*z*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(32.0*sqrt(15.0));

  return mat;
}

double SigmaJ3_1d5_1d3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(w*w*pow(z*w, 3))*exp(-pow(z/w, 2.0))*(15.0*pow(w, 6)*(w*w + 2.0) + 2*pow(w, 4)*(2 + w*w)*(5 + 2*w*w)*z*z + 4*w*w*(2 + w*w)*pow(z, 4) + 4*pow(z, 6));
  mat += sqrt(M_PI)/(z*z*z*z)*exp(w*w)*w*w*(w*w + 2.0)*(exp(2*z)*(-15.0 + 30*z - 24*z*z + 8*z*z*z)*gsl_sf_erfc(w + z/w) + (15 + 30*z + 24*z*z + 8*z*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(40.0*sqrt(30.0));

  return mat;
}

double SigmaJ3_1d5_1d5(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/(w*w*pow(z*w, 3))*exp(-pow(z/w, 2.0))*(15.0*pow(w, 6)*(2*w*w + 9.0) + 2*pow(w, 4)*(45 + 28*w*w + 4*pow(w, 4))*z*z + 4*w*w*(9 + 2*w*w)*pow(z, 4) + 8*pow(z, 6));
  mat -= sqrt(M_PI)/(z*z*z*z)*exp(w*w)*w*w*(2*w*w + 9.0)*(exp(2*z)*(-15.0 + 30*z - 24*z*z + 8*z*z*z)*gsl_sf_erfc(w + z/w) + (15 + 30*z + 24*z*z + 8*z*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(240.0*sqrt(5.0));

  return mat;
}

double SigmaJ5_1d5_1d5(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(pow(z*w, 5))*exp(-pow(z/w, 2.0))*(945*pow(w, 8) + 210*pow(w, 6)*(3 + 2*w*w)*z*z + 4*pow(w, 4)*(63 + 28*w*w + 4*pow(w, 4))*pow(z, 4) + 8*w*w*(9 + 2*w*w)*pow(z, 6) + 16*pow(z, 8));
  mat += sqrt(M_PI)/(pow(z, 6))*exp(w*w)*pow(w, 4)*(exp(2*z)*(-945 + 1890*z - 1680*z*z + 840*z*z*z - 240*pow(z, 4) + 32*pow(z, 5))*gsl_sf_erfc(w + z/w) + (945 + 1890*z + 1680*z*z + 840*z*z*z + 240*pow(z, 4) + 32*pow(z, 5))*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(192.0*sqrt(7.0));

  return mat;
}

double SigmaJ0_2s1_1p1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = exp(-pow(z/w, 2.0))/pow(w, 4)*(w*w - 2*pow(w, 4) - 2*z*z);
  mat -= sqrt(M_PI)/z*exp(w*w)*w*(1 + w*w)*(exp(2*z)*gsl_sf_erfc(w + z/w) - exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *= 1.0/(12.0*sqrt(2.0));

  return mat;
}

double SigmaJ0_1d3_1p3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0*exp(-pow(z/w, 2.0))/pow(w, 4)*(w*w + pow(w, 4) +z*z);
  mat += sqrt(M_PI)/z*exp(w*w)*w*(5 + 2*w*w)*(exp(2*z)*gsl_sf_erfc(w + z/w) - exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *= 1.0/(12.0*sqrt(10.0));

  return mat;
}

double SigmaJ1_2s1_1s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0*exp(-pow(z/w, 2.0))/(w*w*w*z)*(pow(w, 4) + z*z);
  mat -= sqrt(M_PI)/(z*z)*exp(w*w)*w*w*(exp(2*z)*(2*z-1)*gsl_sf_erfc(w + z/w) + exp(-2*z)*(2*z + 1)*gsl_sf_erfc(w - z/w));

  mat *= 1.0/(16.0*sqrt(3.0));

  return mat;
}

double SigmaJ1_1d3_1s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0*exp(-pow(z/w, 2.0))/(w*w*w*z)*(pow(w, 4) + z*z);
  mat -= sqrt(M_PI)/(z*z)*exp(w*w)*w*w*(exp(2*z)*(2*z-1)*gsl_sf_erfc(w + z/w) + exp(-2*z)*(2*z + 1)*gsl_sf_erfc(w - z/w));

  mat *= 1.0/(8.0*sqrt(15.0));

  return mat;
}

double SigmaJ3_1d5_1s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(pow(z*w, 3))*exp(-pow(z/w, 2.0))*(10*w*w*z*z + 4*pow(z, 4) + 15*pow(w, 4) + 4*pow(w, 4)*z*z);
  mat += sqrt(M_PI)/(z*z*z*z)*exp(w*w)*w*w*(exp(2*z)*(-15.0 + 30*z - 24*z*z + 8*z*z*z)*gsl_sf_erfc(w + z/w) + (15 + 30*z + 24*z*z + 8*z*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(32.0*sqrt(10.0));

  return mat;
}


double SigmaJ2_2s1_1p3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*z*pow(w, 4))*exp(-pow(z/w, 2.0))*(3*pow(w, 4)*(1 + w*w) + 2*w*w*(1+ w*w)*z*z + 2*pow(z, 4));
  mat += sqrt(M_PI)/(z*z*z)*exp(w*w)*w*(1+ w*w)*(exp(2*z)*(-3.0 + 6*z -4*z*z)*gsl_sf_erfc(w + z/w) + (3.0 + 6*z + 4*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(48.0);

  return mat;
}

double SigmaJ2_1d3_1p1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*z*pow(w, 4))*exp(-pow(z/w, 2.0))*(6.0*pow(w, 6) - 2*w*w*z*z + 4*pow(z, 4) - 3*pow(w, 4) + 4*pow(w, 4)*z*z);
  mat += sqrt(M_PI)/(z*z*z)*exp(w*w)*w*(-1+ 2*w*w)*(exp(2*z)*(-3.0 + 6*z -4*z*z)*gsl_sf_erfc(w + z/w) + (3.0 + 6*z + 4*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(48.0*sqrt(10));

  return mat;
}

double SigmaJ2_1d3_1p3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*z*pow(w, 4))*exp(-pow(z/w, 2.0))*(3*pow(w, 4)*(1 + w*w) + 2*w*w*(1+ w*w)*z*z + 2*pow(z, 4));
  mat += sqrt(M_PI)/(z*z*z)*exp(w*w)*w*(1+ w*w)*(exp(2*z)*(-3.0 + 6*z -4*z*z)*gsl_sf_erfc(w + z/w) + (3.0 + 6*z + 4*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(24.0*sqrt(10));

  return mat;
}

double SigmaJ2_1d5_1p1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/(z*z*pow(w, 4))*exp(-pow(z/w, 2.0))*(3*pow(w, 4)*(2 + w*w) + 2*w*w*(2 + w*w)*z*z + 2*pow(z, 4));
  mat -= sqrt(M_PI)/(z*z*z)*exp(w*w)*w*(2 + w*w)*(exp(2*z)*(-3.0 + 6*z -4*z*z)*gsl_sf_erfc(w + z/w) + (3.0 + 6*z + 4*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(16.0*sqrt(15.0));

  return mat;
}

double SigmaJ2_1d5_1p3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*z*pow(w, 4))*exp(-pow(z/w, 2.0))*(6*pow(w, 6) + 14*w*w*z*z + 4*pow(z, 4) + 21*pow(w, 4) + 4*pow(w, 4)*z*z);
  mat += sqrt(M_PI)/(z*z*z)*exp(w*w)*w*(7 + 2*w*w)*(exp(2*z)*(-3.0 + 6*z -4*z*z)*gsl_sf_erfc(w + z/w) + (3.0 + 6*z + 4*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(16.0*sqrt(210.0));

  return mat;
}

double SigmaJ3_1p3_1p3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/(pow(z*w, 3))*exp(-pow(z/w, 2.0))*(10*w*w*z*z + 4*pow(z, 4) + 15*pow(w, 4) + 4*pow(w, 4)*z*z);
  mat -= sqrt(M_PI)/(z*z*z*z)*exp(w*w)*w*w*(exp(2*z)*(-15.0 + 30*z - 24*z*z + 8*z*z*z)*gsl_sf_erfc(w + z/w) + (15 + 30*z + 24*z*z + 8*z*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(32.0*sqrt(5.0));

  return mat;
}

double SigmaJ4_1d5_1p3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/(pow(z*w, 4))*exp(-pow(z/w, 2.0))*(28*w*w*pow(z, 4) + 8*pow(z, 6) + 5*pow(w, 6)*(21 + 8*z*z) + pow(w, 4)*(70*z*z + 8*pow(z, 4)));
  mat -= sqrt(M_PI)/pow(z, 5)*exp(w*w)*w*w*w*(exp(2*z)*(-105 + 210*z - 180*z*z + 80*z*z*z - 16*pow(z, 4))*gsl_sf_erfc(w + z/w) + (105 + 210*z + 180*z*z + 80*z*z*z + 16*pow(z, 4))*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(32.0*sqrt(35.0));

  return mat;
}


double b_osc(int a_nuc) {
  double b = sqrt(0.9*pow(a_nuc, 1.0/3.0) + 0.7);
  
  return b;
}
