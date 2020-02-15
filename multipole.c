#include "multipole.h"

void compute_nuclear_multipoles() {
  FILE *in_file;
  int j_op = 0;
  int t_op = 0;
  float density;
  int num_shells = 6;
  double *rhoJ0T0 = (double*) calloc(pow(num_shells, 2), sizeof(double));
  double *rhoJ0T1 = (double*) calloc(pow(num_shells, 2), sizeof(double));
  double *rhoJ1T0 = (double*) calloc(pow(num_shells, 2), sizeof(double));
  double *rhoJ1T1 = (double*) calloc(pow(num_shells, 2), sizeof(double));
  double *rhoJ2T0 = (double*) calloc(pow(num_shells, 2), sizeof(double));
  double *rhoJ2T1 = (double*) calloc(pow(num_shells, 2), sizeof(double));
  double *rhoJ3T0 = (double*) calloc(pow(num_shells, 2), sizeof(double));
  double *rhoJ3T1 = (double*) calloc(pow(num_shells, 2), sizeof(double));
  double *rhoJ4T0 = (double*) calloc(pow(num_shells, 2), sizeof(double));
  double *rhoJ4T1 = (double*) calloc(pow(num_shells, 2), sizeof(double));
  double *rhoJ5T0 = (double*) calloc(pow(num_shells, 2), sizeof(double));
  double *rhoJ5T1 = (double*) calloc(pow(num_shells, 2), sizeof(double));

  int in, ij, inp, ijp;

  // Compute Al27 5/2+ (gs) -> 5/2+ (gs) transition

  in_file = fopen("isotope_data/al27/density/al27-al27_core_1bdy_J0_T0_0_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ0T0[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_core_1bdy_J0_T1_0_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ0T1[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J1_T0_0_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ1T0[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J1_T1_0_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ1T1[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J2_T0_0_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ2T0[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J2_T1_0_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ2T1[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);
  
  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J3_T0_0_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ3T0[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J3_T1_0_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ3T1[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J4_T0_0_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ4T0[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J4_T1_0_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ4T1[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J5_T0_0_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ5T0[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J5_T1_0_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ5T1[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  FILE *out_file;
  double y_mu = pow(105.0*b_osc(27)/(HBARC*2.0), 2.0);

  out_file = fopen("MJ_q_al27_0_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, m0_tot(y, rhoJ0T0, num_shells), m1_tot(y, rhoJ1T0, num_shells), m2_tot(y, rhoJ2T0, num_shells), m3_tot(y, rhoJ3T0, num_shells), m4_tot(y, rhoJ4T0, num_shells), 0.0);
  }
  fclose(out_file);

  printf("M0(m_mu): %g M1(m_mu): %g M2(m_mu): %g M3(m_mu): %g M4(m_mu): %g\n", m0_tot(y_mu, rhoJ0T0, num_shells), m1_tot(y_mu, rhoJ1T0, num_shells), m2_tot(y_mu, rhoJ2T0, num_shells), m3_tot(y_mu, rhoJ3T0, num_shells), m4_tot(y_mu, rhoJ4T0, num_shells));


  out_file = fopen("DPJ_q_al27_0_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, 0.0, deltap1_tot(y, rhoJ1T0, num_shells), deltap2_tot(y, rhoJ2T0, num_shells), deltap3_tot(y, rhoJ3T0, num_shells), 0.0, 0.0);
  }
  fclose(out_file);

  printf("DP1(m_mu): %g DP2(m_mu): %g DP3(m_mu): %g\n", deltap1_tot(y_mu, rhoJ1T0, num_shells), deltap2_tot(y_mu, rhoJ2T0, num_shells), deltap3_tot(y_mu, rhoJ3T0, num_shells));


  out_file = fopen("SJ_q_al27_0_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, 0.0, sigma1_tot(y, rhoJ1T0, num_shells), sigma2_tot(y, rhoJ2T0, num_shells), sigma3_tot(y, rhoJ3T0, num_shells), sigma4_tot(y, rhoJ4T0, num_shells), 0.0);
  }
  fclose(out_file);

  printf("S1(m_mu): %g S2(m_mu): %g S3(m_mu): %g S4(m_mu): %g\n", sigma1_tot(y_mu, rhoJ1T0, num_shells), sigma2_tot(y_mu, rhoJ2T0, num_shells), sigma3_tot(y_mu, rhoJ3T0, num_shells), sigma4_tot(y_mu, rhoJ4T0, num_shells));


  out_file = fopen("DJ_q_al27_0_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, 0.0, delta1_tot(y, rhoJ1T0, num_shells), delta2_tot(y, rhoJ2T0, num_shells), delta3_tot(y, rhoJ3T0, num_shells), 0.0, 0.0);
  }
  fclose(out_file);

  printf("D1(m_mu): %g D2(m_mu): %g D3(m_mu): %g\n", delta1_tot(y_mu, rhoJ1T0, num_shells), delta2_tot(y_mu, rhoJ2T0, num_shells), delta3_tot(y_mu, rhoJ3T0, num_shells));


  out_file = fopen("SPJ_q_al27_0_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, 0.0, sigmap1_tot(y, rhoJ1T0, num_shells), sigmap2_tot(y, rhoJ2T0, num_shells), sigmap3_tot(y, rhoJ3T0, num_shells), sigmap4_tot(y, rhoJ4T0, num_shells), sigmap5_tot(y, rhoJ5T0, num_shells));
  }
  fclose(out_file);

  printf("SP1(m_mu): %g SP2(m_mu): %g SP3(m_mu): %g SP4(m_mu): %g SP5(m_mu): %g\n", sigmap1_tot(y_mu, rhoJ1T0, num_shells), sigmap2_tot(y_mu, rhoJ2T0, num_shells), sigmap3_tot(y_mu, rhoJ3T0, num_shells), sigmap4_tot(y_mu, rhoJ4T0, num_shells), sigmap5_tot(y_mu, rhoJ5T0, num_shells));

  out_file = fopen("SPPJ_q_al27_0_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, sigmapp0_tot(y, rhoJ0T0, num_shells), sigmapp1_tot(y, rhoJ1T0, num_shells), sigmapp2_tot(y, rhoJ2T0, num_shells), sigmapp3_tot(y, rhoJ3T0, num_shells), sigmapp4_tot(y, rhoJ4T0, num_shells), sigmapp5_tot(y, rhoJ5T0, num_shells));
  }
  fclose(out_file);

  printf("SPP0(m_mu): %g SPP1(m_mu): %g SPP2(m_mu): %g SPP3(m_mu): %g SPP4(m_mu): %g SPP5(m_mu): %g\n", sigmapp0_tot(y_mu, rhoJ0T0, num_shells), sigmapp1_tot(y_mu, rhoJ1T0, num_shells), sigmapp2_tot(y_mu, rhoJ2T0, num_shells), sigmapp3_tot(y_mu, rhoJ3T0, num_shells), sigmapp4_tot(y_mu, rhoJ4T0, num_shells), sigmapp5_tot(y_mu, rhoJ5T0, num_shells));


  out_file = fopen("OPJ_q_al27_0_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, omegap0_tot(y, rhoJ0T0, num_shells), omegap1_tot(y, rhoJ1T0, num_shells), omegap2_tot(y, rhoJ2T0, num_shells), omegap3_tot(y, rhoJ3T0, num_shells), 0.0, 0.0);
  }
  fclose(out_file);

  // Compute Al27 1/2+ (1e) -> 5/2+ (gs) transition
  // Only J = 2, 3 operators allowed
  printf("\n");
  for (int i = 0; i < num_shells*num_shells; i++) {
    rhoJ0T0[i] = 0.0;
    rhoJ0T1[i] = 0.0;
    rhoJ1T0[i] = 0.0;
    rhoJ1T1[i] = 0.0;
    rhoJ2T0[i] = 0.0;
    rhoJ2T1[i] = 0.0;
    rhoJ3T0[i] = 0.0;
    rhoJ3T1[i] = 0.0;
    rhoJ4T0[i] = 0.0;
    rhoJ4T1[i] = 0.0;
    rhoJ5T0[i] = 0.0;
    rhoJ5T1[i] = 0.0;
  }

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J2_T0_1_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ2T0[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J2_T1_1_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ2T1[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);
  
  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J3_T0_1_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ3T0[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J3_T1_1_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ3T1[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  out_file = fopen("MJ_q_al27_1_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, m0_tot(y, rhoJ0T0, num_shells), m1_tot(y, rhoJ1T0, num_shells), m2_tot(y, rhoJ2T0, num_shells), m3_tot(y, rhoJ3T0, num_shells), m4_tot(y, rhoJ4T0, num_shells), 0.0);
  }
  fclose(out_file);

  printf("M0(m_mu): %g M1(m_mu): %g M2(m_mu): %g M3(m_mu): %g M4(m_mu): %g\n", m0_tot(y_mu, rhoJ0T0, num_shells), m1_tot(y_mu, rhoJ1T0, num_shells), m2_tot(y_mu, rhoJ2T0, num_shells), m3_tot(y_mu, rhoJ3T0, num_shells), m4_tot(y_mu, rhoJ4T0, num_shells));


  out_file = fopen("DPJ_q_al27_1_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, 0.0, deltap1_tot(y, rhoJ1T0, num_shells), deltap2_tot(y, rhoJ2T0, num_shells), deltap3_tot(y, rhoJ3T0, num_shells), 0.0, 0.0);
  }
  fclose(out_file);

  printf("DP1(m_mu): %g DP2(m_mu): %g DP3(m_mu): %g\n", deltap1_tot(y_mu, rhoJ1T0, num_shells), deltap2_tot(y_mu, rhoJ2T0, num_shells), deltap3_tot(y_mu, rhoJ3T0, num_shells));


  out_file = fopen("SJ_q_al27_1_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, 0.0, sigma1_tot(y, rhoJ1T0, num_shells), sigma2_tot(y, rhoJ2T0, num_shells), sigma3_tot(y, rhoJ3T0, num_shells), sigma4_tot(y, rhoJ4T0, num_shells), 0.0);
  }
  fclose(out_file);

  printf("S1(m_mu): %g S2(m_mu): %g S3(m_mu): %g S4(m_mu): %g\n", sigma1_tot(y_mu, rhoJ1T0, num_shells), sigma2_tot(y_mu, rhoJ2T0, num_shells), sigma3_tot(y_mu, rhoJ3T0, num_shells), sigma4_tot(y_mu, rhoJ4T0, num_shells));


  out_file = fopen("DJ_q_al27_1_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, 0.0, delta1_tot(y, rhoJ1T0, num_shells), delta2_tot(y, rhoJ2T0, num_shells), delta3_tot(y, rhoJ3T0, num_shells), 0.0, 0.0);
  }
  fclose(out_file);

  printf("D1(m_mu): %g D2(m_mu): %g D3(m_mu): %g\n", delta1_tot(y_mu, rhoJ1T0, num_shells), delta2_tot(y_mu, rhoJ2T0, num_shells), delta3_tot(y_mu, rhoJ3T0, num_shells));


  out_file = fopen("SPJ_q_al27_1_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, 0.0, sigmap1_tot(y, rhoJ1T0, num_shells), sigmap2_tot(y, rhoJ2T0, num_shells), sigmap3_tot(y, rhoJ3T0, num_shells), sigmap4_tot(y, rhoJ4T0, num_shells), sigmap5_tot(y, rhoJ5T0, num_shells));
  }
  fclose(out_file);

  printf("SP1(m_mu): %g SP2(m_mu): %g SP3(m_mu): %g SP4(m_mu): %g SP5(m_mu): %g\n", sigmap1_tot(y_mu, rhoJ1T0, num_shells), sigmap2_tot(y_mu, rhoJ2T0, num_shells), sigmap3_tot(y_mu, rhoJ3T0, num_shells), sigmap4_tot(y_mu, rhoJ4T0, num_shells), sigmap5_tot(y_mu, rhoJ5T0, num_shells));

  out_file = fopen("SPPJ_q_al27_1_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, sigmapp0_tot(y, rhoJ0T0, num_shells), sigmapp1_tot(y, rhoJ1T0, num_shells), sigmapp2_tot(y, rhoJ2T0, num_shells), sigmapp3_tot(y, rhoJ3T0, num_shells), sigmapp4_tot(y, rhoJ4T0, num_shells), sigmapp5_tot(y, rhoJ5T0, num_shells));
  }
  fclose(out_file);

  printf("SPP0(m_mu): %g SPP1(m_mu): %g SPP2(m_mu): %g SPP3(m_mu): %g SPP4(m_mu): %g SPP5(m_mu): %g\n", sigmapp0_tot(y_mu, rhoJ0T0, num_shells), sigmapp1_tot(y_mu, rhoJ1T0, num_shells), sigmapp2_tot(y_mu, rhoJ2T0, num_shells), sigmapp3_tot(y_mu, rhoJ3T0, num_shells), sigmapp4_tot(y_mu, rhoJ4T0, num_shells), sigmapp5_tot(y_mu, rhoJ5T0, num_shells));


  out_file = fopen("OPJ_q_al27_1_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, omegap0_tot(y, rhoJ0T0, num_shells), omegap1_tot(y, rhoJ1T0, num_shells), omegap2_tot(y, rhoJ2T0, num_shells), omegap3_tot(y, rhoJ3T0, num_shells), 0.0, 0.0);
  }
  fclose(out_file);

  // Compute Al27 3/2+ (2e) -> 5/2+ (gs) transition
  // Only J = 1,2,3,4 operators allowed
  printf("\n");
  for (int i = 0; i < num_shells*num_shells; i++) {
    rhoJ0T0[i] = 0.0;
    rhoJ0T1[i] = 0.0;
    rhoJ1T0[i] = 0.0;
    rhoJ1T1[i] = 0.0;
    rhoJ2T0[i] = 0.0;
    rhoJ2T1[i] = 0.0;
    rhoJ3T0[i] = 0.0;
    rhoJ3T1[i] = 0.0;
    rhoJ4T0[i] = 0.0;
    rhoJ4T1[i] = 0.0;
    rhoJ5T0[i] = 0.0;
    rhoJ5T1[i] = 0.0;
  }

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J1_T0_2_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ1T0[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J1_T1_2_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ1T1[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);


  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J2_T0_2_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ2T0[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J2_T1_2_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ2T1[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);
  
  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J3_T0_2_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ3T0[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J3_T1_2_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ3T1[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J4_T0_2_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ4T0[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  in_file = fopen("isotope_data/al27/density/al27-al27_usdb_1bdy_J4_T1_2_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ4T1[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);


  out_file = fopen("MJ_q_al27_2_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, m0_tot(y, rhoJ0T0, num_shells), m1_tot(y, rhoJ1T0, num_shells), m2_tot(y, rhoJ2T0, num_shells), m3_tot(y, rhoJ3T0, num_shells), m4_tot(y, rhoJ4T0, num_shells), 0.0);
  }
  fclose(out_file);

  printf("M0(m_mu): %g M1(m_mu): %g M2(m_mu): %g M3(m_mu): %g M4(m_mu): %g\n", m0_tot(y_mu, rhoJ0T0, num_shells), m1_tot(y_mu, rhoJ1T0, num_shells), m2_tot(y_mu, rhoJ2T0, num_shells), m3_tot(y_mu, rhoJ3T0, num_shells), m4_tot(y_mu, rhoJ4T0, num_shells));


  out_file = fopen("DPJ_q_al27_2_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, 0.0, deltap1_tot(y, rhoJ1T0, num_shells), deltap2_tot(y, rhoJ2T0, num_shells), deltap3_tot(y, rhoJ3T0, num_shells), 0.0, 0.0);
  }
  fclose(out_file);

  printf("DP1(m_mu): %g DP2(m_mu): %g DP3(m_mu): %g\n", deltap1_tot(y_mu, rhoJ1T0, num_shells), deltap2_tot(y_mu, rhoJ2T0, num_shells), deltap3_tot(y_mu, rhoJ3T0, num_shells));


  out_file = fopen("SJ_q_al27_2_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, 0.0, sigma1_tot(y, rhoJ1T0, num_shells), sigma2_tot(y, rhoJ2T0, num_shells), sigma3_tot(y, rhoJ3T0, num_shells), sigma4_tot(y, rhoJ4T0, num_shells), 0.0);
  }
  fclose(out_file);

  printf("S1(m_mu): %g S2(m_mu): %g S3(m_mu): %g S4(m_mu): %g\n", sigma1_tot(y_mu, rhoJ1T0, num_shells), sigma2_tot(y_mu, rhoJ2T0, num_shells), sigma3_tot(y_mu, rhoJ3T0, num_shells), sigma4_tot(y_mu, rhoJ4T0, num_shells));


  out_file = fopen("DJ_q_al27_2_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, 0.0, delta1_tot(y, rhoJ1T0, num_shells), delta2_tot(y, rhoJ2T0, num_shells), delta3_tot(y, rhoJ3T0, num_shells), 0.0, 0.0);
  }
  fclose(out_file);

  printf("D1(m_mu): %g D2(m_mu): %g D3(m_mu): %g\n", delta1_tot(y_mu, rhoJ1T0, num_shells), delta2_tot(y_mu, rhoJ2T0, num_shells), delta3_tot(y_mu, rhoJ3T0, num_shells));


  out_file = fopen("SPJ_q_al27_2_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, 0.0, sigmap1_tot(y, rhoJ1T0, num_shells), sigmap2_tot(y, rhoJ2T0, num_shells), sigmap3_tot(y, rhoJ3T0, num_shells), sigmap4_tot(y, rhoJ4T0, num_shells), sigmap5_tot(y, rhoJ5T0, num_shells));
  }
  fclose(out_file);

  printf("SP1(m_mu): %g SP2(m_mu): %g SP3(m_mu): %g SP4(m_mu): %g SP5(m_mu): %g\n", sigmap1_tot(y_mu, rhoJ1T0, num_shells), sigmap2_tot(y_mu, rhoJ2T0, num_shells), sigmap3_tot(y_mu, rhoJ3T0, num_shells), sigmap4_tot(y_mu, rhoJ4T0, num_shells), sigmap5_tot(y_mu, rhoJ5T0, num_shells));

  out_file = fopen("SPPJ_q_al27_2_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, sigmapp0_tot(y, rhoJ0T0, num_shells), sigmapp1_tot(y, rhoJ1T0, num_shells), sigmapp2_tot(y, rhoJ2T0, num_shells), sigmapp3_tot(y, rhoJ3T0, num_shells), sigmapp4_tot(y, rhoJ4T0, num_shells), sigmapp5_tot(y, rhoJ5T0, num_shells));
  }
  fclose(out_file);

  printf("SPP0(m_mu): %g SPP1(m_mu): %g SPP2(m_mu): %g SPP3(m_mu): %g SPP4(m_mu): %g SPP5(m_mu): %g\n", sigmapp0_tot(y_mu, rhoJ0T0, num_shells), sigmapp1_tot(y_mu, rhoJ1T0, num_shells), sigmapp2_tot(y_mu, rhoJ2T0, num_shells), sigmapp3_tot(y_mu, rhoJ3T0, num_shells), sigmapp4_tot(y_mu, rhoJ4T0, num_shells), sigmapp5_tot(y_mu, rhoJ5T0, num_shells));


  out_file = fopen("OPJ_q_al27_2_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, omegap0_tot(y, rhoJ0T0, num_shells), omegap1_tot(y, rhoJ1T0, num_shells), omegap2_tot(y, rhoJ2T0, num_shells), omegap3_tot(y, rhoJ3T0, num_shells), 0.0, 0.0);
  }
  fclose(out_file);


  // Compute Si28 0+ (gs) -> 0+ (gs) transition
  // Only J = 1,2,3,4 operators allowed
  printf("\n");
  for (int i = 0; i < num_shells*num_shells; i++) {
    rhoJ0T0[i] = 0.0;
    rhoJ0T1[i] = 0.0;
    rhoJ1T0[i] = 0.0;
    rhoJ1T1[i] = 0.0;
    rhoJ2T0[i] = 0.0;
    rhoJ2T1[i] = 0.0;
    rhoJ3T0[i] = 0.0;
    rhoJ3T1[i] = 0.0;
    rhoJ4T0[i] = 0.0;
    rhoJ4T1[i] = 0.0;
    rhoJ5T0[i] = 0.0;
    rhoJ5T1[i] = 0.0;
  }

  in_file = fopen("isotope_data/si28/density/si28-si28_core_1bdy_J0_T0_0_0.dens", "r");
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &inp, &ijp, &in, &ij, &density) == 5) {
    int ind_i = get_shell_index(in, ij);
    int ind_f = get_shell_index(inp, ijp);
    rhoJ0T0[ind_f + num_shells*ind_i] = density;
  }
  fclose(in_file);

  out_file = fopen("MJ_q_si28_0_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, m0_tot(y, rhoJ0T0, num_shells), m1_tot(y, rhoJ1T0, num_shells), m2_tot(y, rhoJ2T0, num_shells), m3_tot(y, rhoJ3T0, num_shells), m4_tot(y, rhoJ4T0, num_shells), 0.0);
  }
  fclose(out_file);

  printf("M0(m_mu): %g M1(m_mu): %g M2(m_mu): %g M3(m_mu): %g M4(m_mu): %g\n", m0_tot(y_mu, rhoJ0T0, num_shells), m1_tot(y_mu, rhoJ1T0, num_shells), m2_tot(y_mu, rhoJ2T0, num_shells), m3_tot(y_mu, rhoJ3T0, num_shells), m4_tot(y_mu, rhoJ4T0, num_shells));


  out_file = fopen("DPJ_q_si28_0_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, 0.0, deltap1_tot(y, rhoJ1T0, num_shells), deltap2_tot(y, rhoJ2T0, num_shells), deltap3_tot(y, rhoJ3T0, num_shells), 0.0, 0.0);
  }
  fclose(out_file);

  printf("DP1(m_mu): %g DP2(m_mu): %g DP3(m_mu): %g\n", deltap1_tot(y_mu, rhoJ1T0, num_shells), deltap2_tot(y_mu, rhoJ2T0, num_shells), deltap3_tot(y_mu, rhoJ3T0, num_shells));


  out_file = fopen("SJ_q_si28_0_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, 0.0, sigma1_tot(y, rhoJ1T0, num_shells), sigma2_tot(y, rhoJ2T0, num_shells), sigma3_tot(y, rhoJ3T0, num_shells), sigma4_tot(y, rhoJ4T0, num_shells), 0.0);
  }
  fclose(out_file);

  printf("S1(m_mu): %g S2(m_mu): %g S3(m_mu): %g S4(m_mu): %g\n", sigma1_tot(y_mu, rhoJ1T0, num_shells), sigma2_tot(y_mu, rhoJ2T0, num_shells), sigma3_tot(y_mu, rhoJ3T0, num_shells), sigma4_tot(y_mu, rhoJ4T0, num_shells));


  out_file = fopen("DJ_q_si28_0_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, 0.0, delta1_tot(y, rhoJ1T0, num_shells), delta2_tot(y, rhoJ2T0, num_shells), delta3_tot(y, rhoJ3T0, num_shells), 0.0, 0.0);
  }
  fclose(out_file);

  printf("D1(m_mu): %g D2(m_mu): %g D3(m_mu): %g\n", delta1_tot(y_mu, rhoJ1T0, num_shells), delta2_tot(y_mu, rhoJ2T0, num_shells), delta3_tot(y_mu, rhoJ3T0, num_shells));


  out_file = fopen("SPJ_q_si28_0_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, 0.0, sigmap1_tot(y, rhoJ1T0, num_shells), sigmap2_tot(y, rhoJ2T0, num_shells), sigmap3_tot(y, rhoJ3T0, num_shells), sigmap4_tot(y, rhoJ4T0, num_shells), sigmap5_tot(y, rhoJ5T0, num_shells));
  }
  fclose(out_file);

  printf("SP1(m_mu): %g SP2(m_mu): %g SP3(m_mu): %g SP4(m_mu): %g SP5(m_mu): %g\n", sigmap1_tot(y_mu, rhoJ1T0, num_shells), sigmap2_tot(y_mu, rhoJ2T0, num_shells), sigmap3_tot(y_mu, rhoJ3T0, num_shells), sigmap4_tot(y_mu, rhoJ4T0, num_shells), sigmap5_tot(y_mu, rhoJ5T0, num_shells));

  out_file = fopen("SPPJ_q_si28_0_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, sigmapp0_tot(y, rhoJ0T0, num_shells), sigmapp1_tot(y, rhoJ1T0, num_shells), sigmapp2_tot(y, rhoJ2T0, num_shells), sigmapp3_tot(y, rhoJ3T0, num_shells), sigmapp4_tot(y, rhoJ4T0, num_shells), sigmapp5_tot(y, rhoJ5T0, num_shells));
  }
  fclose(out_file);

  printf("SPP0(m_mu): %g SPP1(m_mu): %g SPP2(m_mu): %g SPP3(m_mu): %g SPP4(m_mu): %g SPP5(m_mu): %g\n", sigmapp0_tot(y_mu, rhoJ0T0, num_shells), sigmapp1_tot(y_mu, rhoJ1T0, num_shells), sigmapp2_tot(y_mu, rhoJ2T0, num_shells), sigmapp3_tot(y_mu, rhoJ3T0, num_shells), sigmapp4_tot(y_mu, rhoJ4T0, num_shells), sigmapp5_tot(y_mu, rhoJ5T0, num_shells));


  out_file = fopen("OPJ_q_si28_0_0.dat", "w");
  for (int i = 0; i < 150; i++) {
    double y = pow(i*b_osc(27)/(HBARC*2.0), 2.0);
    fprintf(out_file, "%g, %g, %g, %g, %g, %g, %g\n", i*1.0, omegap0_tot(y, rhoJ0T0, num_shells), omegap1_tot(y, rhoJ1T0, num_shells), omegap2_tot(y, rhoJ2T0, num_shells), omegap3_tot(y, rhoJ3T0, num_shells), 0.0, 0.0);
  }
  fclose(out_file);


  return;
}

double omegap3_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(2, 1);
  int ind_f = get_shell_index(2, 5);
  mat += omegap3_1d5_2s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 5);
  mat += omegap3_1d5_1d3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  return mat;
}

double omegap2_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(1, 3);
  int ind_f = get_shell_index(2, 1);
  mat += omegap2_2s1_1p3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 3);
  mat += omegap2_1d3_1p1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 3);
  mat += omegap2_1d3_1p3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 5);
  mat += omegap2_1d5_1p1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  return mat;
}

double omegap1_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(1, 1);
  int ind_f = get_shell_index(1, 3);
  mat += omegap1_1p3_1p1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(2, 1);
  mat += omegap1_2s1_1s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(2, 3);
  mat += omegap1_1d3_1s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 1);
  ind_f = get_shell_index(2, 3);
  mat += omegap1_1d3_2s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 5);
  mat += omegap1_1d5_1d3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  return mat;
}

double omegap0_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(0, 1);
  int ind_f = get_shell_index(1, 1);
  mat += omegap0_1p1_1s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 1);
  mat += omegap0_2s1_1p1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 3);
  mat += omegap0_1d3_1p3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  return mat;
}


double sigmapp5_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(2, 5);
  int ind_f = get_shell_index(2, 5);
  mat += sigmapp5_1d5_1d5(y)*rho[ind_f + num_shells*ind_i];

  return mat;
}


double sigmapp4_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(1, 3);
  int ind_f = get_shell_index(2, 5);
  mat += sigmapp4_1d5_1p3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  return mat;
}

double sigmapp3_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(1, 3);
  int ind_f = get_shell_index(1, 3);
  mat += sigmapp3_1p3_1p3(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(2, 5);
  mat += sigmapp3_1d5_1s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 3);
  mat += sigmapp3_1d3_1d3(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(2, 1);
  ind_f = get_shell_index(2, 5);
  mat += sigmapp3_1d5_2s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 5);
  mat += sigmapp3_1d5_1d3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 5);
  ind_f = get_shell_index(2, 5);
  mat += sigmapp3_1d5_1d5(y)*rho[ind_f + num_shells*ind_i];

  return mat;
}

double sigmapp2_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(0, 1);
  int ind_f = get_shell_index(1, 3);
  mat += sigmapp2_1p3_1s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 1);
  mat += sigmapp2_2s1_1p3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 3);
  mat += sigmapp2_1d3_1p1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 3);
  mat += sigmapp2_1d3_1p3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 5);
  mat += sigmapp2_1d5_1p1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 5);
  mat += sigmapp2_1d5_1p3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  return mat;
}

double sigmapp1_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(0, 1);
  int ind_f = get_shell_index(0, 1);
  mat += sigmapp1_1s1_1s1(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(1, 1);
  mat += sigmapp1_1p1_1p1(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(1, 3);
  mat += sigmapp1_1p3_1p1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(1, 3);
  mat += sigmapp1_1p3_1p3(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(2, 1);
  mat += sigmapp1_2s1_1s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(2, 3);
  mat += sigmapp1_1d3_1s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 1);
  ind_f = get_shell_index(2, 1);
  mat += sigmapp1_2s1_2s1(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(2, 1);
  ind_f = get_shell_index(2, 3);
  mat += sigmapp1_1d3_2s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 3);
  mat += sigmapp1_1d3_1d3(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 5);
  mat += sigmapp1_1d5_1d3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 5);
  ind_f = get_shell_index(2, 5);
  mat += sigmapp1_1d5_1d5(y)*rho[ind_f + num_shells*ind_i];

  return mat;
}

double sigmapp0_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(0, 1);
  int ind_f = get_shell_index(1, 1);
  mat += sigmapp0_1p1_1s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 1);
  mat += sigmapp0_2s1_1p1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 3);
  mat += sigmapp0_1d3_1p3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  return mat;
}

double sigmap5_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(2, 5);
  int ind_f = get_shell_index(2, 5);
  mat += sigmap5_1d5_1d5(y)*rho[ind_f + num_shells*ind_i];

  return mat;
}


double sigmap4_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(1, 3);
  int ind_f = get_shell_index(2, 5);
  mat += sigmap4_1d5_1p3(y)*rho[ind_f + num_shells*ind_i];

  return mat;
}

double sigmap3_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(1, 3);
  int ind_f = get_shell_index(1, 3);
  mat += sigmap3_1p3_1p3(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(2, 5);
  mat += sigmap3_1d5_1s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 3);
  mat += sigmap3_1d3_1d3(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(2, 1);
  ind_f = get_shell_index(2, 5);
  mat += sigmap3_1d5_2s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 5);
  mat += sigmap3_1d5_1d3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 5);
  ind_f = get_shell_index(2, 5);
  mat += sigmap3_1d5_1d5(y)*rho[ind_f + num_shells*ind_i];
 
  return mat;
}

double sigmap2_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(0, 1);
  int ind_f = get_shell_index(1, 3);
  mat += sigmap2_1p3_1s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 1);
  mat += sigmap2_2s1_1p3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 3);
  mat += sigmap2_1d3_1p1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 3);
  mat += sigmap2_1d3_1p3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 5);
  mat += sigmap2_1d5_1p1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 5);
  mat += sigmap2_1d5_1p3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  return mat;
}

double sigmap1_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(0, 1);
  int ind_f = get_shell_index(0, 1);
  mat += sigmap1_1s1_1s1(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(1, 1);
  mat += sigmap1_1p1_1p1(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(1, 3);
  mat += sigmap1_1p3_1p1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(1, 3);
  mat += sigmap1_1p3_1p3(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(2, 1);
  mat += sigmap1_2s1_1s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(2, 3);
  mat += sigmap1_1d3_1s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 1);
  ind_f = get_shell_index(2, 1);
  mat += sigmap1_2s1_2s1(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(2, 1);
  ind_f = get_shell_index(2, 3);
  mat += sigmap1_1d3_2s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 3);
  mat += sigmap1_1d3_1d3(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 5);
  mat += sigmap1_1d5_1d3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 5);
  ind_f = get_shell_index(2, 5);
  mat += sigmap1_1d5_1d5(y)*rho[ind_f + num_shells*ind_i];

  return mat;
}

double delta3_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(2, 3);
  int ind_f = get_shell_index(2, 3);
  mat += delta3_1d3_1d3(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 5);
  mat += delta3_1d5_1d3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 5);
  ind_f = get_shell_index(2, 5);
  mat += delta3_1d5_1d5(y)*rho[ind_f + num_shells*ind_i];

  return mat;
}

double delta2_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(1, 1);
  int ind_f = get_shell_index(2, 3);
  mat += delta2_1d3_1p1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 3);
  mat += delta2_1d3_1p3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 5);
  mat += delta2_1d5_1p1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 5);
  mat += delta2_1d5_1p3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  return mat;
}

double delta1_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(1, 1);
  int ind_f = get_shell_index(1, 1);
  mat += delta1_1p1_1p1(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(1, 3);
  mat += delta1_1p3_1p1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(1, 3);
  mat += delta1_1p3_1p3(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 3);
  mat += delta1_1d3_1d3(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 5);
  mat += delta1_1d5_1d3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 5);
  ind_f = get_shell_index(2, 5);
  mat += delta1_1d5_1d5(y)*rho[ind_f + num_shells*ind_i];

  return mat;
}

double sigma4_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(2, 3);
  int ind_f = get_shell_index(2, 5);
  mat += sigma4_1d5_1d3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);
  
  return mat;
}

double sigma3_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(1, 3);
  int ind_f = get_shell_index(2, 3);
  mat += sigma3_1d3_1p3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 5);
  mat += sigma3_1d5_1p1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 5);
  mat += sigma3_1d5_1p3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  return mat;
}

double sigma2_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(1, 1);
  int ind_f = get_shell_index(1, 3);
  mat += sigma2_1p3_1p1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(2, 3);
  mat += sigma2_1d3_1s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(2, 5);
  mat += sigma2_1d5_1s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 1);
  ind_f = get_shell_index(2, 3);
  mat += sigma2_1d3_2s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 1);
  ind_f = get_shell_index(2, 5);
  mat += sigma2_1d5_2s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 5);
  mat += sigma2_1d5_1d3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  return mat;
}

double sigma1_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(0, 1);
  int ind_f = get_shell_index(1, 1);
  mat += sigma1_1p1_1s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(1, 3);
  mat += sigma1_1p3_1s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 1);
  mat += sigma1_2s1_1p1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 1);
  mat += sigma1_2s1_1p3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 3);
  mat += sigma1_1d3_1p1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 3);
  mat += sigma1_1d3_1p3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 5);
  mat += sigma1_1d5_1p3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  return mat;
}

double deltap3_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(1, 3);
  int ind_f = get_shell_index(2, 3);
  mat += deltap3_1d3_1p3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 5);
  mat += deltap3_1d5_1p1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 5);
  mat += deltap3_1d5_1p3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  return mat;
}

double deltap2_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(0, 1);
  int ind_f = get_shell_index(2, 3);
  mat += deltap2_1d3_1s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(2, 5);
  mat += deltap2_1d5_1s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 1);
  ind_f = get_shell_index(2, 3);
  mat += deltap2_1d3_2s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 1);
  ind_f = get_shell_index(2, 5);
  mat += deltap2_1d5_2s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  return mat;
}

double deltap1_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(0, 1);
  int ind_f = get_shell_index(1, 1);
  mat += deltap1_1p1_1s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(1, 3);
  mat += deltap1_1p3_1s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);
  
  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 1);
  mat += deltap1_2s1_1p1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 1);
  mat += deltap1_2s1_1p3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 3);
  mat += deltap1_1d3_1p1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 3);
  mat += deltap1_1d3_1p3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 5);
  mat += deltap1_1d5_1p3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);
 
  return mat;
}

double m4_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(2, 3);
  int ind_f = get_shell_index(2, 5);
  mat += m4_1d5_1d3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 5);
  ind_f = get_shell_index(2, 5);
  mat += m4_1d5_1d5(y)*rho[ind_f + num_shells*ind_i];

  return mat;
}


double m3_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(1, 3);
  int ind_f = get_shell_index(2, 3);
  mat += m3_1d3_1p3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 5);
  mat += m3_1d5_1p1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 5);
  mat += m3_1d5_1p3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  return mat;
}
 

double m2_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(1, 1);
  int ind_f = get_shell_index(1, 3);
  mat += m2_1p3_1p1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(1, 3);
  mat += m2_1p3_1p3(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(2, 3);
  mat += m2_1d3_1s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(2, 5);
  mat += m2_1d5_1s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 1);
  ind_f = get_shell_index(2, 3);
  mat += m2_1d3_2s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 3);
  mat += m2_1d3_1d3(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(2, 1);
  ind_f = get_shell_index(2, 5);
  mat += m2_1d5_2s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 5);
  mat += m2_1d5_1d3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 5);
  ind_f = get_shell_index(2, 5);
  mat += m2_1d5_1d5(y)*rho[ind_f + num_shells*ind_i];

  return mat;
}

double m1_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(0, 1);
  int ind_f = get_shell_index(1, 1);
  mat += m1_1p1_1s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(1, 3);
  mat += m1_1p3_1s1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 1);
  mat += m1_2s1_1p1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 1);
  mat += m1_2s1_1p3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(2, 3);
  mat += m1_1d3_1p1(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 3);
  mat += m1_1d3_1p3(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(2, 5);
  mat += m1_1d5_1p3(y)*(rho[ind_f + num_shells*ind_i] - rho[ind_i + num_shells*ind_f]);

  return mat;
}

double m0_tot(double y, double *rho, int num_shells) {
  double mat = 0.0; 

  int ind_i = get_shell_index(0, 1);
  int ind_f = get_shell_index(0, 1);
  mat += m0_1s1_1s1(y)*rho[ind_f + num_shells*ind_i];
  
  ind_i = get_shell_index(1, 1);
  ind_f = get_shell_index(1, 1);
  mat += m0_1p1_1p1(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(1, 3);
  ind_f = get_shell_index(1, 3);
  mat += m0_1p3_1p3(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(0, 1);
  ind_f = get_shell_index(2, 1);
  mat += m0_2s1_1s1(y)*(rho[ind_f + num_shells*ind_i] + rho[ind_i + num_shells*ind_f]);

  ind_i = get_shell_index(2, 1);
  ind_f = get_shell_index(2, 1);
  mat += m0_2s1_2s1(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(2, 3);
  ind_f = get_shell_index(2, 3);
  mat += m0_1d3_1d3(y)*rho[ind_f + num_shells*ind_i];

  ind_i = get_shell_index(2, 5);
  ind_f = get_shell_index(2, 5);
  mat += m0_1d5_1d5(y)*rho[ind_f + num_shells*ind_i];
 
  return mat;
}
 

double m0_1s1_1s1(double y) {
  // Verified working
  double m0 = 1.0/sqrt(2.0*M_PI)*exp(-y);

  return m0;
}

double sigmap1_1s1_1s1(double y) {
  // Verified working
  double s1 = 2.0/sqrt(4.0*M_PI)*exp(-y);

  return s1;
}

double sigmapp1_1s1_1s1(double y) {
  // Verified working
  double s1 = 1.0/sqrt(2.0*M_PI)*exp(-y);

  return s1;
}

double m1_1p1_1s1(double y) {
  // Verified working
  double m1 = -1.0/sqrt(3.0*M_PI)*sqrt(y)*exp(-y);

  return m1;
}

double deltap1_1p1_1s1(double y) {
  // Verified working
  double d1 = 1.0/sqrt(6.0*4.0*M_PI)*exp(-y)/sqrt(y);

  return d1;
}

double sigma1_1p1_1s1(double y) {
  // Verified working
  double s1 = 2.0/sqrt(6.0*M_PI)*sqrt(y)*exp(-y);

  return s1;
}

double m1_1p3_1s1(double y) {
  // Verified working
  double m1 = 2.0/sqrt(6.0*M_PI)*sqrt(y)*exp(-y);

  return m1;
}

double deltap1_1p3_1s1(double y) {
  // Verified working
  double d1 = -1.0/sqrt(12.0*M_PI)*exp(-y)/sqrt(y);

  return d1;
}

double sigma1_1p3_1s1(double y) {
  // Verified working
  double s1 = 1.0/sqrt(3.0*M_PI)*exp(-y)*sqrt(y);

  return s1;
}

double sigmapp0_1p1_1s1(double y) {
  // Verified working
  double s0 = 1.0/sqrt(3.0*M_PI)*exp(-y)*sqrt(y);

  return s0;
}

double omegap0_1p1_1s1(double y) {
  // Verified working
  double o0 = 1.0/4.0*sqrt(3.0/M_PI)*exp(-y)/sqrt(y);

  return o0;
}

double sigmap2_1p3_1s1(double y) {
  // Verified working
  double s2 = 1.0/sqrt(M_PI)*exp(-y)*sqrt(y);

  return s2;
}

double sigmapp2_1p3_1s1(double y) {
  // Verified working
  double s2 = 2.0/sqrt(6.0*M_PI)*exp(-y)*sqrt(y);

  return s2;
}

double m0_1p1_1p1(double y) {
  // Verified working
  double m0 = 1.0/sqrt(2.0*M_PI)*exp(-y)*(1.0 - 2.0/3.0*y);

  return m0;
}

double m2_1p3_1p1(double y) {
  // Verified working
  double m2 = -2.0/(3.0*sqrt(M_PI))*y*exp(-y);

  return m2;
}

double sigma2_1p3_1p1(double y) {
  // Verified working
  double s2 = -2.0/sqrt(6.0*M_PI)*exp(-y)*y;

  return s2;
}

double m0_1p3_1p3(double y) {
  // Verified working
  double m0 = 1.0/sqrt(M_PI)*exp(-y)*(1.0 - 2.0/3.0*y);

  return m0;
}

double m2_1p3_1p3(double y) {
  // Verified working
  double m2 = -2.0/(3.0*sqrt(M_PI))*exp(-y)*y;

  return m2;
}

double delta1_1p1_1p1(double y) {
  // Verified working
  double d1 = -1.0/(3.0*sqrt(M_PI))*exp(-y);

  return d1;
}

double sigmap1_1p1_1p1(double y) {
  // Verified working
  double s1 = 1.0/(3.0*sqrt(M_PI))*exp(-y)*(-1.0 + 2.0*y);

  return s1;
}

double sigmapp1_1p1_1p1(double y) {
  // Verified working
  double s1 = -1.0/(3.0*sqrt(2.0*M_PI))*exp(-y)*(1.0 + 2.0*y);

  return s1;
}

double delta1_1p3_1p1(double y) {
  // Verified working
  double d1 = -1.0/(3.0*sqrt(2.0*M_PI))*exp(-y);

  return d1;
}

double sigmap1_1p3_1p1(double y) {
  // Verified working
  double s1 = 2.0*sqrt(2.0)/(3.0*sqrt(M_PI))*exp(-y)*(-1.0 + y/2.0);

  return s1;
}

double sigmapp1_1p3_1p1(double y) {
  // Verified working
  double s1 = 2.0/(3.0*sqrt(M_PI))*exp(-y)*(-1.0 + y);

  return s1;
}

double omegap1_1p3_1p1(double y) {
  // Verified working
  double o1 = -1.0/sqrt(4.0*M_PI)*exp(-y);

  return o1;
}

double delta1_1p3_1p3(double y) {
  // Verified working
  double d1 = -sqrt(10.0/M_PI)/6.0*exp(-y);

  return d1;
}

double sigmap1_1p3_1p3(double y) {
  // Verified working
  double s1 = sqrt(10.0/M_PI)/3.0*exp(-y)*(1.0 - 4.0/5.0*y);

  return s1;
}

double sigmapp1_1p3_1p3(double y) {
  // Verified working
  double s1 = sqrt(5.0/M_PI)/3.0*exp(-y)*(1.0 - 2.0/5.0*y);

  return s1;
}

double sigmap3_1p3_1p3(double y) {
  // Verified working
  double s3 = -4.0/sqrt(15.0*M_PI)*exp(-y)*y;

  return s3;
}

double sigmapp3_1p3_1p3(double y) {
  // Verified working
  double s3 = -2.0/sqrt(5.0*M_PI)*exp(-y)*y;

  return s3;
}

double m0_2s1_1s1(double y) {
  // Verified working
  double m0 = 1.0/sqrt(3.0*M_PI)*exp(-y)*y;

  return m0;
}

double m2_1d3_1s1(double y) {
  // Verified working
  double m2 = -2.0/sqrt(15.0*M_PI)*exp(-y)*y;

  return m2;
}

double deltap2_1d3_1s1(double y) {
  // Verified working
  double d2 = 1.0/sqrt(10.0*M_PI)*exp(-y);

  return d2;
}

double sigma2_1d3_1s1(double y) {
  // Verified working
  double s2 = 2.0/sqrt(10.0*M_PI)*exp(-y)*y;

  return s2;
}

double m2_1d5_1s1(double y) {
  // Verified working
  double m2 = 2.0/sqrt(10.0*M_PI)*exp(-y)*y;

  return m2;
}

double deltap2_1d5_1s1(double y) {
  // Verified working
  double d2 = -3.0/(2.0*sqrt(15.0*M_PI))*exp(-y);

  return d2;
}

double sigma2_1d5_1s1(double y) {
  // Verified working
  double s2 = 2.0/sqrt(15.0*M_PI)*exp(-y)*y;

  return s2;
}

double sigmap1_2s1_1s1(double y) {
  // Verified working
  double s1 = 2.0/sqrt(6.0*M_PI)*exp(-y)*y;

  return s1;
}

double sigmapp1_2s1_1s1(double y) {
  // Verified working
  double s1 = 1.0/sqrt(3.0*M_PI)*exp(-y)*y;

  return s1;
}

double omegap1_2s1_1s1(double y) {
  // Verified working
  double o1 = 1.0/(2.0*sqrt(3.0*M_PI))*exp(-y);

  return o1;
}

double sigmap1_1d3_1s1(double y) {
  // Verified working
  double s1 = -2.0/sqrt(30.0*M_PI)*exp(-y)*y;

  return s1;
}

double sigmapp1_1d3_1s1(double y) {
  // Verified working 
  double s1 = 2.0/sqrt(15.0*M_PI)*exp(-y)*y;

  return s1;
}

double omegap1_1d3_1s1(double y) {
  // Verified working 
  double o1 = 5.0/(2.0*sqrt(15.0*M_PI))*exp(-y);

  return o1;
}

double sigmap3_1d5_1s1(double y) {
  // Verified working 
  double s3 = 4.0/sqrt(30.0*M_PI)*exp(-y)*y;

  return s3;
}

double sigmapp3_1d5_1s1(double y) {
  // Verified working
  double s3 = 2.0/sqrt(10.0*M_PI)*exp(-y)*y;

  return s3;
}

double m1_2s1_1p1(double y) {
  // Verified working
  double m1 = 2.0/(3.0*sqrt(2.0*M_PI))*exp(-y)*(sqrt(y) - pow(y, 3.0/2.0));

  return m1;
}

double deltap1_2s1_1p1(double y) {
  // Verified working
  double d1 = -1.0/(3.0*sqrt(4.0*M_PI))*exp(-y)*(1.0/sqrt(y) + sqrt(y));

  return d1;
}

double sigma1_2s1_1p1(double y) {
  // Verified working
  double s1 = 2.0/(3.0*sqrt(M_PI))*exp(-y)*(sqrt(y) - pow(y, 3.0/2.0));

  return s1;
}

double m1_2s1_1p3(double y) {
  // Verified working
  double m1 = 2.0/(3.0*sqrt(M_PI))*exp(-y)*(sqrt(y) - pow(y, 3.0/2.0));

  return m1;
}

double deltap1_2s1_1p3(double y) {
  // Verified working
  double d1 = -1.0/(3.0*sqrt(2.0*M_PI))*exp(-y)*(1.0/sqrt(y) + sqrt(y));

  return d1;
}

double sigma1_2s1_1p3(double y) {
  // Verified working 
  double s1 = 1.0/3.0*sqrt(2.0/M_PI)*exp(-y)*(-sqrt(y) + pow(y, 3.0/2.0));

  return s1;
}

double m1_1d3_1p1(double y) {
  // Verified working
  double m1 = 1.0/3.0*sqrt(10.0/M_PI)*exp(-y)*(sqrt(y) - 2.0/5.0*pow(y, 3.0/2.0));

  return m1;
}

double deltap1_1d3_1p1(double y) {
  // Verified working
  double d1 = 1.0/6.0*sqrt(5.0/M_PI)*exp(-y)*(-1.0/sqrt(y) + 4.0/5.0*sqrt(y));

  return d1;
}

double sigma1_1d3_1p1(double y) {
  // Verified working
  double s1 = 1.0/3.0*sqrt(5.0/M_PI)*exp(-y)*(-sqrt(y) + 2.0/5.0*pow(y, 3.0/2.0));

  return s1;
}

double m1_1d3_1p3(double y) {
  // Verified working
  double m1 = 1.0/3.0*sqrt(2.0/M_PI)*exp(-y)*(-sqrt(y) + 2.0/5.0*pow(y, 3.0/2.0));

  return m1;
}

double deltap1_1d3_1p3(double y) {
  // Verified working
  double d1 = 1.0/(6.0*sqrt(M_PI))*exp(-y)*(1.0/sqrt(y) - 4.0/5.0*sqrt(y));

  return d1;
}

double sigma1_1d3_1p3(double y) {
  // Verified working 
  double s1 = 4.0/(3.0*sqrt(M_PI))*exp(-y)*(sqrt(y) - 2.0/5.0*pow(y, 3.0/2.0));
  
  return s1;
}

double m3_1d3_1p3(double y) {
  // Verified working
  double m3 = 2.0/5.0*sqrt(2.0/M_PI)*exp(-y)*pow(y, 3.0/2.0);

  return m3;
}

double deltap3_1d3_1p3(double y) {
  // Verified working
  double d3 = -1.0/15.0*sqrt(6.0/M_PI)*exp(-y)*sqrt(y);

  return d3;
}

double sigma3_1d3_1p3(double y) {
  // Verified working
  double s3 = -4.0/15.0*sqrt(6.0/M_PI)*exp(-y)*pow(y, 3.0/2.0);

  return s3;
}

double m3_1d5_1p1(double y) {
  // Verified working
  double m3 = -2.0/sqrt(15.0*M_PI)*exp(-y)*pow(y, 3.0/2.0);

  return m3;
}

double deltap3_1d5_1p1(double y) {
  // Verified working
  double d3 = 1.0/15.0*sqrt(5.0/M_PI)*exp(-y)*sqrt(y);

  return d3;
}

double sigma3_1d5_1p1(double y) {
  // Verified working
  double s3 = -4.0/15.0*sqrt(5.0/M_PI)*exp(-y)*pow(y, 3.0/2.0);

  return s3;
}

double m1_1d5_1p3(double y) {
  // Verified working
  double m1 = sqrt(2.0/M_PI)*exp(-y)*(sqrt(y) - 2.0/5.0*pow(y, 3.0/2.0));

  return m1;
}

double deltap1_1d5_1p3(double y) {
  // Verified working 
  double d1 = 1.0/sqrt(4.0*M_PI)*exp(-y)*(-1.0/sqrt(y) + 4.0/5.0*sqrt(y));

  return d1;
}

double sigma1_1d5_1p3(double y) {
  // Verified working
  double s1 = 1.0/sqrt(M_PI)*exp(-y)*(sqrt(y) - 2.0/5.0*pow(y, 3.0/2.0));

  return s1;
}

double m3_1d5_1p3(double y) {
  // Verified working
  double m3 = -4.0/15.0*sqrt(3.0/M_PI)*exp(-y)*pow(y, 3.0/2.0);

  return m3;
}

double deltap3_1d5_1p3(double y) {
  // Verified working
  double d3 = 2.0/(15.0*sqrt(M_PI))*exp(-y)*sqrt(y);

  return d3;
}

double sigma3_1d5_1p3(double y) {
  // Verified working
  double s3 = -2.0/(15.0*sqrt(M_PI))*exp(-y)*pow(y, 3.0/2.0);

  return s3;
}

double sigmapp0_2s1_1p1(double y) {
  // Verified working
  double s0 = 1.0/3.0*sqrt(2.0/M_PI)*exp(-y)*(-sqrt(y) + pow(y, 3.0/2.0));

  return s0;
}

double omegap0_2s1_1p1(double y) {
  // Verified working
  double o0 = -1.0/4.0*sqrt(2.0/M_PI)*exp(-y)*(1.0/sqrt(y) + 1.0/3.0*sqrt(y));

  return o0;
}

double sigmap2_2s1_1p3(double y) {
  // Verified working
  double s2 = 1.0/3.0*sqrt(6.0/M_PI)*exp(-y)*(sqrt(y) - pow(y, 3.0/2.0));

  return s2;
}

double sigmapp2_2s1_1p3(double y) {
  // Verified working
  double s2 = 2.0/(3.0*sqrt(M_PI))*exp(-y)*(sqrt(y) - pow(y, 3.0/2.0));

  return s2;
}

double omegap2_2s1_1p3(double y) {
  // Verified working
  double o2 = -1.0/(3.0*sqrt(M_PI))*exp(-y)*sqrt(y);

  return o2;
}

double delta2_1d3_1p1(double y) {
  // Verified working
  double d2 = -1.0/sqrt(15.0*M_PI)*sqrt(y)*exp(-y);

  return d2;
}

double sigmap2_1d3_1p1(double y) {
  // Verified working
  double s2 = 1.0/sqrt(15.0*M_PI)*exp(-y)*(-sqrt(y) + 2.0*pow(y, 3.0/2.0));

  return s2;
}

double sigmapp2_1d3_1p1(double y) {
  // Verified working
  double s2 = -1.0/15.0*sqrt(10.0/M_PI)*exp(-y)*(sqrt(y) + 2.0*pow(y, 3.0/2.0));

  return s2;
}

double omegap2_1d3_1p1(double y) {
  // Verified working
  double o2 = -1.0/15.0*sqrt(10.0/M_PI)*exp(-y)*sqrt(y);

  return o2;
}

double sigmapp0_1d3_1p3(double y) {
  // Verified working
  double s0 = 1.0/3.0*sqrt(10.0/M_PI)*exp(-y)*(sqrt(y) - 2.0/5.0*pow(y, 3.0/2.0));

  return s0;
}

double omegap0_1d3_1p3(double y) {
  // Verified working
  double o0 = 1.0/4.0*sqrt(10.0/M_PI)*exp(-y)*(1.0/sqrt(y) - 2.0/3.0*sqrt(y));

  return o0;
}

double delta2_1d3_1p3(double y) {
  // Verified working
  double d2 = 1.0/sqrt(15.0*M_PI)*sqrt(y)*exp(-y);

  return d2;
}

double sigmap2_1d3_1p3(double y) {
  // Verified working
  double s2 = 2.0/sqrt(15.0*M_PI)*exp(-y)*sqrt(y);

  return s2;
}

double sigmapp2_1d3_1p3(double y) {
  // Verified working
  double s2 = 2.0/15.0*sqrt(10.0/M_PI)*exp(-y)*(sqrt(y) - pow(y, 3.0/2.0));

  return s2;
}

double omegap2_1d3_1p3(double y) {
  // Verified working
  double o2 = -1.0/6.0*sqrt(10.0/M_PI)*exp(-y)*sqrt(y);

  return o2;
}

double delta2_1d5_1p1(double y) {
  // Verified working
  double d2 = -1.0/15.0*sqrt(10.0/M_PI)*exp(-y)*sqrt(y);

  return d2;
}

double sigmap2_1d5_1p1(double y) {
  // Verified working
  double s2 = 2.0/5.0*sqrt(10.0/M_PI)*exp(-y)*(-sqrt(y) + 1.0/3.0*pow(y, 3.0/2.0));

  return s2;
}

double sigmapp2_1d5_1p1(double y) {
  // Verified working
  double s2 = 4.0/sqrt(15.0*M_PI)*exp(-y)*(-sqrt(y) + 1.0/2.0*pow(y, 3.0/2.0));

  return s2;
}

double omegap2_1d5_1p1(double y) {
  // Verified working 
  double o2 = -1.0/10.0*sqrt(15.0/M_PI)*exp(-y)*sqrt(y);

  return o2;
}

double delta2_1d5_1p3(double y) {
  // Verified working
  double d2 = -1.0/15.0*sqrt(35.0/M_PI)*exp(-y)*sqrt(y);

  return d2;
}

double sigmap2_1d5_1p3(double y) {
  // Verified working
  double s2 = 1.0/5.0*sqrt(35.0/M_PI)*exp(-y)*(sqrt(y) - 10.0/21.0*pow(y, 3.0/2.0));

  return s2;
}

double sigmapp2_1d5_1p3(double y) {
  // Verified working 
  double s2 = 1.0/15.0*sqrt(210.0/M_PI)*exp(-y)*(sqrt(y) - 2.0/7.0*pow(y, 3.0/2.0));

  return s2;
}

double sigmap4_1d5_1p3(double y) {
  // Verified working 
  double s4 = -2.0/sqrt(7.0*M_PI)*exp(-y)*pow(y, 3.0/2.0);

  return s4;
}

double sigmapp4_1d5_1p3(double y) {
  // Verified working 
  double s4 = -4.0/sqrt(35.0*M_PI)*exp(-y)*pow(y, 3.0/2.0);

  return s4;
}

double m0_2s1_2s1(double y) {
 // Verified working 
 double m0 = 1.0/sqrt(2.0*M_PI)*exp(-y)*(1.0 - 4.0/3.0*y + 2.0/3.0*y*y);

 return m0;
}

double m2_1d3_2s1(double y) {
  // Verified working
  double m2 = 4.0/15.0*sqrt(10.0/M_PI)*exp(-y)*(y - 0.5*y*y);

  return m2;
}

double deltap2_1d3_2s1(double y) {
  // Verified working
  double d2 = 1.0/sqrt(15.0*M_PI)*exp(-y)*y;

  return d2;
}

double sigma2_1d3_2s1(double y) {
  // Verified working
  double s2 = 4.0/sqrt(15.0*M_PI)*exp(-y)*(-y + 1.0/2.0*y*y);

  return s2;
}

double m0_1d3_1d3(double y) {
  // Verified working
  double m0 = 1.0/sqrt(M_PI)*exp(-y)*(1.0 - 4.0/3.0*y + 4.0/15.0*y*y);

  return m0;
}

double m2_1d3_1d3(double y) {
  // Verified working
  double m2 = 14.0/(15.0*sqrt(M_PI))*exp(-y)*(-y + 2.0/7.0*y*y);

  return m2;
}

double m2_1d5_2s1(double y) {
  // Verified working
  double m2 = 4.0/sqrt(15.0*M_PI)*exp(-y)*(-y + 1.0/2.0*y*y);

  return m2;
}

double deltap2_1d5_2s1(double y) {
  // Verified working
  double d2 = -1.0/sqrt(10.0*M_PI)*exp(-y)*y;

  return d2;
}

double sigma2_1d5_2s1(double y) {
  // Verified working
  double s2 = 4.0/15.0*sqrt(10.0/M_PI)*exp(-y)*(-y + 1.0/2.0*y*y);

  return s2;
}

double m2_1d5_1d3(double y) {
  // Verified working
  double m2 = 2.0/15.0*sqrt(21.0/M_PI)*exp(-y)*(-y + 2.0/7.0*y*y);

  return m2;
}

double sigma2_1d5_1d3(double y) {
  // Verified working
  double s2 = 1.0/3.0*sqrt(14.0/M_PI)*exp(-y)*(-y + 2.0/7.0*y*y);

  return s2;
}

double m4_1d5_1d3(double y) {
  // Verified working
  double m4 = 4.0/35.0*sqrt(14.0/M_PI)*exp(-y)*y*y;

  return m4;
}

double sigma4_1d5_1d3(double y) {
  // Verified working 
  double s4 = 2.0/35.0*sqrt(70.0/M_PI)*exp(-y)*y*y;

  return s4;
}

double m0_1d5_1d5(double y) {
  // Verified working
  double m0 = 1.0/2.0*sqrt(6.0/M_PI)*exp(-y)*(1.0 - 4.0/3.0*y + 4.0/15.0*y*y);

  return m0;
}

double m2_1d5_1d5(double y) {
  // Verified working
  double m2 = 4.0/15.0*sqrt(21.0/M_PI)*exp(-y)*(-y + 2.0/7.0*y*y);

  return m2;
}

double m4_1d5_1d5(double y) {
  // Verified working
  double m4 = 4.0/35.0*sqrt(7.0/M_PI)*exp(-y)*y*y;

  return m4;
}

double sigmap1_2s1_2s1(double y) {
  // Verified working  
  double s1 = 1.0/sqrt(M_PI)*exp(-y)*(1.0 - 4.0/3.0*y + 2.0/3.0*y*y);

  return s1;
}

double sigmapp1_2s1_2s1(double y) {
  // Verified working
  double s1 = 1.0/sqrt(2.0*M_PI)*exp(-y)*(1.0 - 4.0/3.0*y + 2.0/3.0*y*y);

  return s1;
}

double sigmap1_1d3_2s1(double y) {
  // Verified working
  double s1 = 4.0/15.0*sqrt(5.0/M_PI)*exp(-y)*(y - 1.0/2.0*y*y);

  return s1;
}

double sigmapp1_1d3_2s1(double y) {
  // Verified working
  double s1 = 4.0/15.0*sqrt(10.0/M_PI)*exp(-y)*(-y + 1.0/2.0*y*y);

  return s1;
}

double omegap1_1d3_2s1(double y) {
  // Verified working
  double o1 = 1.0/sqrt(10.0*M_PI)*exp(-y)*y;

  return o1;
}

double delta1_1d3_1d3(double y) {
  // Verified working
  double d1 = 3.0/sqrt(10.0*M_PI)*exp(-y)*(-1.0 + 2.0/5.0*y);

  return d1;
}

double sigmap1_1d3_1d3(double y) {
  // Verified working
  double s1 = 1.0/5.0*sqrt(10.0/M_PI)*exp(-y)*(-1.0 + 34.0/15.0*y - 8.0/15.0*y*y);

  return s1;
}

double sigmapp1_1d3_1d3(double y) {
  // Verified working
  double s1 = 1.0/sqrt(5.0*M_PI)*exp(-y)*(-1.0 - 8.0/15.0*y + 4.0/15.0*y*y);

  return s1;
}

double delta3_1d3_1d3(double y) {
  // Verified working
  double d3 = 4.0/(5.0*sqrt(15*M_PI))*exp(-y)*y;

  return d3;
}

double sigmap3_1d3_1d3(double y) {
  // Verified working
  double s3 = 4.0/(5.0*sqrt(15.0*M_PI))*exp(-y)*(y - 2.0*y*y);

  return s3;
}

double sigmapp3_1d3_1d3(double y) {
  // Verified working
  double s3 = 2.0/(5.0*sqrt(5.0*M_PI))*exp(-y)*(y + 2.0*y*y);

  return s3;
}

double sigmap3_1d5_2s1(double y) {
  // Verified working
  double s3 = 8.0/15.0*sqrt(5.0/M_PI)*exp(-y)*(-y + 1.0/2.0*y*y);

  return s3;
}

double sigmapp3_1d5_2s1(double y) {
  // Verified working
  double s3 = 4.0/sqrt(15.0*M_PI)*exp(-y)*(-y + 1.0/2.0*y*y);

  return s3;
}

double omegap3_1d5_2s1(double y) {
  // Verified working
  double o3 = -1.0/sqrt(15.0*M_PI)*exp(-y)*y;

  return o3;
}

double delta1_1d5_1d3(double y) {
  // Verified working
  double d1 = 1.0/sqrt(10.0*M_PI)*exp(-y)*(-1.0 + 2.0/5.0*y);

  return d1;
}

double sigmap1_1d5_1d3(double y) {
  // Verified working
  double s1 = 2.0/5.0*sqrt(10.0/M_PI)*exp(-y)*(-1.0 + 11.0/10.0*y - 1.0/5.0*y*y);

  return s1;
}

double sigmapp1_1d5_1d3(double y) {
  // Verified working
  double s1 = 2.0/sqrt(5.0*M_PI)*exp(-y)*(-1.0 + 9.0/5.0*y - 2.0/5.0*y*y);

  return s1;
}

double omegap1_1d5_1d3(double y) {
  // Verified working
  double o1 = 1.0/2.0*sqrt(5.0/M_PI)*exp(-y)*(-1.0 + 2.0/5.0*y);

  return o1;
}

double delta3_1d5_1d3(double y) {
  // Verified working
  double d3 = 2.0/25.0*sqrt(10.0/M_PI)*exp(-y)*y;

  return d3;
}

double sigmap3_1d5_1d3(double y) {
  // Verified working
  double s3 = 16.0/75.0*sqrt(10.0/M_PI)*exp(-y)*(y - 1.0/8.0*y*y);

  return s3;
}

double sigmapp3_1d5_1d3(double y) {
  // Verified working
  double s3 = 8.0/75.0*sqrt(30.0/M_PI)*exp(-y)*(y - 1.0/2.0*y*y);

  return s3;
}

double omegap3_1d5_1d3(double y) {
  // Verified working
  double o3 = 1.0/15.0*sqrt(30.0/M_PI)*exp(-y)*y;

  return o3;
}

double delta1_1d5_1d5(double y) {
  // Verified working
  double d1 = 1.0/5.0*sqrt(35.0/M_PI)*exp(-y)*(-1.0 + 2.0/5.0*y);

  return d1;
}

double sigmap1_1d5_1d5(double y) {
  // Verified working
  double s1 = 1.0/5.0*sqrt(35.0/M_PI)*exp(-y)*(1.0 - 8.0/5.0*y + 12.0/35.0*y*y);

  return s1;
}

double sigmapp1_1d5_1d5(double y) {
  // Verified working
  double s1 = 1.0/10.0*sqrt(70.0/M_PI)*exp(-y)*(1.0 - 4.0/5.0*y + 4.0/35.0*y*y);

  return s1;
}

double delta3_1d5_1d5(double y) {
  // Verified working
  double d3 = 2.0/25.0*sqrt(15.0/M_PI)*exp(-y)*y;

  return d3;
}

double sigmap3_1d5_1d5(double y) {
  // Verified working
  double s3 = 8.0/25.0*sqrt(15.0/M_PI)*exp(-y)*(-y + 1.0/3.0*y*y);

  return s3;
}

double sigmapp3_1d5_1d5(double y) {
  // Verified working 
  double s3 = 12.0/25.0*sqrt(5.0/M_PI)*exp(-y)*(-y + 2.0/9.0*y*y);

  return s3;
}

double sigmap5_1d5_1d5(double y) {
  // Verified working
  double s5 = 4.0/105.0*sqrt(210.0/M_PI)*exp(-y)*y*y;

  return s5;
}

double sigmapp5_1d5_1d5(double y) {
  // Verified working
  double s5 = 4.0/(3.0*sqrt(7.0*M_PI))*exp(-y)*y*y;

  return s5;
}

int get_shell_index(int in, int ij) {
  int ind = 0;
  if (in == 0) {
    ind = 0;
  } else if (in == 1) {
    if (ij == 1) {
      ind = 1;
    } else if (ij == 3) {
      ind = 2;
    }
  } else if (in == 2) {
    if (ij == 1) {
      ind = 3;
    } else if (ij == 3) {
      ind = 4;
    } else if (ij == 5) {
      ind = 5;
    }
  }
  return ind;
}

double b_osc(int a_nuc) {
  double b = sqrt(0.9*pow(a_nuc, 1.0/3.0) + 0.7);
  
  return b;
}
