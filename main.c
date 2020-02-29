#include "matrix_element.h"

int main(int argc, char *argv[]) {
/*  int i_eig = 0;
  wfnData *wd = read_binary_wfn_data("al27_usdb.wfn", "al27_usdb.bas", i_eig);
  printf("J: %g\n", wd->j_nuc[i_eig]);
  double b_osc = 1.8439; // [fm]
  double total = 0.0;
  for (int i = 0; i < 1000; i++) {
    double q = i*1.0;
    total += q*q*radial_osc_wfn_q(1, 2, q, b_osc)*radial_osc_wfn_q(1, 2, q, b_osc);
    printf("%g, %g\n", q/M_NEUTRON, pow(radial_osc_wfn_q(0, 2, q, b_osc), 2.0));
  }
  printf("total: %g\n", total);
  */
//  mu2e();
/*  double *y = (double*) malloc(sizeof(double)*4);
  y[0] = pow(0.1*b_osc(27)/(HBARC*2.0), 2.0);
  y[1] = pow(1.0*b_osc(27)/(HBARC*2.0), 2.0);
  y[2] = pow(10.0*b_osc(27)/(HBARC*2.0), 2.0);
  y[3] = pow(100.0*b_osc(27)/(HBARC*2.0), 2.0);

 for (int i = 0; i < 4; i++) {
   printf("%g %g\n", y[i], sigmapp5_1d5_1d5(y[i]));
 }
 */
  compute_isotope_multipoles(105.0);

//  printf("%g\n", compute_rel_potential(1.0, 0.0, 1.0, 0.0, 0, 105.0, 2));
//  printf("%g\n", compute_radial_matrix_element_J_dot_J(0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 105.0));
/*  double q_min = 0.0;
  double q_max = 500.0;
  int num_q_steps = 100.0;
  double delta_q = (q_max - q_min)/((double) num_q_steps);

  for (int i = 0; i < num_q_steps; i++) {
    double q = q_min + i*delta_q;
  }
*/
  //printf("q = 0: %g\n", -compute_matrix_element_sigma_0(2, 3, 2, 3, 0, 2, 3, 2, 3, 0, 2)/(4.0*M_PI));
  printf("Op 5: %g\n", compute_total_matrix_element_sigma_0_finite_q_op5("isotope_data/al27/density/al27-al27_core_J0_T0_0_0.dens", 105.0, 1, 0)/sqrt(2.0));

  printf("New tensor op: %g\n", compute_total_matrix_element_y2_finite_q("isotope_data/al27/density/al27-al27_core_J0_T0_0_0.dens", 105.0, 2, 0)/sqrt(2.0));
  
  printf("Tensor op\n");
  printf("q != 0: %g\n", compute_total_matrix_element_y2_finite_q_alt("isotope_data/al27/density/al27-al27_core_J0_T0_0_0.dens", 105.0, 2, 0)/sqrt(2.0));
  printf("q = 0: %g\n", compute_total_matrix_element_TT("isotope_data/al27/density/al27-al27_core_J0_T0_0_0.dens", 3)/(2.0*sqrt(M_PI)*sqrt(2.0)));

  printf("Sigma op\n");
  printf("q != 0: %g\n", compute_total_matrix_element_sigma_0_finite_q("isotope_data/al27/density/al27-al27_core_J0_T0_0_0.dens", 105.0, 0)*(-1.0/(12.0*sqrt(M_PI)*sqrt(2.0))));
  printf("q = 0: %g\n", compute_total_matrix_element_sigma_0("isotope_data/al27/density/al27-al27_core_J0_T0_0_0.dens")/(24.0*M_PI*sqrt(2.0)*2.0*sqrt(M_PI)));


 // generate_potential_files();
 
  return 0;
}
