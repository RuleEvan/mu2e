#include "wfn.h"

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
  double *y = (double*) malloc(sizeof(double)*4);
  y[0] = pow(0.1*b_osc(27)/(HBARC*2.0), 2.0);
  y[1] = pow(1.0*b_osc(27)/(HBARC*2.0), 2.0);
  y[2] = pow(10.0*b_osc(27)/(HBARC*2.0), 2.0);
  y[3] = pow(100.0*b_osc(27)/(HBARC*2.0), 2.0);

 for (int i = 0; i < 4; i++) {
   printf("%g\n", m1_2s1_1p1(y[i]));
 }
 // generate_potential_files();
 
  return 0;
}
