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
  mu2e();
  return 0;
}
