#include "multipole.h"

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
  double d1 = -1.0/(3.0*sqrt(4.0*M_PI))*exp(-y)*(1.0/sqrt(y) + sqrt(y));

  return d1;
}

double sigma1_2s1_1p1(double y) {
  double s1 = 4.0/(3.0*sqrt(M_PI))*exp(-y)*(sqrt(y) - pow(y, 3.0/2.0));

  return s1;
}

double m1_2s1_1p3(double y) {
  double m1 = 2.0/(3.0*sqrt(M_PI))*exp(-y)*(sqrt(y) - pow(y, 3.0/2.0));

  return m1;
}

double deltap1_2s1_1p3(double y) {
  double d1 = -1.0/(3.0*sqrt(2.0*M_PI))*exp(-y)*(1.0/sqrt(y) + sqrt(y));

  return d1;
}

double sigma1_2s1_1p3(double y) {
  double s1 = 1.0/3.0*sqrt(2.0/M_PI)*exp(-y)*(-sqrt(y) + pow(y, 3.0/2.0));

  return s1;
}

double m1_1d3_1p1(double y) {
  double m1 = 1.0/3.0*sqrt(10.0/M_PI)*exp(-y)*(sqrt(y) - 2.0/5.0*pow(y, 3.0/2.0));

  return m1;
}

double deltap1_1d3_1p1(double y) {
  double d1 = 2.0/3.0*sqrt(5.0/M_PI)*exp(-y)*(-1.0/sqrt(y) + 4.0/5.0*sqrt(y));

  return d1;
}

double sigma1_1d3_1p1(double y) {
  double s1 = 1.0/3.0*sqrt(5.0/M_PI)*exp(-y)*(-sqrt(y) + 2.0/5.0*pow(y, 3.0/2.0));

  return s1;
}

double m1_1d3_1p3(double y) {
  double m1 = 1.0/3.0*sqrt(2.0/M_PI)*exp(-y)*(-sqrt(y) + 2.0/5.0*pow(y, 3.0/2.0));

  return m1;
}

double deltap1_1d3_1p3(double y) {
  double d1 = 1.0/(6.0*sqrt(M_PI))*exp(-y)*(1.0/sqrt(y) - 4.0/5.0*sqrt(y));

  return d1;
}

double sigma1_1d3_1p3(double y) {
  double s1 = 4.0/(3.0*sqrt(M_PI))*exp(-y)*(sqrt(y) - 2.0/5.0*pow(y, 3.0/2.0));
  
  return s1;
}

double m3_1d3_1p3(double y) {
  double m3 = 2.0/5.0*sqrt(2.0/M_PI)*exp(-y)*pow(y, 3.0/2.0);

  return m3;
}

double deltap3_1d3_1p3(double y) {
  double d3 = -1.0/15.0*sqrt(6.0/M_PI)*exp(-y)*sqrt(y);

  return d3;
}

double sigma3_1d3_1p3(double y) {
  double s3 = -4.0/15.0*sqrt(6.0/M_PI)*exp(-y)*pow(y, 3.0/2.0);

  return s3;
}

double m3_1d5_1p1(double y) {
  double m3 = -2.0/sqrt(15.0*M_PI)*exp(-y)*pow(y, 3.0/2.0);

  return m3;
}

double deltap3_1d5_1p1(double y) {
  double d3 = 1.0/15.0*sqrt(5.0/M_PI)*exp(-y)*sqrt(y);

  return d3;
}

double sigma3_1d5_1p1(double y) {
  double s3 = -4.0/15.0*sqrt(15.0/M_PI)*exp(-y)*pow(y, 3.0/2.0);

  return s3;
}

double m1_1d5_1p3(double y) {
  double m1 = sqrt(2.0/M_PI)*exp(-y)*(sqrt(y) - 2.0/5.0*pow(y, 3.0/2.0));

  return m1;
}

double deltap1_1d5_1p3(double y) {
  double d1 = 1.0/sqrt(4.0*M_PI)*exp(-y)*(-1.0/sqrt(y) + 4.0/5.0*sqrt(y));

  return d1;
}

double sigma1_1d5_1p3(double y) {
  double s1 = 1.0/sqrt(M_PI)*exp(-y)*(sqrt(y) - 2.0/5.0*pow(y, 3.0/2.0));

  return s1;
}

double m3_1d5_1p3(double y) {
  double m3 = -4.0/15.0*sqrt(3.0/M_PI)*exp(-y)*pow(y, 3.0/2.0);

  return m3;
}

double deltap3_1d5_1p3(double y) {
  double d3 = 2.0/(15.0*sqrt(M_PI))*exp(-y)*sqrt(y);

  return d3;
}

double sigma3_1d5_1p3(double y) {
  double s3 = -2.0/(15.0*sqrt(M_PI))*exp(-y)*pow(y, 3.0/2.0);

  return s3;
}

double sigmapp0_2s1_1p1(double y) {
  double s0 = 1.0/3.0*sqrt(2.0/M_PI)*exp(-y)*(-sqrt(y) + pow(y, 3.0/2.0));

  return s0;
}

double omegap0_2s1_1p1(double y) {
  double o0 = -1.0/4.0*sqrt(2.0/M_PI)*exp(-y)*(1.0/sqrt(y) + 1.0/3.0*sqrt(y));

  return o0;
}

double sigmap2_2s1_1p3(double y) {
  double s2 = 1.0/3.0*sqrt(6.0/M_PI)*exp(-y)*(sqrt(y) - pow(y, 3.0/2.0));

  return s2;
}

double sigmapp2_2s1_1p3(double y) {
  double s2 = 2.0/(3.0*sqrt(M_PI))*exp(-y)*(sqrt(y) - pow(y, 3.0/2.0));

  return s2;
}

double omegap2_2s1_1p3(double y) {
  double o2 = -1.0/(3.0*sqrt(M_PI))*exp(-y)*sqrt(y);

  return o2;
}

double delta2_1d3_1p1(double y) {
  double d2 = -1.0/sqrt(15.0*M_PI)*sqrt(y)*exp(-y);

  return d2;
}

double sigmap2_1d3_1p1(double y) {
  double s2 = 1.0/sqrt(15.0*M_PI)*exp(-y)*(-sqrt(y) + 2.0*pow(y, 3.0/2.0));

  return s2;
}

double sigmapp2_1d3_1p1(double y) {
  double s2 = -1.0/15.0*sqrt(10.0/M_PI)*exp(-y)*(sqrt(y) + 2.0*pow(y, 3.0/2.0));

  return s2;
}

double omegap2_1d3_1p1(double y) {
  double o2 = -1.0/15.0*sqrt(10.0/M_PI)*exp(-y)*sqrt(y);

  return o2;
}

double sigmapp0_1d3_1p3(double y) {
  double s0 = 1.0/3.0*sqrt(10.0/M_PI)*exp(-y)*(sqrt(y) - 2.0/5.0*pow(y, 3.0/2.0));

  return s0;
}

double omegap0_1d3_1p3(double y) {
  double o0 = 1.0/4.0*sqrt(10.0/M_PI)*exp(-y)*(1.0/sqrt(y) - 2.0/3.0*sqrt(y));

  return o0;
}

double delta2_1d3_1p3(double y) {
  double d2 = 1.0/sqrt(15.0*M_PI)*sqrt(y)*exp(-y);

  return d2;
}

double sigmap2_1d3_1p3(double y) {
  double s2 = 2.0/sqrt(15.0*M_PI)*exp(-y)*sqrt(y);

  return s2;
}

double sigmapp2_1d3_1p3(double y) {
  double s2 = 2.0/15.0*sqrt(10.0/M_PI)*exp(-y)*(sqrt(y) - pow(y, 3.0/2.0));

  return s2;
}

double omegap2_1d3_1p3(double y) {
  double o2 = -2.0/3.0*sqrt(10.0/M_PI)*exp(-y)*sqrt(y);

  return o2;
}

double delta2_1d5_1p1(double y) {
  double d2 = -1.0/15.0*sqrt(10.0/M_PI)*exp(-y)*sqrt(y);

  return d2;
}

double sigmap2_1d5_1p1(double y) {
  double s2 = 2.0/5.0*sqrt(10.0/M_PI)*exp(-y)*(-sqrt(y) + 1.0/3.0*pow(y, 3.0/2.0));

  return s2;
}

double sigmapp2_1d5_1p1(double y) {
  double s2 = 2.0/sqrt(15.0*M_PI)*exp(-y)*(-sqrt(y) + 1.0/2.0*pow(y, 3.0/2.0));

  return s2;
}

double omegap2_1d5_1p1(double y) {
  double o2 = -1.0/10.0*sqrt(15.0/M_PI)*exp(-y)*sqrt(y);

  return o2;
}

double delta2_1d5_1p3(double y) {
  double d2 = -1.0/15.0*sqrt(35.0/M_PI)*exp(-y)*sqrt(y);

  return d2;
}

double sigmap2_1d5_1p3(double y) {
  double s2 = 1.0/5.0*sqrt(35.0/M_PI)*exp(-y)*(sqrt(y) - 10.0/21.0*pow(y, 3.0/2.0));

  return s2;
}

double sigmapp2_1d5_1p3(double y) {
  double s2 = 1.0/15.0*sqrt(210.0/M_PI)*exp(-y)*(sqrt(y) - 2.0/7.0*pow(y, 3.0/2.0));

  return s2;
}

double sigmap4_1d5_1p3(double y) {
  double s4 = -2.0/sqrt(7.0*M_PI)*exp(-y)*pow(y, 3.0/2.0);

  return s4;
}

double sigmapp4_1d5_1p3(double y) {
  double s4 = -4.0/sqrt(35.0*M_PI)*exp(-y)*pow(y, 3.0/2.0);

  return s4;
}

double m0_2s1_2s1(double y) {
 double m0 = 1.0/sqrt(2.0*M_PI)*exp(-y)*(1.0 - 4.0/3.0*y + 2.0/3.0*y*y);

 return m0;
}

double m2_1d3_2s1(double y) {
  double m2 = 4.0/15.0*sqrt(10.0/M_PI)*exp(-y)*(y - 0.5*y*y);

  return m2;
}

double deltap2_1d3_2s1(double y) {
  double d2 = 1.0/sqrt(15.0*M_PI)*exp(-y)*y;

  return d2;
}

double sigma2_1d3_2s1(double y) {
  double s2 = 4.0/sqrt(15.0*M_PI)*exp(-y)*(-y + 1.0/2.0*y*y);

  return s2;
}

double m0_1d3_1d3(double y) {
  double m0 = 1.0/sqrt(M_PI)*exp(-y)*(1.0 - 4.0/3.0*y + 4.0/15.0*y*y);

  return m0;
}

double m2_1d3_1d3(double y) {
  double m2 = 14.0/(35.0*sqrt(M_PI))*exp(-y)*(-y + 2.0/7.0*y*y);

  return m2;
}

double m2_1d5_2s1(double y) {
  double m2 = 4.0/sqrt(15.0*M_PI)*exp(-y)*(-y + 1.0/2.0*y*y);

  return m2;
}

double deltap2_1d5_2s1(double y) {
  double d2 = -1.0/sqrt(10.0*M_PI)*exp(-y)*y;

  return d2;
}

double sigma2_1d5_2s1(double y) {
  double s2 = 4.0/15.0*sqrt(10.0/M_PI)*exp(-y)*(-y + 1.0/2.0*y*y);

  return s2;
}

double m2_1d5_1d3(double y) {
  double m2 = 2.0/15.0*sqrt(21.0/M_PI)*exp(-y)*(-y + 2.0/7.0*y*y);

  return m2;
}

double sigma2_1d5_1d3(double y) {
  double s2 = 1.0/3.0*sqrt(14.0/M_PI)*exp(-y)*(-y + 2.0/7.0*y*y);

  return s2;
}

double m4_1d5_1d3(double y) {
  double m4 = 4.0/35.0*sqrt(14.0/M_PI)*exp(-y)*y*y;

  return m4;
}

double sigma4_1d5_1d3(double y) {
  double s4 = 2.0/35.0*sqrt(70.0/M_PI)*exp(-y)*y*y;

  return s4;
}

double m0_1d5_1d5(double y) {
  double m0 = 1.0/2.0*sqrt(6.0/M_PI)*exp(-y)*(1.0 - 4.0/3.0*y + 4.0/15.0*y*y);

  return m0;
}

double m2_1d5_1d5(double y) {
  double m2 = 2.0/35.0*sqrt(21.0/M_PI)*exp(-y)*(-y + 2.0/7.0*y*y);

  return m2;
}

double m4_1d5_1d5(double y) {
  double m4 = 2.0/35.0*sqrt(7.0/M_PI)*exp(-y)*y*y;

  return m4;
}

double sigmap1_2s1_2s1(double y) {
  double s1 = 1.0/sqrt(M_PI)*exp(-y)*(1.0 - 4.0/3.0*y + 2.0/3.0*y*y);

  return s1;
}

double sigmapp1_2s1_2s1(double y) {
  double s1 = 1.0/sqrt(2.0*M_PI)*exp(-y)*(1.0 - 4.0/3.0*y + 2.0/3.0*y*y);

  return s1;
}

double sigmap1_1d3_2s1(double y) {
  double s1 = 4.0/15.0*sqrt(5.0/M_PI)*exp(-y)*(y - 1.0/2.0*y*y);

  return s1;
}

double sigmapp1_1d3_2s1(double y) {
  double s1 = 4.0/15.0*sqrt(10.0/M_PI)*exp(-y)*(-y + 1.0/2.0*y*y);

  return s1;
}

double omegap1_1d3_2s1(double y) {
  double o1 = 1.0/sqrt(10.0*M_PI)*exp(-y)*y;

  return o1;
}

double delta1_1d3_1d3(double y) {
  double d1 = 3.0/sqrt(10.0*M_PI)*exp(-y)*(-1.0 + 2.0/5.0*y);

  return d1;
}

double sigmap1_1d3_1d3(double y) {
  double s1 = 1.0/5.0*sqrt(10.0/M_PI)*exp(-y)*(-1.0 + 34.0/15.0*y - 8.0/15.0*y*y);

  return s1;
}

double sigmapp1_1d3_1d3(double y) {
  double s1 = 1.0/sqrt(5.0*M_PI)*exp(-y)*(-1.0 - 8.0/15.0*y + 4.0/15.0*y*y);

  return s1;
}

double delta3_1d3_1d3(double y) {
  double d3 = 4.0/(5.0*sqrt(15*M_PI))*exp(-y)*y;

  return d3;
}

double sigmap3_1d3_1d3(double y) {
  double s3 = 4.0/(5.0*sqrt(15.0*M_PI))*exp(-y)*(y - 2.0*y*y);

  return s3;
}

double sigmapp3_1d3_1d3(double y) {
  double s3 = 2.0/(5.0*sqrt(5.0*M_PI))*exp(-y)*(y + 2.0*y*y);

  return s3;
}

double sigmap3_1d5_2s1(double y) {
  double s3 = 8.0/15.0*sqrt(5.0/M_PI)*exp(-y)*(-y + 1.0/2.0*y*y);

  return s3;
}

double sigmapp3_1d5_2s1(double y) {
  double s3 = 4.0/sqrt(15.0*M_PI)*exp(-y)*(-y + 1.0/2.0*y*y);

  return s3;
}

double omegap3_1d5_2s1(double y) {
  double o3 = -1.0/sqrt(15.0*M_PI)*exp(-y)*y;

  return o3;
}

double delta1_1d5_1d3(double y) {
  double d1 = 1.0/sqrt(10.0*M_PI)*exp(-y)*(-1.0 + 2.0/5.0*y);

  return d1;
}

double sigmap1_1d5_1d3(double y) {
  double s1 = 2.0/5.0*sqrt(10.0/M_PI)*exp(-y)*(-1.0 + 11.0/10.0*y - 1.0/5.0*y*y);

  return s1;
}

double sigmapp1_1d5_1d3(double y) {
  double s1 = 2.0/sqrt(5.0*M_PI)*exp(-y)*(-1.0 + 9.0/5.0*y - 2.0/5.0*y*y);

  return s1;
}

double omegap1_1d5_1d3(double y) {
  double o1 = 1.0/2.0*sqrt(5.0/M_PI)*exp(-y)*(-1.0 + 2.0/5.0*y);

  return o1;
}

double delta3_1d5_1d3(double y) {
  double d3 = 2.0/25.0*sqrt(10.0/M_PI)*exp(-y)*y;

  return d3;
}

double sigmap3_1d5_1d3(double y) {
  double s3 = 4.0/15.0*sqrt(10.0/M_PI)*exp(-y)*(y - 1.0/4.0*y*y);

  return s3;
}

double sigmapp3_1d5_1d3(double y) {
  double s3 = 8.0/75.0*sqrt(30.0/M_PI)*exp(-y)*(y - 1.0/2.0*y*y);

  return s3;
}

double omegap3_1d5_1d3(double y) {
  double o3 = 1.0/15.0*sqrt(30.0/M_PI)*exp(-y)*y;

  return o3;
}

double delta1_1d5_1d5(double y) {
  double d1 = 1.0/5.0*sqrt(35.0/M_PI)*exp(-y)*(-1.0 + 2.0/5.0*y);

  return d1;
}

double sigmap1_1d5_1d5(double y) {
  double s1 = 1.0/5.0*sqrt(35.0/M_PI)*exp(-y)*(1.0 - 4.0/5.0*y + 12.0/35.0*y*y);

  return s1;
}

double sigmapp1_1d5_1d5(double y) {
  double s1 = 1.0/10.0*sqrt(70.0/M_PI)*exp(-y)*(1.0 - 4.0/5.0*y + 4.0/35.0*y*y);

  return s1;
}

double delta3_1d5_1d5(double y) {
  double d3 = 2.0/25.0*sqrt(15.0/M_PI)*exp(-y)*y;

  return d3;
}

double sigmap3_1d5_1d5(double y) {
  double s3 = 8.0/25.0*sqrt(15.0/M_PI)*exp(-y)*(-y + 1.0/3.0*y*y);

  return s3;
}

double sigmapp3_1d5_1d5(double y) {
  double s3 = 12.0/25.0*sqrt(5.0/M_PI)*exp(-y)*(-y + 2.0/9.0*y*y);

  return s3;
}

double sigmap5_1d5_1d5(double y) {
  double s5 = 4.0/105.0*sqrt(210.0/M_PI)*exp(-y)*y*y;

  return s5;
}

double sigmapp5_1d5_1d5(double y) {
  double s5 = 2.0/(3.0*sqrt(7.0*M_PI))*exp(-y)*y*y;

  return s5;
}

