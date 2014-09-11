#ifndef RKF45_H
#define RKF45_H

#include <vector>

namespace NMath {

  typedef double (*pfunc)(const std::vector<double> &vars);

  void rkf45(std::vector< std::vector<double> > &solve, const std::vector<double> &vars, const double &h,
             const double &x_beg, const double &x_end, const std::vector<pfunc> &RightFuncs);
  double Func(const double *wcoeff, const double *ki, unsigned short num);
  double ki(const std::vector<double> &vars, unsigned short i, const double &h, const std::vector<pfunc> &RightFuncs, unsigned j);

}

#endif // RKF45_H
