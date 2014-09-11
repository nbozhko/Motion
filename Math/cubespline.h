#ifndef CUBESPLINE_H
#define CUBESPLINE_H

#include <vector>

//This class is CubeSpline class. It creates Spline
//from points with constant step
class CubeSpline
{
    //table of numbers is 2xN
    std::vector<double> X, Y;
    std::vector<double> q;
    double h;
public:
    CubeSpline() {h = 0;}
    CubeSpline(const std::vector<double> &x, const std::vector<double> &y);
    double GetValueOfFunction(const double &x) const;

    void operator=(const CubeSpline &spline);
    ~CubeSpline();
};

#endif // CUBESPLINE_H
