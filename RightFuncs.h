#ifndef RIGHTFUNCS_H
#define RIGHTFUNCS_H

#include "Math/vector3d.h"
#include "Math/matrix3x3.h"
#include "Math/cubespline.h"
#include <boost/math/special_functions.hpp>
#include <vector>

#define TO_RADS(grad, mins, secs) (((grad + mins / 60.0 + secs / 3600.00) * 3.14159265) / 180.0)

extern CubeSpline *spForce;

extern double alpha, gamma, beta;

extern vector3d    Force,
            radius;
extern matrix3x3   Rotate,
            Inertia;
extern double mass;

vector3d Multiply(const matrix3x3 &mat, const vector3d &vec);

void updateRotateMatrix(matrix3x3 &rot, double l0, double l1, double l2, double l3);

double dx(const std::vector<double> &args);
double dy(const std::vector<double> &args);
double dz(const std::vector<double> &args);

double ddx(const std::vector<double> &args);
double ddy(const std::vector<double> &args);
double ddz(const std::vector<double> &args);

double dwx(const std::vector<double> &args);
double dwy(const std::vector<double> &args);
double dwz(const std::vector<double> &args);

double dlambda0(const std::vector<double> &args);
double dlambda1(const std::vector<double> &args);
double dlambda2(const std::vector<double> &args);
double dlambda3(const std::vector<double> &args);

double KineticIntegral(double wx, double wy, double wz);
double EnergyIntegral(double wx, double wy, double wz);

double AnalyticOmegaX(double wx0, double wy0, double wz0, double t0, double t);
double AnalyticOmegaY(double wx0, double wy0, double wz0, double t0, double t);
double AnalyticOmegaZ(double wx0, double wy0, double wz0, double t0, double t);

CubeSpline Pitch(const std::vector<std::vector<double> > &args, double t0, double h);
CubeSpline Yaw(const std::vector<std::vector<double> > &args, double t0, double h);

#endif // RIGHTFUNCS_H
