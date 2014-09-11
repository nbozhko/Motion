#ifndef ODESOLVEADAMS_H
#define ODESOLVEADAMS_H

#include "cubespline.h"
#include "matrix3x3.h"
#include "vector3d.h"
#include <cmath>

#define TO_RADS(grad, mins, secs) (((grad + mins / 60.0 + secs / 3600.00) * 3.14159265) / 180.0)

struct SolveOfEquation{
    double *X;
    double **Y;
    unsigned size;
    int numOfEquations;
    double step;
    double **dOmega;
    short numOmega;
};

class ODESolveAdams
{
private:
    double *X;
    double **Y;
    double **dOmega;
    double step;
    CubeSpline *splineForce;
    vector3d *vecForce1;
    vector3d *vecForce;
    matrix3x3 tensorOfInertia;
    matrix3x3 matRotate;
    vector3d centerofmass;
    double mass;
    struct vec4 {
        double lambda0, lambda1, lambda2, lambda3;
    } lambda;
    int NumOfEquations;
    unsigned SizeofArrayPoints;
    double force_projX, force_projY, force_projZ;

    void DeleteData();
public:
    ODESolveAdams(const double *Y0, const double &x0, const double &xk, const double &st, int num);
    void SetForce(const TableOfPoints &table);
    void SetTensorOfInertia(const double *tensor);
    void SetTensorOfInertia(const matrix3x3 &tensor);
    void SetMass(const double &m);
    void SetNewInitialCondition(const double *Y0, const double &x0, const double &xk, const double &st);
    void SolveODE();
    void GetSolveODE(SolveOfEquation *solve) const;
    void ClearSolveStruct(SolveOfEquation *solve) const;
    ~ODESolveAdams();
};

#endif // ODESOLVEADAMS_H
