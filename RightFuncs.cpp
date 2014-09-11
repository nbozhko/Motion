#include "RightFuncs.h"

using namespace std;

CubeSpline *spForce;

double alpha, gamma, beta;

vector3d    Force,
            radius;
matrix3x3   Rotate,
            Inertia;
double mass;

vector3d Multiply(const matrix3x3 &mat, const vector3d &vec)
{
    vector3d temp;

    temp.SetX(mat[0] * vec.GetX() + mat[1] * vec.GetY() + mat[2] * vec.GetZ());
    temp.SetY(mat[3] * vec.GetX() + mat[4] * vec.GetY() + mat[5] * vec.GetZ());
    temp.SetZ(mat[6] * vec.GetX() + mat[7] * vec.GetY() + mat[8] * vec.GetZ());

    return temp;
}

void updateRotateMatrix(matrix3x3 &rot, double l0, double l1, double l2, double l3)
{
    Rotate[0] = 2.0 * (l0*l0 + l1*l1) - 1.0; Rotate[1] = 2.0 * (l1*l2 - l0*l3); Rotate[2] = 2.0 * (l1*l3 + l0*l2);
    Rotate[3] = 2.0 * (l1*l2 + l0*l3); Rotate[4] = 2.0 * (l0*l0 + l2*l2) - 1.0; Rotate[5] = 2.0 * (l2*l3 - l0*l1);
    Rotate[6] = 2.0 * (l1*l3 - l0*l2); Rotate[7] = 2.0 * (l2*l3 + l0*l1); Rotate[8] = 2.0 * (l0*l0 + l3*l3) - 1.0;
}

double dx(const vector<double> &args)
{
    return args[3];
}

double dy(const vector<double> &args)
{
    return args[4];
}

double dz(const vector<double> &args)
{
    return args[5];
}

double ddx(const vector<double> &args)
{
    double fx = -spForce->GetValueOfFunction(args[13]) * cos(gamma) * cos(alpha);
    double fy = -spForce->GetValueOfFunction(args[13]) * (sin(alpha) * cos(gamma) * cos(beta)
                                                         + sin(gamma) * sin(beta));
    double fz = -spForce->GetValueOfFunction(args[13]) * (sin(alpha) * cos(gamma) * sin(beta)
                                                         - sin(gamma) * cos(beta));
    updateRotateMatrix(Rotate, args[9], args[10], args[11], args[12]);
    Force.SetNewCoords(fx, fy, fz);

    return Multiply(Rotate, Force).GetX() / mass;
}

double ddy(const vector<double> &args)
{
    double fx = -spForce->GetValueOfFunction(args[13]) * cos(gamma) * cos(alpha);
    double fy = -spForce->GetValueOfFunction(args[13]) * (sin(alpha) * cos(gamma) * cos(beta)
                                                         + sin(gamma) * sin(beta));
    double fz = -spForce->GetValueOfFunction(args[13]) * (sin(alpha) * cos(gamma) * sin(beta)
                                                         - sin(gamma) * cos(beta));
    updateRotateMatrix(Rotate, args[9], args[10], args[11], args[12]);
    Force.SetNewCoords(fx, fy, fz);

    return Multiply(Rotate, Force).GetY() / mass;
}

double ddz(const vector<double> &args)
{
    double fx = -spForce->GetValueOfFunction(args[13]) * cos(gamma) * cos(alpha);
    double fy = -spForce->GetValueOfFunction(args[13]) * (sin(alpha) * cos(gamma) * cos(beta)
                                                         + sin(gamma) * sin(beta));
    double fz = -spForce->GetValueOfFunction(args[13]) * (sin(alpha) * cos(gamma) * sin(beta)
                                                         - sin(gamma) * cos(beta));
    updateRotateMatrix(Rotate, args[9], args[10], args[11], args[12]);
    Force.SetNewCoords(fx, fy, fz);

    return Multiply(Rotate, Force).GetZ() / mass;
}

double dwx(const vector<double> &args)
{
    double fx = -spForce->GetValueOfFunction(args[13]) * cos(gamma) * cos(alpha);
    double fy = -spForce->GetValueOfFunction(args[13]) * (sin(alpha) * cos(gamma) * cos(beta)
                                                         + sin(gamma) * sin(beta));
    double fz = -spForce->GetValueOfFunction(args[13]) * (sin(alpha) * cos(gamma) * sin(beta)
                                                         - sin(gamma) * cos(beta));
    Force.SetNewCoords(fx, fy, fz);

    return ( radius.cross(Force).GetX() / Inertia[0] +
            (Inertia[4] - Inertia[8]) * args[7] * args[8] / Inertia[0] );
}

double dwy(const vector<double> &args)
{
    double fx = -spForce->GetValueOfFunction(args[13]) * cos(gamma) * cos(alpha);
    double fy = -spForce->GetValueOfFunction(args[13]) * (sin(alpha) * cos(gamma) * cos(beta)
                                                         + sin(gamma) * sin(beta));
    double fz = -spForce->GetValueOfFunction(args[13]) * (sin(alpha) * cos(gamma) * sin(beta)
                                                         - sin(gamma) * cos(beta));
    Force.SetNewCoords(fx, fy, fz);

    return ( radius.cross(Force).GetY() / Inertia[4] +
            (Inertia[8] - Inertia[0]) * args[6] * args[8] / Inertia[4] );
}

double dwz(const vector<double> &args)
{
    double fx = -spForce->GetValueOfFunction(args[13]) * cos(gamma) * cos(alpha);
    double fy = -spForce->GetValueOfFunction(args[13]) * (sin(alpha) * cos(gamma) * cos(beta)
                                                         + sin(gamma) * sin(beta));
    double fz = -spForce->GetValueOfFunction(args[13]) * (sin(alpha) * cos(gamma) * sin(beta)
                                                         - sin(gamma) * cos(beta));
    Force.SetNewCoords(fx, fy, fz);

    return ( radius.cross(Force).GetZ() / Inertia[8] +
            (Inertia[0] - Inertia[4]) * args[6] * args[7] / Inertia[8] );
}

double dlambda0(const vector<double> &args)
{
    return ( 0.5 * (-args[6] * args[10] - args[7] * args[11] - args[8] * args[12]) );
}

double dlambda1(const vector<double> &args)
{
    return ( 0.5 * (args[6] * args[9] - args[7] * args[12] + args[8] * args[11]) );
}

double dlambda2(const vector<double> &args)
{
    return ( 0.5 * (args[7] * args[9] - args[8] * args[10] + args[6] * args[12]) );
}

double dlambda3(const vector<double> &args)
{
    return ( 0.5 * (args[8] * args[9] + args[7] * args[10] - args[6] * args[11]) );
}

double KineticIntegral(double wx, double wy, double wz)
{
    return Multiply(Inertia, vector3d(wx, wy, wz)).GetLength();
}

double EnergyIntegral(double wx, double wy, double wz)
{
    return vector3d(wx, wy, wz).dot(Multiply(Inertia, vector3d(wx, wy, wz)));
}

inline double sign(double arg)
{
    return (arg < 0) ? (-1.0) : (1.0);
}

double AnalyticOmegaX(double wx0, double wy0, double wz0, double t0, double t)
{
    using namespace boost::math;

    double T = EnergyIntegral(wx0, wy0, wz0) / 2.0;
    double D = KineticIntegral(wx0, wy0, wz0) * KineticIntegral(wx0, wy0, wz0) / (2 * T);

    double a = sqrt(2 * T * (D - Inertia[0]) / (Inertia[4] * (Inertia[4] - Inertia[0])));
    double b = sqrt(2 * T * (Inertia[8] - D) / (Inertia[4] * (Inertia[8] - Inertia[4])));
    double k = a / b;
    double x0 = ellint_1(k, asin(wy0 / a));
    double tau = (t - t0) * sign(wx0 * wz0) * b * sqrt( (Inertia[8] - Inertia[4]) * (Inertia[4] - Inertia[0]) / (Inertia[0] * Inertia[8]) );


    return ( sign(wx0) * sqrt( 2 * T * (Inertia[8] - D) / (Inertia[0] * (Inertia[8] - Inertia[0])) ) * jacobi_dn(k, tau + x0) );
}

double AnalyticOmegaY(double wx0, double wy0, double wz0, double t0, double t)
{
    using namespace boost::math;

    double T = EnergyIntegral(wx0, wy0, wz0) / 2.0;
    double D = KineticIntegral(wx0, wy0, wz0) * KineticIntegral(wx0, wy0, wz0) / (2 * T);

    double a = sqrt(2 * T * (D - Inertia[0]) / (Inertia[4] * (Inertia[4] - Inertia[0])));
    double b = sqrt(2 * T * (Inertia[8] - D) / (Inertia[4] * (Inertia[8] - Inertia[4])));
    double k = a / b;
    double x0 = ellint_1(k, asin(wy0 / a));
    double tau = (t - t0) * sign(wx0 * wz0) * b * sqrt( (Inertia[8] - Inertia[4]) * (Inertia[4] - Inertia[0]) / (Inertia[0] * Inertia[8]) );


    return ( a * jacobi_sn(k, tau + x0) );
}

double AnalyticOmegaZ(double wx0, double wy0, double wz0, double t0, double t)
{
    using namespace boost::math;

    double T = EnergyIntegral(wx0, wy0, wz0) / 2.0;
    double D = KineticIntegral(wx0, wy0, wz0) * KineticIntegral(wx0, wy0, wz0) / (2 * T);

    double a = sqrt(2 * T * (D - Inertia[0]) / (Inertia[4] * (Inertia[4] - Inertia[0])));
    double b = sqrt(2 * T * (Inertia[8] - D) / (Inertia[4] * (Inertia[8] - Inertia[4])));
    double k = a / b;
    double x0 = ellint_1(k, asin(wy0 / a));
    double tau = (t - t0) * sign(wx0 * wz0) * b * sqrt( (Inertia[8] - Inertia[4]) * (Inertia[4] - Inertia[0]) / (Inertia[0] * Inertia[8]) );


    return ( sign(wz0) * sqrt( 2 * T * (D - Inertia[0]) / (Inertia[8] * (Inertia[8] - Inertia[0])) ) * jacobi_cn(k, tau + x0) );
}

CubeSpline Pitch(const vector<vector<double> > &args, double t0, double h)
{
    unsigned k = static_cast<int>(t0 / h) + 1;
    vector<double> X, Y;

    double L = KineticIntegral(args[6][k], args[7][k], args[8][k]);

    double t = t0;
    for(unsigned i = k; i < args[13].size(); i++)
    {
        X.push_back(t);
        Y.push_back(acos(Inertia[0] / L * args[6][i]));
        t += h;
    }

    return CubeSpline(X, Y);
}

CubeSpline Yaw(const vector<vector<double> > &args, double t0, double h)
{
    unsigned k = static_cast<int>(t0 / h) + 1;
    vector<double> X, Y;
    vector<double> temp(args.size());

    double L = KineticIntegral(args[6][k], args[7][k], args[8][k]);

    for (unsigned i = k; i < args[0].size(); i++)
    {
        for (unsigned j = 0; j < args.size(); j++)
            temp[j] = args[j][i];

        X.push_back(temp[temp.size() - 1]);
        double dphi = Inertia[8] / Inertia[4] * (pow(Inertia[4] * temp[7], 2) / (pow(Inertia[4] * temp[7], 2) + pow(Inertia[8] * temp[8], 2))) *
                (dwz(temp) * temp[7] - dwy(temp) * temp[8]) / pow(temp[7], 2);
        Y.push_back(L / Inertia[0] * (1.0 - dphi / temp[6]));
    }

    return CubeSpline(X, Y);
}
