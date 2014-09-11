#include "cubespline.h"

#include <cmath>
using namespace std;

CubeSpline::CubeSpline(const vector<double> &x, const vector<double> &y)
{
    X = x;
    Y = y;
    h = abs(X[0] - X[1]);

    double a = h / 6.0;
    double b = 2 * h / 3.0;
    double c = h / 6.0;
    q.resize(x.size());
    vector<double> d(x.size());
    vector<double> A(x.size());
    vector<double> B(x.size());

    for (unsigned i = 0; i < q.size(); i++)
        d[i] = q[i] = A[i] = B[i] = 0;

    for (unsigned i = 1; i < d.size() - 1; i++)
        d[i] = (Y[i + 1] - Y[i]) / h - (Y[i] - Y[i - 1]) / h;

    A[1] = - c / b;
    B[1] = d[1] / b;
    for (unsigned i = 2; i < A.size() - 1; i++)
    {
        A[i] = - c / (b + a * A[i - 1]);
        B[i] = (d[i] - a * B[i - 1]) / (b + a * A[i - 1]);
    }

    for (unsigned i = q.size() - 2; i > 0; i--)
    {
        q[i] = q[i + 1] * A[i] + B[i];
    }

}

double CubeSpline::GetValueOfFunction(const double &x) const
{
    unsigned k = 0;

    for (unsigned i = 1; i < X.size(); i++)
    {
        if (x <= X[i] && x >= X[i - 1])
        {
            k = i;
            break;
        }
    }
    if (!k)
        return 0.0;

    double res;
    res = q[k - 1] * (X[k] - x) * (X[k] - x) * (X[k] - x) / (6 * h) +
            q[k] * (x - X[k - 1]) * (x - X[k - 1]) * (x - X[k - 1]) / (6 * h) +
            (Y[k - 1] / h - q[k - 1] * h / 6) * (X[k] - x) +
            (Y[k] / h - q[k] * h / 6) * (x - X[k - 1]);

    return res;
}

void CubeSpline::operator=(const CubeSpline &spline)
{
    X = spline.X;
    Y = spline.Y;
    h = spline.h;
    q = spline.q;
}

CubeSpline::~CubeSpline()
{
    q.clear();
    X.clear();
    Y.clear();
}
