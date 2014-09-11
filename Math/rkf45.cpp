#include "rkf45.h"
#include <cmath>
#include <fstream>
#include <QTimer>

using namespace std;

namespace NMath {

double Func(const double *wcoeff, const double *ki, unsigned short num)
{
    double result = 0.0;

    for (unsigned short i = 0; i < num; i++)
        result += wcoeff[i] * ki[i];

    return result;
}

double ki(const vector<double> &vars, unsigned short i, const double &h, const vector<pfunc> &RightFuncs, unsigned j)
{
    switch (i) {
    case 1:
        return RightFuncs[j](vars);
    case 2:
    {
        vector<double> temp_vars = vars;
        for (unsigned i = 0; i < vars.size() - 1; i++)
            temp_vars[i] += 0.5 * h * ki(vars, 1, h, RightFuncs, i);
        temp_vars[vars.size() - 1] += 0.5 * h;

        return RightFuncs[j](temp_vars);
    }
    case 3:
    {
        vector<double> temp_vars = vars;
        for (unsigned i = 0; i < vars.size() - 1; i++)
            temp_vars[i] += 0.5 * h * ki(vars, 2, h, RightFuncs, i);
        temp_vars[vars.size() - 1] += 0.5 * h;

        return RightFuncs[j](temp_vars);
    }
    case 4:
    {
        vector<double> temp_vars = vars;
        for (unsigned i = 0; i < vars.size() - 1; i++)
            temp_vars[i] += h * ki(vars, 3, h, RightFuncs, i);
        temp_vars[vars.size() - 1] += h;

        return RightFuncs[j](temp_vars);
    }

    }
}

void rkf45(vector< vector<double> > &solve, const vector<double> &vars, const double &h, const double &x_beg, const double &x_end, const vector<pfunc> &RightFuncs)
{
    double kcoeff[4];
    unsigned num = static_cast<unsigned>( (x_end - x_beg) / h ) + 1;
    vector<double> temp(vars.size());
    const double w[] = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};

    for (unsigned j = 0; j < solve.size(); j++)
    {
        solve[j].resize(num + 1);
        solve[j][0] = vars[j];
    }

    for (unsigned i = 1; i <= num; i++)
    {
        for (unsigned k = 0; k < solve.size(); k++)
            temp[k] = solve[k][i - 1];
        for (unsigned j = 0; j < solve.size() - 1; j++)
        {
            for (unsigned short k = 1; k <= 4; k++)
                kcoeff[k - 1] = ki(temp, k, h, RightFuncs, j);

            solve[j][i] = solve[j][i - 1] + h * Func(w, kcoeff, 4);
        }
        solve[solve.size() - 1][i] = solve[solve.size() - 1][i - 1] + h;
    }
}

}
