#include "graphdata.h"
#include <cmath>

using namespace std;

GraphData::GraphData(const vector<double> &X, const vector<double> &Y,
                     int width, int height)
{
    nPoints = X.size();
    x = X;
    y = Y;

    w = width; h = height;

    fMin = FindMinimum(y); /*fMin = (fMin < 0) ? (floor(fMin) - 1.0) : (ceil(fMin) + 1.0);*/
    fMax = FindMaximum(y);
    aMin = FindMinimum(x); /*aMin = (aMin < 0) ? (floor(aMin - 1.0)) : (ceil(aMin) + 1.0);*/
    aMax = FindMaximum(x);

    coeffX = static_cast<double>(w) / abs(aMax - aMin);
    coeffY = static_cast<double>(h) / abs(fMax - fMin);
}

void GraphData::setWidthHeight(int width, int height)
{
    w = width; h = height;
    coeffX = static_cast<double>(w) / abs(aMax - aMin);
    coeffY = static_cast<double>(h) / abs(fMax - fMin);
}

GraphData::~GraphData()
{
}

double GraphData::getArg(unsigned i) const
{
    if (i < nPoints)
        return (x[i] - aMin) * coeffX;
    else
        return (x[nPoints - 1] - aMin) * coeffX;
}

double GraphData::getVal(unsigned i) const
{
    if (i < nPoints)
        return ( h - (y[i] - fMin) * coeffY );
    else
        return (h - (y[nPoints - 1] - fMin) * coeffY);
}

double GraphData::FindMinimum(const vector<double> &arg)
{
    double temp = arg[0];

    for (unsigned i = 0; i < nPoints; i++)
        if (temp > arg[i])
            temp = arg[i];

    return temp;
}

double GraphData::FindMaximum(const vector<double> &arg)
{
    double temp = arg[0];

    for (unsigned i = 0; i < nPoints; i++)
        if (temp < arg[i])
            temp = arg[i];

    return temp;
}

double GraphData::getArgMax() const
{
    return aMax;
}

double GraphData::getArgMin() const
{
    return aMin;
}

double GraphData::getValMax() const
{
    return fMax;
}

double GraphData::getValMin() const
{
    return fMin;
}

double GraphData::getCoeffX() const
{
    return coeffX;
}

double GraphData::getCoeffY() const
{
    return coeffY;
}

int GraphData::getWidth() const
{
    return w;
}

int GraphData::getHeight() const
{
    return h;
}

unsigned GraphData::numPoints() const
{
    return nPoints;
}
