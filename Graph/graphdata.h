#ifndef GRAPHDATA_H
#define GRAPHDATA_H

#include <vector>

class GraphData
{
    std::vector<double> x;
    std::vector<double> y;
    unsigned nPoints;

    double coeffX, coeffY;
    int w, h;
    double aMin, aMax;
    double fMin, fMax;

private://methods
    double FindMinimum(const std::vector<double> &arg);
    double FindMaximum(const std::vector<double> &arg);
public:
    explicit GraphData(const std::vector<double> &X, const std::vector<double> &Y,
                       int w, int h);
    double getArg(unsigned i) const;
    double getVal(unsigned i) const;

    double getArgMax() const;
    double getArgMin() const;

    double getValMax() const;
    double getValMin() const;

    double getCoeffX() const;
    double getCoeffY() const;

    int getWidth() const;
    int getHeight() const;

    void setWidthHeight(int width, int height);

    unsigned numPoints() const;

    ~GraphData();
};

#endif // GRAPHDATA_H
