#include <QApplication>
#include <QtWidgets>
#include "RightFuncs.h"
#include "Math/allmath.h"
#include <fstream>
#include <vector>
#include <iomanip>
#include "Graph/graphdata.h"
#include "Graph/graphicwidget.h"
using namespace std;
using namespace NMath;

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    app.setStyle("fusion");

    alpha = TO_RADS(20.0, 0.0, 0.0);
    beta = TO_RADS(47.0, 0.0, 0.0);
    gamma = TO_RADS(4.0, 19.0, 0.0);

    vector<double> X(11), Y(11);

    X[0] = 0.0;       Y[0] = 3000.0;
    X[1] = 5.0;       Y[1] = 2727.56981;
    X[2] = 10.0;      Y[2] = 2109.77421;
    X[3] = 15.0;      Y[3] = 1145.66250;
    X[4] = 20.0;      Y[4] = 784.87819;
    X[5] = 25.0;      Y[5] = 567.91443;
    X[6] = 30.0;      Y[6] = 401.27748;
    X[7] = 35.0;      Y[7] = 259.83392;
    X[8] = 40.0;      Y[8] = 160.75460;
    X[9] = 45.0;      Y[9] = 69.57813;
    X[10] = 50.0;      Y[10] = 0.0;

    spForce = new CubeSpline(X, Y);

    mass = 3000.0;
    radius = vector3d(0.9, 1.33 * cos(beta), 1.33 * sin(beta));

    Inertia = matrix3x3(3000.0, 15000.0, 15500.0);

    vector< vector<double> > solves(14);
    vector<double> beg_cond(14);
    for (unsigned i = 0; i < beg_cond.size(); i++)
        beg_cond[i] = 0.0;
    beg_cond[9] = 1.0;
    vector<pfunc> funcs;
    funcs.push_back(dx);
    funcs.push_back(dy);
    funcs.push_back(dz);

    funcs.push_back(ddx);
    funcs.push_back(ddy);
    funcs.push_back(ddz);

    funcs.push_back(dwx);
    funcs.push_back(dwy);
    funcs.push_back(dwz);

    funcs.push_back(dlambda0);
    funcs.push_back(dlambda1);
    funcs.push_back(dlambda2);
    funcs.push_back(dlambda3);

    rkf45(solves, beg_cond, 0.1, 0.0, 60.0, funcs);

    CubeSpline spPitch, spYaw;
    spPitch = Pitch(solves, 50.0, 0.1);
    spYaw = Yaw(solves, 50.0, 0.1);
    vector<double> T, P, Ya, anWx, anWy, anWz, numWx, numWy, numWz, numT;
    double wx0 = solves[6][static_cast<int>(50 / 0.1) + 1];
    double wy0 = solves[7][static_cast<int>(50 / 0.1) + 1];
    double wz0 = solves[8][static_cast<int>(50 / 0.1) + 1];
    for (double t = 50.000001; t <= 60.0; t += 0.01)
    {
        T.push_back(t);
        P.push_back(spPitch.GetValueOfFunction(t));
        Ya.push_back(spYaw.GetValueOfFunction(t));
        anWx.push_back(AnalyticOmegaX(wx0, wy0, wz0, 50.0, t));
        anWy.push_back(AnalyticOmegaY(wx0, wy0, wz0, 50.0, t));
        anWz.push_back(AnalyticOmegaZ(wx0, wy0, wz0, 50.0, t));
    }
    for (int i = static_cast<int>(50 / 0.1) + 1; i < solves[13].size(); ++i) {
        numT.push_back(solves[13][i]);
        numWx.push_back(solves[6][i]);
        numWy.push_back(solves[7][i]);
        numWz.push_back(solves[8][i]);
    }

    QTabWidget *pTabWidget = new QTabWidget;
    pTabWidget->resize(600, 600);
    QGridLayout *pGridLayouts[4];
    QWidget *pWidgets[4];
    GraphicWidget *pGraphWidgets[13];

    GraphData gd[] = {
        GraphData(solves[13], solves[0], 100, 100), GraphData(solves[13], solves[1], 100, 100), GraphData(solves[13], solves[2], 100, 100),
        GraphData(solves[13], solves[3], 100, 100), GraphData(solves[13], solves[4], 100, 100), GraphData(solves[13], solves[5], 100, 100),
        GraphData(solves[13], solves[6], 100, 100), GraphData(solves[13], solves[7], 100, 100), GraphData(solves[13], solves[8], 100, 100),
        GraphData(solves[13], solves[9], 100, 100), GraphData(solves[13], solves[10], 100, 100), GraphData(solves[13], solves[11], 100, 100), GraphData(solves[13], solves[12], 100, 100)
                     };

    for (int i = 0; i < 13; i++)
    {
        pGraphWidgets[i] = new GraphicWidget;
        pGraphWidgets[i]->setColor(Qt::blue);
        pGraphWidgets[i]->setLineWidth(2);
        pGraphWidgets[i]->AddGraphic(gd[i]);
    }

    for (int i = 0; i < 4; i++)
        pGridLayouts[i] = new QGridLayout;
    pGridLayouts[0]->addWidget(pGraphWidgets[0], 0, 0);
    pGridLayouts[0]->addWidget(pGraphWidgets[1], 0, 1);
    pGridLayouts[0]->addWidget(pGraphWidgets[2], 1, 0);

    pGridLayouts[1]->addWidget(pGraphWidgets[3], 0, 0);
    pGridLayouts[1]->addWidget(pGraphWidgets[4], 0, 1);
    pGridLayouts[1]->addWidget(pGraphWidgets[5], 1, 0);

    pGridLayouts[2]->addWidget(pGraphWidgets[6], 0, 0);
    pGridLayouts[2]->addWidget(pGraphWidgets[7], 0, 1);
    pGridLayouts[2]->addWidget(pGraphWidgets[8], 1, 0);

    pGridLayouts[3]->addWidget(pGraphWidgets[9], 0, 0);
    pGridLayouts[3]->addWidget(pGraphWidgets[10], 0, 1);
    pGridLayouts[3]->addWidget(pGraphWidgets[11], 1, 0);
    pGridLayouts[3]->addWidget(pGraphWidgets[12], 1, 1);

    for (int i = 0; i < 4; i++)
    {
        pWidgets[i] = new QWidget;
        pWidgets[i]->setLayout(pGridLayouts[i]);
    }

    GraphData pitchData(T, P, 100, 100);
    GraphData yawData(T, Ya, 100, 100);
    GraphicWidget *pPitchGraph = new GraphicWidget;
    GraphicWidget *pYawGraph = new GraphicWidget;

    pPitchGraph->setColor(Qt::red);
    pPitchGraph->setLineWidth(2);
    pPitchGraph->AddGraphic(pitchData);

    pYawGraph->setColor(Qt::red);
    pYawGraph->setLineWidth(2);
    pYawGraph->AddGraphic(yawData);

    QGridLayout *pGridLayoutsPitchYaw = new QGridLayout;
    pGridLayoutsPitchYaw->addWidget(pPitchGraph, 0, 0);
    pGridLayoutsPitchYaw->addWidget(pYawGraph, 0, 1);
    QWidget *pWidgetsPitchYaw = new QWidget;
    pWidgetsPitchYaw->setLayout(pGridLayoutsPitchYaw);


    GraphData analytGD[] = {
        GraphData(T, anWx, 100, 100),
        GraphData(T, anWy, 100, 100),
        GraphData(T, anWz, 100, 100),
        GraphData(numT, numWx, 100, 100),
        GraphData(numT, numWy, 100, 100),
        GraphData(numT, numWz, 100, 100)
                            };
    GraphicWidget *pAnalyticGraph[3];
    for (int i = 0; i < 3; i++)
    {
        pAnalyticGraph[i] = new GraphicWidget;
        pAnalyticGraph[i]->setColor(Qt::green);
        pAnalyticGraph[i]->setLineWidth(2);
        pAnalyticGraph[i]->AddGraphic(analytGD[i]);
    }
    for (int i = 3; i < 6; i++)
    {
        pAnalyticGraph[i - 3]->AddGraphic(analytGD[i]);
    }
    QGridLayout *pAnalyticLayout = new QGridLayout;
    pAnalyticLayout->addWidget(pAnalyticGraph[0], 0, 0);
    pAnalyticLayout->addWidget(pAnalyticGraph[1], 0, 1);
    pAnalyticLayout->addWidget(pAnalyticGraph[2], 1, 0);
    QWidget *pWidgetAnalytic = new QWidget;
    pWidgetAnalytic->setLayout(pAnalyticLayout);

    pTabWidget->addTab(pWidgets[0], "Coordinates center of mass");
    pTabWidget->addTab(pWidgets[1], "Velocities center of mass");
    pTabWidget->addTab(pWidgets[2], "Angular velocities");
    pTabWidget->addTab(pWidgets[3], "Quaternion params");
    pTabWidget->addTab(pWidgetsPitchYaw, "Pitch and Precession");
    pTabWidget->addTab(pWidgetAnalytic, "Analytic solves");

    pTabWidget->show();

    delete spForce;

    return app.exec();
}
