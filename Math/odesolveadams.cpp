#include "odesolveadams.h"

#include <fstream>
#include <ctime>

vector3d Multiply(const matrix3x3 &mat, const vector3d &vec)
{
    vector3d temp;

    temp.SetX(mat[0] * vec.GetX() + mat[1] * vec.GetY() + mat[2] * vec.GetZ());
    temp.SetY(mat[3] * vec.GetX() + mat[4] * vec.GetY() + mat[5] * vec.GetZ());
    temp.SetZ(mat[6] * vec.GetX() + mat[7] * vec.GetY() + mat[8] * vec.GetZ());

    return temp;
}

ODESolveAdams::ODESolveAdams(const double *Y0, const double &x0, const double &xk, const double &st, int num)
{
    //Init constants
    NumOfEquations = num;
    step = st;
    SizeofArrayPoints = unsigned(float((xk - x0) / step));
    SizeofArrayPoints = SizeofArrayPoints + 1;

    //Init vector of solves
    X = new double[SizeofArrayPoints];
    vecForce = new vector3d[SizeofArrayPoints];
    vecForce1 = new vector3d[SizeofArrayPoints];
    Y = new double*[NumOfEquations];
    for (int i = 0; i < NumOfEquations; i++)
        Y[i] = new double[SizeofArrayPoints];
    for (int i = 0; i < NumOfEquations; i++)
        Y[i][0] = Y0[i];
    for (unsigned i = 0; i < SizeofArrayPoints; i++)
        X[i] = i * step;

    dOmega = new double*[3];
    for (int i = 0; i < 3; i++)
        dOmega[i] = new double[SizeofArrayPoints];
    for (int i = 0; i < 3; i++)
        dOmega[i][0] = 0;

    //Quaternion parameters
    lambda.lambda0 = Y0[3];
    lambda.lambda1 = Y0[4];
    lambda.lambda2 = Y0[5];
    lambda.lambda3 = Y0[6];

    //make rotate matrix
    double m[9];
    m[0] = 2 * (lambda.lambda0 * lambda.lambda0 + lambda.lambda1 * lambda.lambda1) - 1;
    m[4] = 2 * (lambda.lambda0 * lambda.lambda0 + lambda.lambda2 * lambda.lambda2) - 1;
    m[8] = 2 * (lambda.lambda0 * lambda.lambda0 + lambda.lambda3 * lambda.lambda3) - 1;

    m[1] = 2 * (lambda.lambda1 * lambda.lambda2 - lambda.lambda0 * lambda.lambda3);
    m[2] = 2 * (lambda.lambda1 * lambda.lambda3 + lambda.lambda0 * lambda.lambda2);

    m[3] = 2 * (lambda.lambda1 * lambda.lambda2 + lambda.lambda0 * lambda.lambda3);
    m[5] = 2 * (lambda.lambda2 * lambda.lambda3 - lambda.lambda0 * lambda.lambda1);

    m[6] = 2 * (lambda.lambda1 * lambda.lambda3 - lambda.lambda0 * lambda.lambda2);
    m[7] = 2 * (lambda.lambda2 * lambda.lambda3 + lambda.lambda0 * lambda.lambda1);

    matRotate = m;

    //Center of mass of the system
    centerofmass.SetX(0.9);
    centerofmass.SetY(1.33 * cos(TO_RADS(47, 0, 0)));
    centerofmass.SetZ(1.33 * sin(TO_RADS(47, 0, 0)));

    //Force projection on related coordinate system
    force_projX = cos(TO_RADS(4, 19, 0)) * cos(TO_RADS(20, 0, 0));
    force_projY = sin(TO_RADS(20, 0, 0)) * cos(TO_RADS(4, 19, 0)) * cos(TO_RADS(47, 0, 0)) + sin(TO_RADS(4, 19, 0)) * sin(TO_RADS(47, 0, 0));
    force_projZ = sin(TO_RADS(20, 0, 0)) * cos(TO_RADS(4, 19, 0)) * sin(TO_RADS(47, 0, 0)) - sin(TO_RADS(4, 19, 0)) * cos(TO_RADS(47, 0, 0));
}

void ODESolveAdams::DeleteData()
{
    delete [] X;
    X = 0;

    step = 0;
    SizeofArrayPoints = 0;

    delete [] vecForce;
    vecForce = 0;
    delete [] vecForce1;
    vecForce1 = 0;

    for (int i = 0; i < NumOfEquations; i++)
    {
        delete [] Y[i];
        Y[i] = 0;
    }
    delete [] Y;
    Y = 0;
    delete splineForce;
}

ODESolveAdams::~ODESolveAdams()
{
    DeleteData();
}

void ODESolveAdams::SetForce(const TableOfPoints &table)
{
    splineForce = new CubeSpline(table);
    splineForce->CreateSpline();
}

void ODESolveAdams::SetTensorOfInertia(const double *tensor)
{
    for (int i = 0; i < 9; i++)
        tensorOfInertia[i] = tensor[i];
}

void ODESolveAdams::SetTensorOfInertia(const matrix3x3 &tensor)
{
    for (int i = 0; i < 9; i++)
        tensorOfInertia[i] = tensor[i];
}

void ODESolveAdams::SetNewInitialCondition(const double *Y0, const double &x0, const double &xk, const double &st)
{
    //Delete old data
    DeleteData();

    //rewrite all
    step = st;
    SizeofArrayPoints = int((xk - x0) / step);
    X = new double[SizeofArrayPoints];
    vecForce = new vector3d[SizeofArrayPoints];
    vecForce1 = new vector3d[SizeofArrayPoints];
    Y = new double*[NumOfEquations];
    for (int i = 0; i < NumOfEquations; i++)
        Y[i] = new double[SizeofArrayPoints];

    lambda.lambda0 = Y0[3];
    lambda.lambda1 = Y0[4];
    lambda.lambda2 = Y0[5];
    lambda.lambda3 = Y0[6];

    double m[9];
    m[0] = 2 * (lambda.lambda0 * lambda.lambda0 + lambda.lambda1 * lambda.lambda1) - 1;
    m[4] = 2 * (lambda.lambda0 * lambda.lambda0 + lambda.lambda2 * lambda.lambda2) - 1;
    m[8] = 2 * (lambda.lambda0 * lambda.lambda0 + lambda.lambda3 * lambda.lambda3) - 1;

    m[1] = 2 * (lambda.lambda1 * lambda.lambda2 - lambda.lambda0 * lambda.lambda3);
    m[2] = 2 * (lambda.lambda1 * lambda.lambda3 + lambda.lambda0 * lambda.lambda2);

    m[3] = 2 * (lambda.lambda1 * lambda.lambda2 + lambda.lambda0 * lambda.lambda3);
    m[5] = 2 * (lambda.lambda2 * lambda.lambda3 - lambda.lambda0 * lambda.lambda1);

    m[6] = 2 * (lambda.lambda1 * lambda.lambda3 - lambda.lambda0 * lambda.lambda2);
    m[7] = 2 * (lambda.lambda2 * lambda.lambda3 + lambda.lambda0 * lambda.lambda1);

    matRotate = m;

    for (int i = 0; i < NumOfEquations; i++)
        Y[i][0] = Y0[i];
}

void ODESolveAdams::SetMass(const double &m)
{
    mass = m;
}

void ODESolveAdams::SolveODE()
{
    std::ofstream log("log.txt");
    log << "Welcome to log file" << "\n";
    log << "Here we output computation strings\n";
    time_t tme = time(NULL);
    char *tm_str = asctime(localtime(&tme));
    log << "[" << tm_str << "]" << " : Iteration " << 1 << " |Lambda| is " <<
           sqrt(lambda.lambda0 * lambda.lambda0 + lambda.lambda1 * lambda.lambda1 + lambda.lambda2 * lambda.lambda2 + lambda.lambda3 * lambda.lambda3) << std::endl;

    //function in right side of vector dw = {dwx, dwy, dwz}
    vector3d *Fw = new vector3d[SizeofArrayPoints];

    //if k = 0 - Euler method y(i+1) = y(i) + h f(i)
    for (int i = 0; i < 3; i++)
    {
        //x, y, z (dx/dt = vx, etc)
        Y[i][1] = Y[i][0] + step * Y[7 + i][0];
    }

    //l0, l1, l2, l3 (dl0/dt = -wx*l1.., etc) Quaternions parameters
    Y[3][1] = Y[3][0] + step * 0.5 * (- Y[10][0] * Y[4][0] - Y[11][0] * Y[5][0] - Y[12][0] * Y[6][0]);
    Y[4][1] = Y[4][0] + step * 0.5 * (Y[10][0] * Y[3][0] + Y[12][0] * Y[5][0] - Y[11][0] * Y[6][0]);
    Y[5][1] = Y[5][0] + step * 0.5 * (Y[11][0] * Y[3][0] - Y[12][0] * Y[4][0] + Y[10][0] * Y[6][0]);
    Y[6][1] = Y[6][0] + step * 0.5 * (Y[12][0] * Y[3][0] + Y[11][0] * Y[4][0] - Y[10][0] * Y[5][0]);

    //Force in related frame
    vecForce1[0].SetX( - splineForce->GetValueOfFunction(X[0]) * force_projX);
    vecForce1[0].SetY(- splineForce->GetValueOfFunction(X[0]) * force_projY);
    vecForce1[0].SetZ( - splineForce->GetValueOfFunction(X[0]) * force_projZ);

    //matRotate.Transpose();
    vecForce[0] = Multiply(matRotate, vecForce1[0]);
    //vx, vy, vz (dvx/dt = Fx / m, etc)
    Y[7][1] = Y[7][0] + step * vecForce[0].GetX() / mass;
    Y[8][1] = Y[8][0] + step * vecForce[0].GetY() / mass;
    Y[9][1] = Y[9][0] + step * vecForce[0].GetZ() / mass;

    //wx, wy, wz (dwx/dt = f(wx, wy, wz, Mx), etc)
    vector3d temp_cross = centerofmass.cross(vecForce1[0]);
    Fw[0].SetX((temp_cross.GetX() +
            (tensorOfInertia[4] - tensorOfInertia[8]) * Y[11][0] * Y[12][0]) / tensorOfInertia[0]);
    Fw[0].SetY((temp_cross.GetY() +
            (tensorOfInertia[8] - tensorOfInertia[0]) * Y[10][0] * Y[12][0]) / tensorOfInertia[4]);
    Fw[0].SetZ((temp_cross.GetZ() +
            (tensorOfInertia[0] - tensorOfInertia[4]) * Y[11][0] * Y[10][0]) / tensorOfInertia[8]);
    dOmega[0][1] = Fw[0].GetX();
    dOmega[1][1] = Fw[0].GetY();
    dOmega[2][1] = Fw[0].GetZ();
    Y[10][1] = Y[10][0] + step * Fw[0].GetX();
    Y[11][1] = Y[11][0] + step * Fw[0].GetY();
    Y[12][1] = Y[12][0] + step * Fw[0].GetZ();

//------------------------------------------------------------------------------------------------------//

    //if k = 1 y(i+1) = y(i) + h (3/2f(i)-1/2f(i-1))
    lambda.lambda0 = Y[3][1];
    lambda.lambda1 = Y[4][1];
    lambda.lambda2 = Y[5][1];
    lambda.lambda3 = Y[6][1];

    log << "[" << tm_str << "]" << " : Iteration " << 2 << "|Lambda| is " <<
           sqrt(lambda.lambda0 * lambda.lambda0 + lambda.lambda1 * lambda.lambda1 + lambda.lambda2 * lambda.lambda2 + lambda.lambda3 * lambda.lambda3) << std::endl;

    double m[9];
    m[0] = 2 * (lambda.lambda0 * lambda.lambda0 + lambda.lambda1 * lambda.lambda1) - 1;
    m[4] = 2 * (lambda.lambda0 * lambda.lambda0 + lambda.lambda2 * lambda.lambda2) - 1;
    m[8] = 2 * (lambda.lambda0 * lambda.lambda0 + lambda.lambda3 * lambda.lambda3) - 1;

    m[1] = 2 * (lambda.lambda1 * lambda.lambda2 - lambda.lambda0 * lambda.lambda3);
    m[2] = 2 * (lambda.lambda1 * lambda.lambda3 + lambda.lambda0 * lambda.lambda2);

    m[3] = 2 * (lambda.lambda1 * lambda.lambda2 + lambda.lambda0 * lambda.lambda3);
    m[5] = 2 * (lambda.lambda2 * lambda.lambda3 - lambda.lambda0 * lambda.lambda1);

    m[6] = 2 * (lambda.lambda1 * lambda.lambda3 - lambda.lambda0 * lambda.lambda2);
    m[7] = 2 * (lambda.lambda2 * lambda.lambda3 + lambda.lambda0 * lambda.lambda1);

    matRotate = m;

    vecForce1[1].SetX( - splineForce->GetValueOfFunction(X[1]) * force_projX);
    vecForce1[1].SetY(- splineForce->GetValueOfFunction(X[1]) * force_projY);
    vecForce1[1].SetZ( - splineForce->GetValueOfFunction(X[1]) * force_projZ);

    //matRotate.Transpose();
    vecForce[1] = Multiply(matRotate, vecForce1[1]);

    for (int i = 0; i < 3; i++)
    {
        //x, y, z (dx/dt = vx, etc)
        Y[i][2] = Y[i][1] + step * (1.5 * Y[7 + i][1] - 0.5 * Y[7 + i][0]);
    }

    //vx, vy, vz (dvx/dt = Fx / m, etc)
    Y[7][2] = Y[7][1] + step * (1.5 * vecForce[1].GetX() / mass - 0.5 * vecForce[0].GetX() / mass);
    Y[8][2] = Y[8][1] + step * (1.5 * vecForce[1].GetY() / mass - 0.5 * vecForce[0].GetY() / mass);
    Y[9][2] = Y[9][1] + step * (1.5 * vecForce[1].GetZ() / mass - 0.5 * vecForce[0].GetZ() / mass);

    //wx, wy, wz (dwx/dt = f(wx, wy, wz, Mx), etc)
    temp_cross = centerofmass.cross(vecForce1[1]);
    Fw[1].SetX((temp_cross.GetX() +
            (tensorOfInertia[4] - tensorOfInertia[8]) * Y[11][1] * Y[12][1]) / tensorOfInertia[0]);
    Fw[1].SetY((temp_cross.GetY() +
            (tensorOfInertia[8] - tensorOfInertia[0]) * Y[10][1] * Y[12][1]) / tensorOfInertia[4]);
    Fw[1].SetZ((temp_cross.GetZ() +
            (tensorOfInertia[0] - tensorOfInertia[4]) * Y[11][1] * Y[10][1]) / tensorOfInertia[8]);
    dOmega[0][2] = Fw[1].GetX();
    dOmega[1][2] = Fw[1].GetY();
    dOmega[2][2] = Fw[1].GetZ();
    Y[10][2] = Y[10][1] + step * (1.5 * Fw[1].GetX() - 0.5 * Fw[0].GetX());
    Y[11][2] = Y[11][1] + step * (1.5 * Fw[1].GetY() - 0.5 * Fw[0].GetY());
    Y[12][2] = Y[12][1] + step * (1.5 * Fw[1].GetZ() - 0.5 * Fw[0].GetZ());

    //l0, l1, l2, l3 (dl0/dt = -wx*l1.., etc) Quaternions parameters
    Y[3][2] = Y[3][1] + step * 0.5 * (1.5 *(- Y[10][1] * Y[4][1] - Y[11][1] * Y[5][1] - Y[12][1] * Y[6][1]) -
              0.5 * (- Y[10][0] * Y[4][0] - Y[11][0] * Y[5][0] - Y[12][0] * Y[6][0]));

    Y[4][2] = Y[4][1] + step * 0.5 * (1.5 * (Y[10][1] * Y[3][1] + Y[12][1] * Y[5][1] - Y[11][1] * Y[6][1]) -
              0.5 * (Y[10][0] * Y[3][0] + Y[12][0] * Y[5][0] - Y[11][0] * Y[6][0]));

    Y[5][2] = Y[5][1] + step * 0.5 * (1.5 * (Y[11][1] * Y[3][1] - Y[12][1] * Y[4][1] + Y[10][1] * Y[6][1]) -
            0.5 * (Y[11][0] * Y[3][0] - Y[12][0] * Y[4][0] + Y[10][0] * Y[6][0]));

    Y[6][2] = Y[6][1] + step * 0.5 * (1.5 * (Y[12][1] * Y[3][1] + Y[11][1] * Y[4][1] - Y[10][1] * Y[5][1]) -
              0.5 * (Y[12][0] * Y[3][0] + Y[11][0] * Y[4][0] - Y[10][0] * Y[5][0]));

    //-------------------------------------------------------------------------------------------------------//

    //if k = 2 y(i+1) = y(i) + h / 12 (23f(i)-16f(i-1) + 5f(i-2))
    lambda.lambda0 = Y[3][2];
    lambda.lambda1 = Y[4][2];
    lambda.lambda2 = Y[5][2];
    lambda.lambda3 = Y[6][2];

    log << "[" << tm_str << "]" << " : Iteration " << 3 << "|Lambda| is " <<
           sqrt(lambda.lambda0 * lambda.lambda0 + lambda.lambda1 * lambda.lambda1 + lambda.lambda2 * lambda.lambda2 + lambda.lambda3 * lambda.lambda3) << std::endl;


    m[0] = 2 * (lambda.lambda0 * lambda.lambda0 + lambda.lambda1 * lambda.lambda1) - 1;
    m[4] = 2 * (lambda.lambda0 * lambda.lambda0 + lambda.lambda2 * lambda.lambda2) - 1;
    m[8] = 2 * (lambda.lambda0 * lambda.lambda0 + lambda.lambda3 * lambda.lambda3) - 1;

    m[1] = 2 * (lambda.lambda1 * lambda.lambda2 - lambda.lambda0 * lambda.lambda3);
    m[2] = 2 * (lambda.lambda1 * lambda.lambda3 + lambda.lambda0 * lambda.lambda2);

    m[3] = 2 * (lambda.lambda1 * lambda.lambda2 + lambda.lambda0 * lambda.lambda3);
    m[5] = 2 * (lambda.lambda2 * lambda.lambda3 - lambda.lambda0 * lambda.lambda1);

    m[6] = 2 * (lambda.lambda1 * lambda.lambda3 - lambda.lambda0 * lambda.lambda2);
    m[7] = 2 * (lambda.lambda2 * lambda.lambda3 + lambda.lambda0 * lambda.lambda1);

    matRotate = m;

    vecForce1[2].SetX( - splineForce->GetValueOfFunction(X[2]) * force_projX);
    vecForce1[2].SetY( - splineForce->GetValueOfFunction(X[2]) * force_projY);
    vecForce1[2].SetZ( - splineForce->GetValueOfFunction(X[2]) * force_projZ);

    //matRotate.Transpose();
    vecForce[2] = Multiply(matRotate, vecForce1[2]);

    for (int i = 0; i < 3; i++)
    {
        //x, y, z (dx/dt = vx, etc)
        Y[i][3] = Y[i][2] + step / 12.0 * (23.0 * Y[7 + i][2] - 16.0 * Y[7 + i][1] + 5.0 * Y[7 + i][0]);
    }

    //vx, vy, vz (dvx/dt = Fx / m, etc)
    Y[7][3] = Y[7][2] + step / 12.0 * (23.0 * vecForce[2].GetX() - 16.0 * vecForce[1].GetX() + 5.0 * vecForce[0].GetX()) / mass;
    Y[8][3] = Y[8][2] + step / 12.0 * (23.0 * vecForce[2].GetY() - 16.0 * vecForce[1].GetY() + 5.0 * vecForce[0].GetY()) / mass;
    Y[9][3] = Y[9][2] + step / 12.0 * (23.0 * vecForce[2].GetZ() - 16.0 * vecForce[1].GetZ() + 5.0 * vecForce[0].GetZ()) / mass;

    //wx, wy, wz (dwx/dt = f(wx, wy, wz, Mx), etc)
    temp_cross = centerofmass.cross(vecForce1[2]);
    Fw[2].SetX((temp_cross.GetX() +
            (tensorOfInertia[4] - tensorOfInertia[8]) * Y[11][2] * Y[12][2]) / tensorOfInertia[0]);
    Fw[2].SetY((temp_cross.GetY() +
            (tensorOfInertia[8] - tensorOfInertia[0]) * Y[10][2] * Y[12][2]) / tensorOfInertia[4]);
    Fw[2].SetZ((temp_cross.GetZ() +
            (tensorOfInertia[0] - tensorOfInertia[4]) * Y[11][2] * Y[10][2]) / tensorOfInertia[8]);
    dOmega[0][3] = Fw[2].GetX();
    dOmega[1][3] = Fw[2].GetY();
    dOmega[2][3] = Fw[2].GetZ();
    Y[10][3] = Y[10][2] + step / 12.0 * (23.0 * Fw[2].GetX() - 16.0 * Fw[1].GetX() + 5.0 * Fw[0].GetX());
    Y[11][3] = Y[11][2] + step / 12.0 * (23.0 * Fw[2].GetY() - 16.0 * Fw[1].GetY() + 5.0 * Fw[0].GetY());
    Y[12][3] = Y[12][2] + step / 12.0 * (23.0 * Fw[2].GetZ() - 16.0 * Fw[1].GetZ() + 5.0 * Fw[0].GetZ());

    //l0, l1, l2, l3 (dl0/dt = -wx*l1.., etc) Quaternions parameters
    Y[3][3] = Y[3][2] + step / 12.0 * 0.5 *
            (23.0 *(- Y[10][2] * Y[4][2] - Y[11][2] * Y[5][2] - Y[12][2] * Y[6][2]) -
              16.0 * (- Y[10][1] * Y[4][1] - Y[11][1] * Y[5][1] - Y[12][1] * Y[6][1])
            + 5.0 * (- Y[10][0] * Y[4][0] - Y[11][0] * Y[5][0] - Y[12][0] * Y[6][0]));

    Y[4][3] = Y[4][2] + step / 12.0 * 0.5 *
            (23.0 * (Y[10][2] * Y[3][2] + Y[12][2] * Y[5][2] - Y[11][2] * Y[6][2]) -
              16.0 * (Y[10][1] * Y[3][1] + Y[12][1] * Y[5][1] - Y[11][1] * Y[6][1])
            + 5.0 * (Y[10][0] * Y[3][0] + Y[12][0] * Y[5][0] - Y[11][0] * Y[6][0]));

    Y[5][3] = Y[5][2] + step / 12.0 * 0.5 *
            (23.0 * (Y[11][2] * Y[3][2] - Y[12][2] * Y[4][2] + Y[10][2] * Y[6][2]) -
            16.0 * (Y[11][1] * Y[3][1] - Y[12][1] * Y[4][1] + Y[10][1] * Y[6][1])
            + 5.0 * (Y[11][0] * Y[3][0] - Y[12][0] * Y[4][0] + Y[10][0] * Y[6][0]));

    Y[6][3] = Y[6][2] + step / 12.0 * 0.5 *
            (23.0 * (Y[12][2] * Y[3][2] + Y[11][2] * Y[4][2] - Y[10][2] * Y[5][2]) -
            16.0 * (Y[12][1] * Y[3][1] + Y[11][1] * Y[4][1] - Y[10][1] * Y[5][1])
            + 5.0 * (Y[12][0] * Y[3][0] + Y[11][0] * Y[4][0] - Y[10][0] * Y[5][0]));

    //----------------------------------------------------------------------------------------------------------------------------------------//

    //if k = 3 y(i+1)=y(i)+h/24(55f(i)-59f(i-1)+37f(i-2)-9f(i-3)
    for (unsigned i = 4; i < SizeofArrayPoints; i++)
    {
        lambda.lambda0 = Y[3][i - 1];
        lambda.lambda1 = Y[4][i - 1];
        lambda.lambda2 = Y[5][i - 1];
        lambda.lambda3 = Y[6][i - 1];

        log << "[" << tm_str << "]" << " : Iteration " << i << "|Lambda| is " <<
               sqrt(lambda.lambda0 * lambda.lambda0 + lambda.lambda1 * lambda.lambda1 + lambda.lambda2 * lambda.lambda2 + lambda.lambda3 * lambda.lambda3) << std::endl;


        m[0] = 2 * (lambda.lambda0 * lambda.lambda0 + lambda.lambda1 * lambda.lambda1) - 1;
        m[4] = 2 * (lambda.lambda0 * lambda.lambda0 + lambda.lambda2 * lambda.lambda2) - 1;
        m[8] = 2 * (lambda.lambda0 * lambda.lambda0 + lambda.lambda3 * lambda.lambda3) - 1;

        m[1] = 2 * (lambda.lambda1 * lambda.lambda2 - lambda.lambda0 * lambda.lambda3);
        m[2] = 2 * (lambda.lambda1 * lambda.lambda3 + lambda.lambda0 * lambda.lambda2);

        m[3] = 2 * (lambda.lambda1 * lambda.lambda2 + lambda.lambda0 * lambda.lambda3);
        m[5] = 2 * (lambda.lambda2 * lambda.lambda3 - lambda.lambda0 * lambda.lambda1);

        m[6] = 2 * (lambda.lambda1 * lambda.lambda3 - lambda.lambda0 * lambda.lambda2);
        m[7] = 2 * (lambda.lambda2 * lambda.lambda3 + lambda.lambda0 * lambda.lambda1);

        matRotate = m;

        vecForce1[i - 1].SetX( - splineForce->GetValueOfFunction(X[i - 1]) * force_projX);
        vecForce1[i - 1].SetY( - splineForce->GetValueOfFunction(X[i - 1]) * force_projY);
        vecForce1[i - 1].SetZ( - splineForce->GetValueOfFunction(X[i - 1]) * force_projZ);

        //matRotate.Transpose();
        vecForce[i - 1] = Multiply(matRotate, vecForce1[i - 1]);

        for (int j = 0; j < 3; j++)
        {
            //x, y, z (dx/dt = vx, etc)
            Y[j][i] = Y[j][i - 1] + step / 24.0 * (55.0 * Y[7 + j][i - 1] -
                    59.0 * Y[7 + j][i - 2] + 37.0 * Y[7 + j][i - 3] - 9.0 * Y[7 + j][i - 4]);
        }
        log << "[" << tm_str << "]" << " : Iteration " << i << "|Distance| is " <<
               sqrt(Y[0][i] * Y[0][i] + Y[1][i] * Y[1][i] + Y[2][i] * Y[2][i]) << std::endl;

        //Quaternion paramets
        //lambda0
        double fi   = - Y[10][i - 1] * Y[4][i - 1] - Y[11][i - 1] * Y[5][i - 1] - Y[12][i - 1] * Y[6][i - 1];
        double fi_1 = - Y[10][i - 2] * Y[4][i - 2] - Y[11][i - 2] * Y[5][i - 2] - Y[12][i - 2] * Y[6][i - 2];
        double fi_2 = - Y[10][i - 3] * Y[4][i - 3] - Y[11][i - 3] * Y[5][i - 3] - Y[12][i - 3] * Y[6][i - 3];
        double fi_3 = - Y[10][i - 4] * Y[4][i - 4] - Y[11][i - 4] * Y[5][i - 4] - Y[12][i - 4] * Y[6][i - 4];
        Y[3][i] = Y[3][i - 1] + step / 24.0 * 0.5 * (55.0 * fi - 59.0 * fi_1 + 37.0 * fi_2 - 9.0 * fi_3);

        //lambda1
        fi      = Y[10][i - 1] * Y[3][i - 1] + Y[12][i - 1] * Y[5][i - 1] - Y[11][i - 1] * Y[6][i - 1];
        fi_1    = Y[10][i - 2] * Y[3][i - 2] + Y[12][i - 2] * Y[5][i - 2] - Y[11][i - 2] * Y[6][i - 2];
        fi_2    = Y[10][i - 3] * Y[3][i - 3] + Y[12][i - 3] * Y[5][i - 3] - Y[11][i - 3] * Y[6][i - 3];
        fi_3    = Y[10][i - 4] * Y[3][i - 4] + Y[12][i - 4] * Y[5][i - 4] - Y[11][i - 4] * Y[6][i - 4];
        Y[4][i] = Y[4][i - 1] + step / 24.0 * 0.5 * (55.0 * fi - 59.0 * fi_1 + 37.0 * fi_2 - 9.0 * fi_3);

        //lambda2
        fi      = Y[11][i - 1] * Y[3][i - 1] - Y[12][i - 1] * Y[4][i - 1] + Y[10][i - 1] * Y[6][i - 1];
        fi_1    = Y[11][i - 2] * Y[3][i - 2] - Y[12][i - 2] * Y[4][i - 2] + Y[10][i - 2] * Y[6][i - 2];
        fi_2    = Y[11][i - 3] * Y[3][i - 3] - Y[12][i - 3] * Y[4][i - 3] + Y[10][i - 3] * Y[6][i - 3];
        fi_3    = Y[11][i - 4] * Y[3][i - 4] - Y[12][i - 4] * Y[4][i - 4] + Y[10][i - 4] * Y[6][i - 4];
        Y[5][i] = Y[5][i - 1] + step / 24.0 * 0.5 * (55.0 * fi - 59.0 * fi_1 + 37.0 * fi_2 - 9.0 * fi_3);

        //lambda3
        fi      = Y[12][i - 1] * Y[3][i - 1] + Y[11][i - 1] * Y[4][i - 1] - Y[10][i - 1] * Y[5][i - 1];
        fi_1    = Y[12][i - 2] * Y[3][i - 2] + Y[11][i - 2] * Y[4][i - 2] - Y[10][i - 2] * Y[5][i - 2];
        fi_2    = Y[12][i - 3] * Y[3][i - 3] + Y[11][i - 3] * Y[4][i - 3] - Y[10][i - 3] * Y[5][i - 3];
        fi_3    = Y[12][i - 4] * Y[3][i - 4] + Y[11][i - 4] * Y[4][i - 4] - Y[10][i - 4] * Y[5][i - 4];
        Y[6][i] = Y[6][i - 1] + step / 24.0 * 0.5 * (55.0 * fi - 59.0 * fi_1 + 37.0 * fi_2 - 9.0 * fi_3);

        //linear velocities
        Y[7][i] = Y[7][i - 1] + step / 24.0 * (55.0 * vecForce[i - 1].GetX() -
                59.0 * vecForce[i - 2].GetX() + 37.0 * vecForce[i - 3].GetX() -
                9.0 * vecForce[i - 4].GetX()) / mass;
        Y[8][i] = Y[8][i - 1] + step / 24.0 * (55.0 * vecForce[i - 1].GetY() -
                59.0 * vecForce[i - 2].GetY() + 37.0 * vecForce[i - 3].GetY() -
                9.0 * vecForce[i - 4].GetY()) / mass;
        Y[9][i] = Y[9][i - 1] + step / 24.0 * (55.0 * vecForce[i - 1].GetZ() -
                59.0 * vecForce[i - 2].GetZ() + 37.0 * vecForce[i - 3].GetZ() -
                9.0 * vecForce[i - 4].GetZ()) / mass;

        //angular velocities
        //wx, wy, wz (dwx/dt = f(wx, wy, wz, Mx), etc)
        temp_cross = centerofmass.cross(vecForce1[i - 1]);
        Fw[i - 1].SetX((temp_cross.GetX() +
            (tensorOfInertia[4] - tensorOfInertia[8]) * Y[11][i - 1] * Y[12][i - 1]) / tensorOfInertia[0]);
        Fw[i - 1].SetY((temp_cross.GetY() +
            (tensorOfInertia[8] - tensorOfInertia[0]) * Y[10][i - 1] * Y[12][i - 1]) / tensorOfInertia[4]);
        Fw[i - 1].SetZ((temp_cross.GetZ() +
            (tensorOfInertia[0] - tensorOfInertia[4]) * Y[11][i - 1] * Y[10][i - 1]) / tensorOfInertia[8]);
        dOmega[0][i] = Fw[i - 1].GetX();
        dOmega[1][i] = Fw[i - 1].GetY();
        dOmega[2][i] = Fw[i - 1].GetZ();
        Y[10][i] = Y[10][i - 1] + step / 24.0 * (55.0 * Fw[i - 1].GetX() - 59.0 * Fw[i - 2].GetX() + 37.0 * Fw[i - 3].GetX()
                - 9.0 * Fw[i - 4].GetX());
        Y[11][i] = Y[11][i - 1] + step / 24.0 * (55.0 * Fw[i - 1].GetY() - 59.0 * Fw[i - 2].GetY() + 37.0 * Fw[i - 3].GetY()
                - 9.0 * Fw[i - 4].GetY());
        Y[12][i] = Y[12][i - 1] + step / 24.0 * (55.0 * Fw[i - 1].GetZ() - 59.0 * Fw[i - 2].GetZ() + 37.0 * Fw[i - 3].GetZ()
                - 9.0 * Fw[i - 4].GetZ());
    }

    log.flush();
    log.close();

    delete [] Fw;
    Fw = 0;
}

void ODESolveAdams::GetSolveODE(SolveOfEquation *solve) const
{
    ClearSolveStruct(solve);

    solve->numOfEquations = NumOfEquations;
    solve->size = SizeofArrayPoints;
    solve->step = step;
    solve->numOmega = 3;
    solve->X = new double[solve->size];
    solve->Y = new double*[solve->numOfEquations];
    for (int i = 0; i < solve->numOfEquations; i++)
        solve->Y[i] = new double[solve->size];
    for (unsigned i = 0; i < solve->size; i++)
    {
        solve->X[i] = X[i];
        for (int j = 0; j < solve->numOfEquations; j++)
            solve->Y[j][i] = Y[j][i];
    }
    solve->dOmega = new double*[solve->numOmega];
    for (int i = 0; i < solve->numOmega; i++)
        solve->dOmega[i] = new double[solve->size];

    for (unsigned i = 0; i < solve->size; i++)
    {
        for (int j = 0; j < solve->numOmega; j++)
            solve->dOmega[j][i] = dOmega[j][i];
    }
}

void ODESolveAdams::ClearSolveStruct(SolveOfEquation *solve) const
{
    if (solve->size)
    {
        for (int i = 0; i < solve->numOfEquations; i++)
        {
            delete [] (solve->Y[i]);
            solve->Y[i] = 0;
        }
        delete [] (solve->Y);
        solve->Y = 0;

        for (int i = 0; i < solve->numOmega; i++)
        {
            delete [] (solve->dOmega[i]);
            solve->dOmega[i] = 0;
        }
        delete [] (solve->dOmega);
        solve->dOmega = 0;

        delete [] (solve->X);
        solve->X = 0;

        solve->numOfEquations = 0;
        solve->size = 0;
        solve->step = 0;
        solve->numOmega = 0;
    }
}
