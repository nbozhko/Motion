#include "matrix3x3.h"

matrix3x3::matrix3x3() {
	for (int i = 0; i < 9; i++)
		mat[i] = 0;
	transpose = false;
}

matrix3x3::matrix3x3(double a11, double a22, double a33)
{
    for (int i = 0; i < 9; i++)
        mat[i] = 0;
    transpose = false;
    mat[0] = a11;
    mat[4] = a22;
    mat[8] = a33;
}

matrix3x3::matrix3x3(const double *m)
{
	for (int i = 0; i < 9; i++)
		mat[i] = m[i];
	transpose = false;
}

matrix3x3::matrix3x3(const matrix3x3 &m)
{
	for (int i = 0; i < 9; i++)
		mat[i] = m.mat[i];
	transpose = false;
}

matrix3x3 matrix3x3::multiply(const matrix3x3 & m)
{
	matrix3x3 temp;

	temp.mat[0] = mat[0] * m.mat[0] + mat[1] * m.mat[3] + mat[2] * m.mat[6];
	temp.mat[1] = mat[0] * m.mat[1] + mat[1] * m.mat[4] + mat[2] * m.mat[7];
	temp.mat[2] = mat[0] * m.mat[2] + mat[1] * m.mat[5] + mat[2] * m.mat[8];

	temp.mat[3] = mat[3] * m.mat[0] + mat[4] * m.mat[3] + mat[5] * m.mat[6];
	temp.mat[4] = mat[3] * m.mat[1] + mat[4] * m.mat[4] + mat[5] * m.mat[7];
	temp.mat[5] = mat[3] * m.mat[2] + mat[4] * m.mat[5] + mat[5] * m.mat[8];

	temp.mat[6] = mat[6] * m.mat[0] + mat[7] * m.mat[3] + mat[8] * m.mat[6];
	temp.mat[7] = mat[6] * m.mat[1] + mat[7] * m.mat[4] + mat[8] * m.mat[7];
	temp.mat[8] = mat[6] * m.mat[2] + mat[7] * m.mat[5] + mat[8] * m.mat[8];

	return temp;
}

matrix3x3 matrix3x3::multiply(const double *m)
{
	matrix3x3 temp(m);
	return (multiply(temp));
}

matrix3x3 multiply(const double *m1, const matrix3x3 &m2)
{
	matrix3x3 temp;

	temp.mat[0] = m1[0] * m2.mat[0] + m1[1] * m2.mat[3] + m1[2] * m2.mat[6];
	temp.mat[1] = m1[0] * m2.mat[1] + m1[1] * m2.mat[4] + m1[2] * m2.mat[7];
	temp.mat[2] = m1[0] * m2.mat[2] + m1[1] * m2.mat[5] + m1[2] * m2.mat[8];

	temp.mat[3] = m1[3] * m2.mat[0] + m1[4] * m2.mat[3] + m1[5] * m2.mat[6];
	temp.mat[4] = m1[3] * m2.mat[1] + m1[4] * m2.mat[4] + m1[5] * m2.mat[7];
	temp.mat[5] = m1[3] * m2.mat[2] + m1[4] * m2.mat[5] + m1[5] * m2.mat[8];

	temp.mat[6] = m1[6] * m2.mat[0] + m1[7] * m2.mat[3] + m1[8] * m2.mat[6];
	temp.mat[7] = m1[6] * m2.mat[1] + m1[7] * m2.mat[4] + m1[8] * m2.mat[7];
	temp.mat[8] = m1[6] * m2.mat[2] + m1[7] * m2.mat[5] + m1[8] * m2.mat[8];

	return temp;
}

void matrix3x3::MakeIdentity()
{
	for (int i = 0; i < 9; i++)
		mat[i] = 0.0;

	mat[0] = mat[4] = mat[8] = 1.0;
}

void matrix3x3::Transpose()
{
	if (transpose)
		transpose = false;
	else
		transpose = true;

	double temp;
	temp = mat[1];
	mat[1] = mat[3];
	mat[3] = temp;

	temp = mat[2];
	mat[2] = mat[6];
	mat[6] = temp;

	temp = mat[5];
	mat[5] = mat[7];
	mat[7] = temp;
}

const double *matrix3x3::GetMatrix() const
{
	return mat;
}

matrix3x3 &matrix3x3::operator=(const matrix3x3 &m)
{
	for (int j = 0; j < 9; j++)
		this->mat[j] = m.mat[j];

	return (*this);
}

matrix3x3 &matrix3x3::operator=(const double *m)
{
	for (int j = 0; j < 9; j++)
		this->mat[j] = m[j];

	return (*this);
}

matrix3x3 matrix3x3::operator+(const matrix3x3 &m)
{
	matrix3x3 temp;

	for (int j = 0; j < 9; j++)
		temp.mat[j] = mat[j] + m.mat[j];

	return temp;
}

matrix3x3 matrix3x3::operator+(const double *m)
{
	matrix3x3 temp(m);

	return (temp + (*this));
}

void matrix3x3::operator+=(const matrix3x3 &m)
{
	*this = (*this) + m;
}

void matrix3x3::operator+=(const double *m)
{
	matrix3x3 temp(m);

	*this = (*this) + temp;
}

matrix3x3 operator+(const double *m1, const matrix3x3 &m2)
{
	return (matrix3x3(m1) + m2);
}

matrix3x3 matrix3x3::operator-(const matrix3x3 &m)
{
	matrix3x3 temp;

	for (int j = 0; j < 9; j++)
		temp.mat[j] = mat[j] - m.mat[j];

	return temp;
}

matrix3x3 matrix3x3::operator-(const double *m)
{
	matrix3x3 temp(m);

	return (temp - (*this));
}

void matrix3x3::operator-=(const matrix3x3 &m)
{
	*this = (*this) - m;
}

void matrix3x3::operator-=(const double *m)
{
	matrix3x3 temp(m);

	*this = (*this) - temp;
}

matrix3x3 operator-(const double *m1, const matrix3x3 &m2)
{
	return (matrix3x3(m1) - m2);
}

matrix3x3 matrix3x3::operator*(const double &op)
{
	matrix3x3 temp;

	for (int j = 0; j < 9; j++)
		temp.mat[j] = mat[j] * op;

	return temp;
}

matrix3x3 operator*(const double &op, const matrix3x3 &m)
{
	matrix3x3 temp;

	for (int j = 0; j < 9; j++)
		temp.mat[j] = m.mat[j] * op;

	return temp;
}

void matrix3x3::operator*=(const double &op)
{
	*this = (*this) * op;
}

matrix3x3 matrix3x3::operator /(const double & op)
{
	return ((*this) * (1 / op));
}

void matrix3x3::operator /=(const double & op)
{
	*this = (*this) / op;
}

double &matrix3x3::operator[](int i)
{
	if (i >= 9)
		return mat[0];
	return mat[i];
}

double matrix3x3::operator [](int i) const
{
	if (i >= 9)
		return mat[0];
	return mat[i];
}
