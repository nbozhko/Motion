#include "matrix2x2.h"

matrix2x2::matrix2x2()
{
	for (int i = 0; i < 4; i++)
		mat[i] = 0;
	transpose = false;
}

matrix2x2::matrix2x2(const matrix2x2 &m)
{
	for (int i = 0; i < 4; i++)
		mat[i] = m.mat[i];
	transpose = false;
}

matrix2x2::matrix2x2(const double *m)
{
	for (int i = 0; i < 4; i++)
		mat[i] = m[i];

	transpose = false;
}

matrix2x2 matrix2x2::multiply(const matrix2x2 &m)
{
	matrix2x2 temp;

	temp.mat[0] = mat[0] * m.mat[0] + mat[1] * m.mat[2];
	temp.mat[1] = mat[0] * m.mat[1] + mat[1] * m.mat[3];

	temp.mat[2] = mat[2] * m.mat[0] + mat[3] * m.mat[2];
	temp.mat[3] = mat[2] * m.mat[1] + mat[3] * m.mat[3];

	return temp;
}

matrix2x2 multiply(const matrix2x2 &m1, const matrix2x2 &m2)
{
	matrix2x2 temp;

	temp.mat[0] = m1.mat[0] * m2.mat[0] + m1.mat[1] * m2.mat[2];
	temp.mat[1] = m1.mat[0] * m2.mat[1] + m1.mat[1] * m2.mat[3];

	temp.mat[2] = m1.mat[2] * m2.mat[0] + m1.mat[3] * m2.mat[2];
	temp.mat[3] = m1.mat[2] * m2.mat[1] + m1.mat[3] * m2.mat[3];

	return temp;
}

matrix2x2 matrix2x2::multiply(const double *m)
{
	matrix2x2 temp;

	temp.mat[0] = mat[0] * m[0] + mat[1] * m[2];
	temp.mat[1] = mat[0] * m[1] + mat[1] * m[3];

	temp.mat[2] = mat[2] * m[0] + mat[3] * m[2];
	temp.mat[3] = mat[2] * m[1] + mat[3] * m[3];

	return temp;
}

matrix2x2 multiply(const double *m1, const matrix2x2 &m2)
{
	matrix2x2 temp;

	temp.mat[0] = m1[0] * m2.mat[0] + m1[1] * m2.mat[2];
	temp.mat[1] = m1[0] * m2.mat[1] + m1[1] * m2.mat[3];

	temp.mat[2] = m1[2] * m2.mat[0] + m1[3] * m2.mat[2];
	temp.mat[3] = m1[2] * m2.mat[1] + m1[3] * m2.mat[3];

	return temp;
}

void matrix2x2::MakeIdentity()
{
	for (int i = 0; i < 4; i++)
		mat[i] = 0;

	mat[0] = mat[3] = 1.0;
}

const double *matrix2x2::GetMatrix() const
{
	return mat;
}

void matrix2x2::Transpose()
{
	transpose = true;

	double temp = mat[1];

	mat[1] = mat[2];
	mat[2] = temp;
}

void matrix2x2::SetRot2D(const double &angle)
{
	matrix2x2 temp;

	temp.mat[0] = temp.mat[3] = cos(angle);

	temp.mat[1] = - sin(angle);
	temp.mat[2] = sin(angle);

	*this = ::multiply(temp, *this);
}

matrix2x2 &matrix2x2::operator=(const matrix2x2 &m)
{
	for(int i = 0; i < 4; i++)
	{
		mat[i] = m.mat[i];
	}

	return (*this);
}

matrix2x2 &matrix2x2::operator=(const double *m)
{
	for(int i = 0; i < 4; i++)
	{
		mat[i] = m[i];
	}

	return (*this);
}

matrix2x2 matrix2x2::operator +(const double *m)
{
	matrix2x2 temp;

	for (int i = 0; i < 4; i++)
		temp.mat[i] = this->mat[i] + m[i];

	return temp;
}

matrix2x2 matrix2x2::operator+(const matrix2x2 &m)
{
	matrix2x2 temp;

	for (int i = 0; i < 4; i++)
		temp.mat[i] = m.mat[i] + m[i];

	return temp;
}

void matrix2x2::operator +=(const matrix2x2 &m)
{
	for (int i = 0; i < 4; i++)
		this->mat[i] += m.mat[i];
}

void matrix2x2::operator+=(const double *m)
{
	for (int i = 0; i < 4; i++)
		this->mat[i] += m[i];
}

matrix2x2 operator+(const double *m1, const matrix2x2 &m2)
{
	return (matrix2x2(m1) + m2);
}

matrix2x2 matrix2x2::operator -(const double *m)
{
	matrix2x2 temp;

	for (int i = 0; i < 4; i++)
		temp.mat[i] = this->mat[i] - m[i];

	return temp;
}

matrix2x2 matrix2x2::operator-(const matrix2x2 &m)
{
	matrix2x2 temp;

	for (int i = 0; i < 4; i++)
		temp.mat[i] = m.mat[i] - m[i];

	return temp;
}

void matrix2x2::operator -=(const matrix2x2 &m)
{
	for (int i = 0; i < 4; i++)
		this->mat[i] -= m.mat[i];
}

void matrix2x2::operator-=(const double *m)
{
	for (int i = 0; i < 4; i++)
		this->mat[i] -= m[i];
}

matrix2x2 operator-(const double *m1, const matrix2x2 &m2)
{
	return (matrix2x2(m1) - m2);
}

matrix2x2 matrix2x2::operator*(const double &op)
{
	matrix2x2 temp;

	for (int i = 0; i < 4; i++)
		temp.mat[i] = this->mat[i] * op;

	return temp;
}

matrix2x2 operator*(const double &op, const matrix2x2 &m)
{
	matrix2x2 temp;

	for (int i = 0; i < 4; i++)
		temp.mat[i] = m.mat[i] * op;

	return temp;
}

void matrix2x2::operator*=(const double &op)
{
	for (int i = 0; i < 4; i++)
		this->mat[i] *= op;
}

matrix2x2 matrix2x2::operator/(const double &op)
{
	matrix2x2 temp = *this;

	return (temp * (1.0 / op));
}

void matrix2x2::operator/=(const double &op)
{
	for (int i = 0; i < 4; i++)
		this->mat[i] /= op;
}

double &matrix2x2::operator[](int i)
{
	if (i >= 4)
		return mat[0];
	return mat[i];
}

double matrix2x2::operator[](int i) const
{
	if (i >= 4)
		return mat[0];
	return mat[i];
}
