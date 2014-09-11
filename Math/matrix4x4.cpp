#include "matrix4x4.h"

matrix4x4::matrix4x4()
{
	transpose = false;
	for (int i = 0; i < 16; i++)
		state[i] = 0;

	state[0] = state[5] = state[10] = state[15] = 1.0;
}

matrix4x4::matrix4x4(const double *m)
{
	transpose = false;
	for (int i = 0; i < 16; i++)
		state[i] = m[i];
}

matrix4x4::matrix4x4(const matrix4x4 &m)
{
	transpose = m.transpose;
	for (int i = 0; i < 16; i++)
		this->state[i] = m.state[i];
}

matrix4x4 matrix4x4::multiply(const double *m)
{
	matrix4 temp;

	temp[0] = state[0] * m[0] + state[1] * m[4] + state[2] * m[8] + state[3] * m[12];
	temp[1] = state[0] * m[1] + state[1] * m[5] + state[2] * m[9] + state[3] * m[13];
	temp[2] = state[0] * m[2] + state[1] * m[6] + state[2] * m[10] + state[3] * m[14];
	temp[3] = state[0] * m[3] + state[1] * m[7] + state[2] * m[11] + state[3] * m[15];

	temp[4] = state[4] * m[0] + state[5] * m[4] + state[6] * m[8] + state[7] * m[12];
	temp[5] = state[4] * m[1] + state[5] * m[5] + state[6] * m[9] + state[7] * m[13];
	temp[6] = state[4] * m[2] + state[5] * m[6] + state[6] * m[10] + state[7] * m[14];
	temp[7] = state[4] * m[3] + state[5] * m[7] + state[6] * m[11] + state[7] * m[15];

	temp[8] = state[8] * m[0] + state[9] * m[4] + state[10] * m[8] + state[11] * m[12];
	temp[9] = state[8] * m[1] + state[9] * m[5] + state[10] * m[9] + state[11] * m[13];
	temp[10] = state[8] * m[2] + state[9] * m[6] + state[10] * m[10] + state[11] * m[14];
	temp[11] = state[8] * m[3] + state[9] * m[7] + state[10] * m[11] + state[11] * m[15];

	temp[12] = state[12] * m[0] + state[13] * m[4] + state[14] * m[8] + state[15] * m[12];
	temp[13] = state[12] * m[1] + state[13] * m[5] + state[14] * m[9] + state[15] * m[13];
	temp[14] = state[12] * m[2] + state[13] * m[6] + state[14] * m[10] + state[15] * m[14];
	temp[15] = state[12] * m[3] + state[13] * m[7] + state[14] * m[11] + state[15] * m[15];

	return matrix4x4(temp);
}

matrix4x4 matrix4x4::multiply(const matrix4x4 &m)
{
	return (this->multiply(m.GetMatrix()));
}

matrix4x4 multiply(const matrix4x4 &m1, const matrix4x4 &m2)
{
	matrix4x4 temp;

	temp = temp.multiply(m1);

	return temp.multiply(m2);
}

matrix4x4 multiply(const double *m1, const matrix4x4 &m2)
{
	matrix4x4 temp(m1);

	return multiply(temp, m2);
}

void matrix4x4::MakeIdentity()
{
	transpose = false;
	for (int i = 0; i < 16; i++)
		state[i] = 0;

	state[0] = state[5] = state[10] = state[15] = 1.0;
}

const double *matrix4x4::GetMatrix() const
{
	return state;
}

void matrix4x4::Transpose()
{
	double temp;

	if (transpose)
		transpose = false;
	else
		transpose = true;

	temp = state[1];
	state[1] = state[4];
	state[4] = temp;

	temp = state[2];
	state[2] = state[8];
	state[8] = temp;

	temp = state[3];
	state[3] = state[12];
	state[12] = temp;

	temp = state[6];
	state[6] = state[9];
	state[9] = temp;

	temp = state[7];
	state[7] = state[10];
	state[10] = temp;

	temp = state[11];
	state[11] = state[14];
	state[14] = temp;
}

void matrix4x4::SetTranslate2D(const vector2d &vec)
{
	state[12] = vec.GetX();
	state[13] = vec.GetY();
}

void matrix4x4::SetTranslate3D(const vector3d &vec)
{
	state[12] = vec.GetX();
	state[13] = vec.GetY();
	state[14] = vec.GetZ();
}

void matrix4x4::MakeRotateMat(const double &angle, const double &x, const double &y, const double &z)
{
	double temp[16];

	for (int j = 0; j < 16; j++)
		temp[j] = 0.0;
	temp[0] = temp[5] = temp[10] = temp[15] = 1.0;

	(*this) = temp;

	if (x != 0.0)
	{
		temp[5] = temp[10] = cos(angle);
		temp[6] = - sin(angle);
		temp[9] = sin(angle);
		(*this) = this->multiply(temp);
	}

	for (int j = 0; j < 16; j++)
		temp[j] = 0.0;
	temp[0] = temp[5] = temp[10] = temp[15] = 1.0;

	if (y != 0.0)
	{
		temp[0] = temp[10] = cos(angle);
		temp[2] = sin(angle);
		temp[8] = - sin(angle);
		(*this) = this->multiply(temp);
	}

	for (int j = 0; j < 16; j++)
		temp[j] = 0.0;
	temp[0] = temp[5] = temp[10] = temp[15] = 1.0;

	if (z != 0.0)
	{
		temp[0] = temp[5] = cos(angle);
		temp[1] = - sin(angle);
		temp[4] = sin(angle);
		(*this) = this->multiply(temp);
	}
}

matrix4x4 &matrix4x4::operator =(const double *m)
{
	if (!m) {
		return (*this);
	}
	for (int i = 0; i < size(); i++)
		this->state[i] = m[i];

	return (*this);
}

matrix4x4 &matrix4x4::operator =(const matrix4x4 &m)
{
	for (int i = 0; i < size(); i++)
		this->state[i] = m.state[i];

	return (*this);
}

matrix4x4 matrix4x4::operator +(const matrix4x4 &m)
{
	matrix4x4 temp;

	for (int i = 0; i < 16; i++)
		temp.state[i] = this->state[i] + m.state[i];

	return temp;
}

matrix4x4 matrix4x4::operator +(const double *m)
{
	matrix4x4 temp;

	for (int i = 0; i < 16; i++)
		temp.state[i] = this->state[i] + m[i];

	return temp;
}

void matrix4x4::operator +=(const double *m)
{
	*this = *this + m;
}

void matrix4x4::operator+=(const matrix4x4 &m)
{
	*this = *this + m;
}

matrix4x4 operator+(const double *m1, const matrix4x4 &m2)
{
	return (matrix4x4(m1) + m2);
}

matrix4x4 matrix4x4::operator -(const matrix4x4 &m)
{
	matrix4x4 temp;

	for (int i = 0; i < 16; i++)
		temp.state[i] = this->state[i] - m.state[i];

	return temp;
}

matrix4x4 matrix4x4::operator -(const double *m)
{
	matrix4x4 temp;

	for (int i = 0; i < 16; i++)
		temp.state[i] = this->state[i] - m[i];

	return temp;
}

void matrix4x4::operator -=(const double *m)
{
	*this = *this - m;
}

void matrix4x4::operator-=(const matrix4x4 &m)
{
	*this = *this - m;
}

matrix4x4 operator-(const double *m1, const matrix4x4 &m2)
{
	return (matrix4x4(m1) - m2);
}

matrix4x4 matrix4x4::operator *(const double &op)
{
	matrix4x4 temp;

	for(int i = 0; i < size(); i++)
		temp[i] = this->state[i] * op;

	return temp;
}

matrix4x4 operator*(const double &op, const matrix4x4 &m)
{
	matrix4x4 temp;

	for(int i = 0; i < m.size(); i++)
		temp[i] = m[i] * op;

	return temp;
}

void matrix4x4::operator *=(const double &op)
{
	*this = matrix4x4((*this) * op);
}

matrix4x4 matrix4x4::operator /(const double &op)
{
	matrix4x4 temp((*this) * (1.0 / op));
	return (temp);
}

void matrix4x4::operator /=(const double &op)
{
	*this = matrix4x4((*this) / op);
}

double &matrix4x4::operator [](int i)
{
	if (i >= 16)
		return state[0];
	return state[i];
}

double matrix4x4::operator [](int i) const
{
	if (i >= 16)
		return state[0];
	return state[i];
}
