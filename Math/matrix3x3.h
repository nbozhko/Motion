#ifndef MATRIX3X3_H_
#define MATRIX3X3_H_

#include "vector3d.h"
#include <cmath>

class matrix3x3 {
	double mat[9];
	bool transpose;
public:
	matrix3x3();
	matrix3x3(const double *m);
	matrix3x3(const matrix3x3 &m);
    //diagonal matrix
    matrix3x3(double a11, double a22, double a33);

	short size() const {return 9;}

	matrix3x3 multiply(const matrix3x3 &m);
	matrix3x3 multiply(const double *m);
	friend matrix3x3 multiply(const double *m1, const matrix3x3 &m2);

	void MakeIdentity();
	void Transpose();
	const double *GetMatrix() const;

	//operators
	matrix3x3 &operator=(const matrix3x3 &m);
	matrix3x3 &operator=(const double *m);

	matrix3x3 operator+(const matrix3x3 &m);
	matrix3x3 operator+(const double *m);
	void operator+=(const matrix3x3 &m);
	void operator+=(const double *m);
	friend matrix3x3 operator+(const double *m1, const matrix3x3 &m2);

	matrix3x3 operator-(const matrix3x3 &m);
	void operator-=(const matrix3x3 &m);
	matrix3x3 operator-(const double *m);
	void operator-=(const double *m);
	friend matrix3x3 operator-(const double *m1, const matrix3x3 &m2);

	matrix3x3 operator*(const double &op);
	friend matrix3x3 operator*(const double &op, const matrix3x3 &m);
	void operator*=(const double &op);
	matrix3x3 operator/(const double &op);
	void operator/=(const double &op);

	double &operator[](int i);
	double operator[](int i) const;
};

#endif /* MATRIX3X3_H_ */
