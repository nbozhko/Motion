#ifndef MATRIX_H_
#define MATRIX_H_

#include <cmath>
#include "vector3d.h"
#include "matrix3x3.h"
#include "matrix2x2.h"

typedef double matrix4[16];

class matrix4x4 {
	matrix4 state;
	bool transpose;
public:
	matrix4x4();
	matrix4x4(const double *m);
	matrix4x4(const matrix4x4 &m);

	short size() const {return 16;}
	//matrix operations
	matrix4x4 multiply(const matrix4x4 &m);
	friend matrix4x4 multiply(const matrix4x4 &m1, const matrix4x4 &m2);
	matrix4x4 multiply(const double *m);
	friend matrix4x4 multiply(const double *m1, const matrix4x4 &m2);

	void MakeIdentity();
	const double *GetMatrix() const;
	void Transpose();

	//set positions
	void SetTranslate2D(const vector2d &vec);
	void SetTranslate3D(const vector3d &vec);

	void MakeRotateMat(const double &angle,const double &x, const double &y, const double &z);//Not Release yet

	//operators
	matrix4x4 &operator=(const matrix4x4 &m);
	matrix4x4 &operator=(const double *m);

	matrix4x4 operator+(const matrix4x4 &m);
	matrix4x4 operator+(const double *m);
	void operator+=(const matrix4x4 &m);
	void operator+=(const double *m);
	friend matrix4x4 operator+(const double *m1, const matrix4x4 &m2);

	matrix4x4 operator-(const matrix4x4 &m);
	void operator-=(const matrix4x4 &m);
	matrix4x4 operator-(const double *m);
	void operator-=(const double *m);
	friend matrix4x4 operator-(const double *m1, const matrix4x4 &m2);

	matrix4x4 operator*(const double &op);
	friend matrix4x4 operator*(const double &op, const matrix4x4 &m);
	void operator*=(const double &op);
	matrix4x4 operator/(const double &op);
	void operator/=(const double &op);

	double &operator[](int i);
	double operator[](int i) const;
};

#endif /* MATRIX_H_ */
