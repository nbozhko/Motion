#ifndef MATRIX2X2_H_
#define MATRIX2X2_H_

#include <cmath>

class matrix2x2 {
	double mat[4];
	bool transpose;
public:
	matrix2x2();
	matrix2x2(const double *m);
	matrix2x2(const matrix2x2 &m);

	short size() const {return 4;}

	//matrix operations
	matrix2x2 multiply(const matrix2x2 &m);
	friend matrix2x2 multiply(const matrix2x2 &m1, const matrix2x2 &m2);
	matrix2x2 multiply(const double *m);
	friend matrix2x2 multiply(const double *m1, const matrix2x2 &m2);

	void MakeIdentity();
	const double *GetMatrix() const;
	void Transpose();

	void SetRot2D(const double &angle);

	//operators
	matrix2x2 &operator=(const matrix2x2 &m);
	matrix2x2 &operator=(const double *m);

	matrix2x2 operator+(const matrix2x2 &m);
	matrix2x2 operator+(const double *m);
	void operator+=(const matrix2x2 &m);
	void operator+=(const double *m);
	friend matrix2x2 operator+(const double *m1, const matrix2x2 &m2);

	matrix2x2 operator-(const matrix2x2 &m);
	void operator-=(const matrix2x2 &m);
	matrix2x2 operator-(const double *m);
	void operator-=(const double *m);
	friend matrix2x2 operator-(const double *m1, const matrix2x2 &m2);//not release

	matrix2x2 operator*(const double &op);
	friend matrix2x2 operator*(const double &op, const matrix2x2 &m);
	void operator*=(const double &op);
	matrix2x2 operator/(const double &op);
	void operator/=(const double &op);

	double &operator[](int i);
	double operator[](int i) const;

};

#endif /* MATRIX2X2_H_ */
