#ifndef VECTOR3D_H_
#define VECTOR3D_H_

#include <cmath>
#include "vector2d.h"

class vector3d {
	double x, y, z;
	double length;
public:
	vector3d();
	vector3d(const double &_x, const double &_y, const double &_z);
	vector3d(const vector3d &v);
	vector3d(const vector2d &v);

	void SetNewCoords(const double &_x, const double &_y, const double &_z);
	void SetX(const double &_x);
	void SetY(const double &_y);
	void SetZ(const double &_z);

	double GetX() const;
	double GetY() const;
	double GetZ() const;
	double GetLength() const;
	vector3d GetNormVector() const;

	double dot(const vector3d &v);
	vector3d cross(const vector3d &v);
	vector3d cross(const vector2d &v);
    vector3d cross(const vector2d &v1, const vector2d &v2);
	vector3d operator+(const vector3d &v);
	vector3d operator+(const vector2d &v);
	friend vector3d operator+(const vector2d &v1, const vector3d &v2);
	vector3d operator-(const vector3d &v);
	vector3d operator-(const vector2d &v);
	friend vector3d operator-(const vector2d &v1, const vector3d &v2);
	void operator+=(const vector3d &v);
	void operator-=(const vector3d &v);
	void operator+=(const vector2d &v);
	void operator-=(const vector2d &v);
	vector3d operator*(const double & num);
	friend vector3d operator*(const double & op, const vector3d &v);
	void operator*=(const double & num);
	vector3d operator/(const double & num);
	void operator/=(const double & num);
	vector3d &operator=(const vector3d &v);
	vector3d &operator=(const vector2d &v);
};


#endif /* VECTOR3D_H_ */
