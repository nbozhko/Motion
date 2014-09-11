#ifndef VECTOR_H_
#define VECTOR_H_

#include <cmath>

typedef const double & cref_d;

class vector2d {
	double x, y;
	double length;
	double angle;
public:
	vector2d();
	vector2d(cref_d _x, cref_d y);
	vector2d(const vector2d &v);

	//set operations
	void setNewPair(cref_d _x, cref_d _y);

	//gets operations
	double GetX() const;
	double GetY() const;
	double GetLength() const;
	double GetAngle() const;
	vector2d GetNorm() const;

	//vector operations
	double dot(const vector2d &v) const;
	friend double dot(const vector2d &v1, const vector2d &v2);
	friend double AngleV1V2(const vector2d &v1, const vector2d &v2);
	//vector2d cross(const vector2d &v) const;
	//friend vector2d cross(const vector2d &v);
	vector2d operator+(const vector2d &v);
	friend vector2d operator+(const vector2d &v1, const vector2d &v2);
	vector2d operator-(const vector2d &v);
	friend vector2d operator-(const vector2d &v1, const vector2d &v2);
	void operator+=(const vector2d &v);
	void operator-=(const vector2d &v);
	vector2d operator*(cref_d num);
	friend vector2d operator*(cref_d op, const vector2d &v);
	void operator*=(cref_d num);
	vector2d operator/(cref_d num);
	void operator/=(cref_d num);
	vector2d &operator=(const vector2d &v);
};


#endif /* VECTOR_H_ */
