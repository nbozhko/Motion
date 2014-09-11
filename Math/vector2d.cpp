#include "vector2d.h"

vector2d::vector2d()
{
	x = y = length = angle = 0;
}

vector2d::vector2d(cref_d _x, cref_d _y)
{
	x = _x;
	y = _y;
	angle = atan2(y, x);
	length = sqrt(x*x + y*y);
}

vector2d::vector2d(const vector2d &v)
{
	x = v.x;
	y = v.y;
	angle = v.angle;
	length = v.length;
}

double vector2d::GetX() const
{
	return x;
}

double vector2d::GetY() const
{
	return y;
}

double vector2d::GetAngle() const
{
	return angle;
}

double vector2d::GetLength() const
{
	return length;
}

void vector2d::setNewPair(cref_d _x, cref_d _y)
{
	x = _x;
	y = _y;
	angle = atan2(y, x);
	length = sqrt(x*x + y*y);
}

vector2d vector2d::GetNorm() const
{
	return vector2d(x / length, y / length);
}

double vector2d::dot(const vector2d &v) const
{
	return (this->x * v.x + this->y * v.y);
}

double dot(const vector2d &v1, const vector2d &v2)
{
	return (v1.x * v2.x + v1.y * v2.y);
}

double AngleV1V2(const vector2d &v1, const vector2d &v2)
{
	return ( (v1.x * v2.x + v1.y * v2.y) / (v1.GetLength() * v2.GetLength()) );
}

vector2d &vector2d::operator =(const vector2d &v)
{
	this->x = v.x;
	this->y = v.y;
	this->angle = atan2(y, x);
	this->length = sqrt(x*x + y*y);

	return (*this);
}

vector2d vector2d::operator +(const vector2d &v)
{
	return vector2d(this->x + v.x, this->y + v.y);
}

vector2d operator+(const vector2d &v1, const vector2d &v2)
{
	return vector2d(v1.x + v2.x, v1.y + v2.y);
}

void vector2d::operator +=(const vector2d &v)
{
	*this = vector2d(this->x + v.x, this->y + v.y);
}

vector2d vector2d::operator -(const vector2d &v)
{
	return vector2d(this->x - v.x, this->y - v.y);
}

vector2d operator-(const vector2d &v1, const vector2d &v2)
{
	return vector2d(v1.x - v2.x, v1.y - v2.y);
}

void vector2d::operator -=(const vector2d &v)
{
	*this = vector2d(this->x - v.x, this->y - v.y);
}

vector2d vector2d::operator *(cref_d op)
{
	return vector2d(this->x * op, this->y * op);
}

vector2d operator*(cref_d op, const vector2d &v)
{
	return vector2d(v.x * op, v.y * op);
}

void vector2d::operator *=(cref_d op)
{
	*this = vector2d(this->x * op, this->y * op);
}

vector2d vector2d::operator/(cref_d num)
{
	return vector2d(this->x / num, this->y / num);
}

void vector2d::operator/=(cref_d num)
{
	*this = vector2d(this->x / num, this->y / num);
}
