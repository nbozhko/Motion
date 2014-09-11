#include "vector3d.h"

vector3d::vector3d()
{
	x = y = z = 0;
	length = 0;
}

vector3d::vector3d(const double & _x, const double & _y, const double & _z)
{
	x = _x;
	y = _y;
	z = _z;
	length = this->GetLength();
}

vector3d::vector3d(const vector3d &v)
{
	x = v.x;
	y = v.y;
	z = v.z;
	length = v.length;
}

vector3d::vector3d(const vector2d & v)
{
	x = v.GetX();
	y = v.GetY();
	z = 0;
	length = v.GetLength();
}

void vector3d::SetX(const double & _x)
{
	x = _x;
}

void vector3d::SetY(const double & _y)
{
	y = _y;
}

void vector3d::SetZ(const double &_z)
{
	z = _z;
}

void vector3d::SetNewCoords(const double &_x, const double &_y, const double &_z)
{
	x = _x;
	y = _y;
	z = _z;
}

double vector3d::GetX() const
{
	return x;
}

double vector3d::GetY() const
{
	return y;
}

double vector3d::GetZ() const
{
	return z;
}

double vector3d::GetLength() const
{
	return (sqrt(x * x + y * y + z * z));
}

vector3d vector3d::GetNormVector() const
{
	vector3d temp;

	temp.x = x / GetLength();
	temp.y = y / GetLength();
	temp.z = z / GetLength();

	return temp;
}

double vector3d::dot(const vector3d & v)
{
	return (x * v.x + y * v.y + z * v.z);
}

vector3d vector3d::cross(const vector3d & v)
{
	return vector3d(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - v.x * y);
}

vector3d vector3d::cross(const vector2d & v)
{
	return vector3d(0, 0, x * v.GetY() - v.GetX() * y);
}

vector3d vector3d::cross(const vector2d &v1, const vector2d & v2)
{
	return vector3d(0, 0, v1.GetX() * v2.GetY() - v2.GetX() * v1.GetY());
}

vector3d vector3d::operator+(const vector3d &v)
{
	return vector3d(x + v.x, y + v.y, z + v.z);
}

vector3d vector3d::operator +(const vector2d & v)
{
	return vector3d(x + v.GetX(), y + v.GetY(), z);
}

void vector3d::operator +=(const vector3d & v)
{
	this->x += v.x;
	this->y += v.y;
	this->z += v.z;
}

void vector3d::operator +=(const vector2d & v)
{
	this->x += v.GetX();
	this->y += v.GetY();
}

vector3d vector3d::operator-(const vector3d &v)
{
	return vector3d(x - v.x, y - v.y, z - v.z);
}

vector3d vector3d::operator -(const vector2d & v)
{
	return vector3d(x - v.GetX(), y - v.GetY(), z);
}

void vector3d::operator -=(const vector3d & v)
{
	this->x -= v.x;
	this->y -= v.y;
	this->z -= v.z;
}

void vector3d::operator -=(const vector2d & v)
{
	this->x -= v.GetX();
	this->y -= v.GetY();
}

vector3d vector3d::operator*(const double & num)
{
	vector3d temp;

	temp.x = x * num;
	temp.y = y * num;
	temp.z = z * num;

	return temp;
}

vector3d operator*(const double & op, const vector3d &v)
{
	vector3d temp;

	temp.x = v.x * op;
	temp.y = v.y * op;
	temp.z = v.z * op;

	return temp;
}

void vector3d::operator*=(const double & num)
{
	this->x *= num;
	this->y *= num;
	this->z *= num;
}

vector3d vector3d::operator/(const double & num)
{
	vector3d temp;

	temp.x = x / num;
	temp.y = y / num;
	temp.z = z / num;

	return temp;
}

void vector3d::operator/=(const double & num)
{
	if (num != 0)
	{
		this->x /= num;
		this->y /= num;
		this->z /= num;
	}
}

vector3d &vector3d::operator=(const vector3d &v)
{
	this->x = v.x;
	this->y = v.y;
	this->z = v.z;
	this->length = v.GetLength();

	return (*this);
}

vector3d &vector3d::operator=(const vector2d &v)
{
	this->x = v.GetX();
	this->y = v.GetY();
	this->z = 0;
	this->length = v.GetLength();

	return (*this);
}
