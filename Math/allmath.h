/*
 * allmath.h
 *
 *  Created on: 31 авг. 2013 г.
 *      Author: nikita
 */

#ifndef ALLMATH_H_
#define ALLMATH_H_

#include "vector2d.h"
#include "vector3d.h"
#include "matrix2x2.h"
#include "matrix3x3.h"
#include "matrix4x4.h"
#include "cubespline.h"
#include "rkf45.h"

#define TO_RADS(grad, mins, secs) (((grad + mins / 60.0 + secs / 3600.00) * 3.14159265) / 180.0)

#endif /* ALLMATH_H_ */
