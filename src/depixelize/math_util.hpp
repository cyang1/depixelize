#ifndef DEPIXELIZE_MATH_UTIL_HPP
#define DEPIXELIZE_MATH_UTIL_HPP

#include <cmath>

#include "geometry/point.hpp"

namespace Depixelize {

// acos(dot(v1, v2) / (||v1|| * ||v2||))
inline double vector_angle(Point v1, Point v2)
{
    double v1_mag = sqrt(v1.x * v1.x + v1.y * v1.y);
    double v2_mag = sqrt(v2.x * v2.x + v2.y * v2.y);
    return acos((v1.x * v2.x + v1.y * v2.y) / v1_mag * v2_mag);
}

} /* Depixelize */

#endif /* DEPIXELIZE_MATH_UTIL_HPP */
