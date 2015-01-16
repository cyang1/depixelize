#ifndef DEPIXELIZE_POINT_HPP
#define DEPIXELIZE_POINT_HPP

#include <cmath>
#include <limits>

#include <boost/polygon/voronoi.hpp>

#include "geometry/types.hpp"

namespace Depixelize {

class Point
{
public:

    Point(double x, double y) : x(x), y(y) { }
    Point(const vd_type::vertex_type &v) : x(v.x()), y(v.y()) { }
    virtual ~Point() { }
    bool operator==(const Point &other) const
    {
        return
          fabs(this->x - other.x) < std::numeric_limits<double>::epsilon() &&
          fabs(this->y - other.y) < std::numeric_limits<double>::epsilon();
    }

    bool operator<(const Point &other) const
    {
        return
          this->x < other.x || (
            fabs(this->x - other.x) < std::numeric_limits<double>::epsilon() &&
            this->y < other.y
          );
    }

    Point operator+(const Point &other) const
    {
        return Point(this->x + other.x, this->y + other.y);
    }

    Point operator-(const Point &other) const
    {
        return Point(this->x - other.x, this->y - other.y);
    }

    Point operator*(double m) const
    {
        return Point(this->x * m, this->y * m);
    }

    Point operator/(double m) const
    {
        return Point(this->x / m, this->y / m);
    }

    double x, y;
    // Spline *owner;
};

} /* Depixelize */

namespace std {

    template <>
    struct hash<Depixelize::Point>
    {
        std::size_t operator()(const Depixelize::Point& p) const
        {
            return std::hash<double>()(p.x) ^ (std::hash<double>()(p.y) << 1);
        }
    };

}

#endif /* DEPIXELIZE_POINT_HPP */
