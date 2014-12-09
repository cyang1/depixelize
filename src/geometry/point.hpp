#ifndef DEPIXELIZE_POINT_HPP
#define DEPIXELIZE_POINT_HPP

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
        return this->x == other.x && this->y == other.y;
    }

    bool operator<(const Point &other) const
    {
        return this->x < other.x || (this->x == other.x && this->y < other.y);
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
