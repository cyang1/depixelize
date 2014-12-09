#ifndef DEPIXELIZE_EDGE_HPP
#define DEPIXELIZE_EDGE_HPP

#include <boost/polygon/voronoi.hpp>

#include "geometry/types.hpp"
#include "geometry/point.hpp"

namespace Depixelize {

class Point;

// Undirected edge
class Edge
{
public:

    Edge(Point p1, Point p2) : low(std::min(p1, p2)), high(std::max(p1, p2)) { }
    Edge(const vd_type::edge_type &e) : Edge(Point(*e.vertex0()), Point(*e.vertex1())) { }
    virtual ~Edge() { }
    bool operator==(const Edge &other) const
    {
        return low == other.low && high == other.high;
    }

    Point low, high;
};

} /* Depixelize */

namespace std {

    template <>
    struct hash<Depixelize::Edge>
    {
        std::size_t operator()(const Depixelize::Edge& e) const
        {
            // No bitshifting because the edge is undirected
            using Depixelize::Point;
            return std::hash<Point>()(e.low) ^ std::hash<Point>()(e.high);
        }
    };

} /* std */

#endif /* DEPIXELIZE_EDGE_HPP */
