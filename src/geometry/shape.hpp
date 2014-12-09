#ifndef DEPIXELIZE_SHAPE_HPP
#define DEPIXELIZE_SHAPE_HPP

#include <cmath>
#include <deque>
#include <vector>

#include <opencv2/opencv.hpp>

#include "geometry/point.hpp"

namespace Depixelize {

class ColorPoint
{
public:

    ColorPoint(Point cen, cv::Vec3b c) : centroid(cen), color(c) { }
    virtual ~ColorPoint() { }

    Point centroid;
    cv::Vec3b color;
};

class Shape
{
public:

    Shape(const std::deque<Point> &edge, const std::vector<ColorPoint> &colors)
        : edge(edge), colors(colors), area(0)
    {
        // All credits to http://alienryderflex.com/polygon_area/
        uint32_t j = edge.size() - 1;
        for (uint32_t i = 0; i < edge.size(); i++) {
            this->area += (edge[j].x + edge[i].x) * (edge[j].y - edge[i].y);
            j = i;
        }

        this->area = fabs(this->area * 0.5);
    }
    virtual ~Shape() { }

    // Sorting functions
    bool operator<(const Shape &other) const
    {
        return this->area < other.area;
    }

    friend void swap(Shape &a, Shape &b)
    {
        using std::swap;

        swap(a.edge, b.edge);
        swap(a.colors, b.colors);
        swap(a.area, b.area);
    }

    std::deque<Point> edge;
    std::vector<ColorPoint> colors;
    double area;

};

} /* Depixelize */

#endif /* DEPIXELIZE_SHAPE_HPP */
