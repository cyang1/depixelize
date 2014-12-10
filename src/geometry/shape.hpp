#ifndef DEPIXELIZE_SHAPE_HPP
#define DEPIXELIZE_SHAPE_HPP

#include <cmath>
#include <deque>
#include <vector>

#include <opencv2/opencv.hpp>

#include "geometry/point.hpp"
#include "geometry/bspline.hpp"

namespace Depixelize {

class ColorPoint
{
public:

    ColorPoint(const Point &cen, const cv::Vec3b &c) : centroid(cen), color(c) { }
    virtual ~ColorPoint() { }

    Point centroid;
    cv::Vec3b color;
};

class Shape
{
public:

    Shape(const BSpline &edge, const std::vector<ColorPoint> &colors, double area)
        : edge(edge), colors(colors), area(area) { }
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

    BSpline edge;
    std::vector<ColorPoint> colors;
    double area;

};

} /* Depixelize */

#endif /* DEPIXELIZE_SHAPE_HPP */
