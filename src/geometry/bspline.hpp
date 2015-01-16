#ifndef DEPIXELIZE_BSPLINE_HPP
#define DEPIXELIZE_BSPLINE_HPP

#include "geometry/point.hpp"

#include <cmath>
#include <functional>
#include <cstdint>
#include <vector>

namespace Depixelize {

class Point;

// A B-spline has:
//   * n + 1 control points
//   * m + 1 knots
//   * degree p
//   * where m = n + p + 1
class BSpline
{
public:

    BSpline() { }

    // Complete B-spline constructor
    BSpline(const std::vector<double> &knots, const std::vector<Point*> points, uint32_t degree = 2)
        : knots(knots), points(points), degree(degree) { }

    // Closed B-spline from polyline
    BSpline(const std::vector<Point*> &points, uint32_t degree = 2)
        : points(points), degree(degree)
    {
        for (uint32_t i = 0; i < degree; i++) {
            this->points.push_back(this->points[i]);
        }
        uint32_t m = this->points.size() + degree;
        for (uint32_t i = 0; i < m + 1; i++) {
            this->knots.push_back(i * 1.0 / m);
        }
    }

    virtual ~BSpline() { }

    // De Boor's Algorithm. Given a knot u, returns a point on the curve C(u)
    Point operator()(double u) const
    {
        uint32_t i;
        for (i = 0; i < this->knots.size() && knots[i] <= u; i++);
        return de_boor(this->degree, i - 1, u);
    }

    // See Wikipedia for details
    Point de_boor(uint32_t k, uint32_t i, double u) const
    {
        if (k == 0) {
            return *this->points[i % this->points.size()];
        } else {
            int next_knot = (i + this->degree + 1 - k) % this->knots.size();
            if (next_knot < 0) next_knot += this->knots.size();
            double alpha = (u - this->knots[i]) / (this->knots[next_knot] - this->knots[i]);
            uint32_t dec_i = (i == 0 ? this->knots.size() : i) - 1;
            return this->de_boor(k - 1, dec_i, u) * (1 - alpha) + this->de_boor(k - 1, i, u) * alpha;
        }
    }

    // TODO: Currently only works for quadratic curves
    void make_bezier_segments(std::vector<Point> *control_points, std::vector<Point> *curve_points) const
    {
        assert(this->degree == 2);

        control_points->clear();
        curve_points->clear();

        for (uint32_t i = 1; i < this->points.size() - 1; i++) {
            control_points->push_back(*this->points[i]);
        }
        for (uint32_t i = 2; i < this->knots.size() - 2; i++) {
            curve_points->push_back((*this)(this->knots[i]));
        }
    }

    Point domain() const
    {
        return Point(this->knots[this->degree], this->knots[this->knots.size() - this->degree - 1]);
    }

    Point get_span(uint32_t index) const
    {
        double d0 = this->knots[index];
        double d1 = this->knots[index + 1];
        Point dom = this->domain();

        if (d0 < dom.x) {
            d0 = this->knots[index + this->points.size() - this->degree];
            d1 = this->knots[index + this->points.size() - this->degree + 1];
        } else if (d1 > dom.y) {
            d0 = this->knots[index + this->degree - this->points.size()];
            d1 = this->knots[index + this->degree - this->points.size() + 1];
        }
        // Clamp
        return Point(fmax(dom.x, fmin(d0, dom.y)), fmax(dom.x, fmin(d1, dom.y)));
    }

    BSpline derivative() const
    {
        std::vector<Point*> new_points;
        for (uint32_t i = 0; i < this->points.size() - 1; i++) {
            double coeff = this->degree / (this->knots[i + 1 + this->degree] - this->knots[i + 1]);
            new_points.push_back(new Point((*this->points[i + 1] - *this->points[i]) * coeff));
        }
        std::vector<double> new_knots(this->knots.begin() + 1, this->knots.end() - 1);

        return BSpline(new_knots, new_points, this->degree - 1);
    }

    double curvature(double u) const
    {
        BSpline deriv1 = this->derivative();
        BSpline deriv2 = deriv1.derivative();

        Point d1 = deriv1(u);
        Point d2 = deriv2(u);

        double num = d1.x * d2.y - d1.y * d2.x;
        double den = pow(d1.x * d1.x + d1.y * d1.y, 3 / 2.0);

        for (auto &p : deriv1.points) {
            delete p;
        }
        for (auto &p : deriv2.points) {
            delete p;
        }

        double q = num / den;
        return (std::isinf(q) || std::isnan(q))
            ? 0
            : fabs(q);
    }

    double integrate_span(std::function<double(double)> f, const Point &span) const
    {
        // Numerical integral
        const uint32_t NUM_SAMPLES = 16;
        double ret = 0;
        double width = (span.y - span.x) / NUM_SAMPLES;

        for (uint32_t i = 0; i < NUM_SAMPLES; i++) {
            ret += width * f(span.x + width * (i + 0.5));
        }

        return ret;
    }

    double curvature_energy(uint32_t index) const
    {
        auto curvature_fp = std::bind(&BSpline::curvature, this, std::placeholders::_1);
        double ret = 0;
        for (uint32_t i = 0; i < this->degree; i++) {
            ret += this->integrate_span(curvature_fp, this->get_span(index + i));
        }
        return ret;
    }

private:
    std::vector<double> knots;
    std::vector<Point*> points;
    uint32_t degree;
};

} /* Depixelize */

#endif /* DEPIXELIZE_BSPLINE_HPP */
