#include "render/svg_renderer.hpp"
#include "simple_svg_1.0.0.hpp"
#include "depixelize/color_util.hpp"

#include <vector>
#include <opencv2/opencv.hpp>

namespace Depixelize {

std::string SVGRenderer::get_extension() const
{
    return "svg";
}

void SVGRenderer::render(const std::string &filename, const std::vector<Shape> &shapes) const
{
    svg::Dimensions dimensions(this->input_width * this->output_scale, this->input_height * this->output_scale);
    svg::Document doc(filename, svg::Layout(dimensions, svg::Layout::TopLeft));

    for (auto &s : shapes) {
        std::vector<Point> control_points, curve_points;
        s.edge.make_bezier_segments(&control_points, &curve_points);

        std::vector<svg::Point> svg_control_points, svg_curve_points;
        for (auto &p : control_points) {
            svg_control_points.push_back(svg::Point(this->output_scale * p.x, this->output_scale * p.y));
        }
        for (auto &p : curve_points) {
            svg_curve_points.push_back(svg::Point(this->output_scale * p.x, this->output_scale * p.y));
        }

        double r, g, b;
        r = g = b = 0;
        for (auto &cp : s.colors) {
            cv::Vec3b rgb = yuv2rgb(cp.color);
            r += rgb.val[0] / (255.0 * s.colors.size());
            g += rgb.val[1] / (255.0 * s.colors.size());
            b += rgb.val[2] / (255.0 * s.colors.size());
        }

        doc << svg::QuadraticBezier(svg_control_points, svg_curve_points,
                svg::Fill(svg::Color((int)(255 * r), (int)(255 * g), (int)(255 * b))));
    }

    doc.save();
}

} /* Depixelize */
