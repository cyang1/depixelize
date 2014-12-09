#ifndef DEPIXELIZE_RENDERER_HPP
#define DEPIXELIZE_RENDERER_HPP

#include <vector>

#include "geometry/shape.hpp"

namespace Depixelize {

class Renderer
{
public:
    Renderer(int w, int h, double scale) : input_width(w), input_height(h), output_scale(scale) { }
    virtual ~Renderer() { }

    std::string make_full_filename(std::string stem) const
    {
        return stem + "." + this->get_extension();
    }

    virtual std::string get_extension() const = 0;
    virtual void render(const std::string &filename, const std::vector<Shape> &shapes) const = 0;

private:

    int input_width, input_height;
    double output_scale;

};

} /* Depixelize */

#endif /* DEPIXELIZE_RENDERER_HPP */
