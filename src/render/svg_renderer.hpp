#ifndef DEPIXELIZE_SVG_RENDERER_HPP
#define DEPIXELIZE_SVG_RENDERER_HPP

#include "render/renderer.hpp"

namespace Depixelize {

class SVGRenderer : public Renderer
{
public:

    using Renderer::Renderer;

    virtual std::string get_extension() const;
    virtual void render(const std::string &filename, const std::vector<Shape> &shapes) const;

};

} /* Depixelize */

#endif /* DEPIXELIZE_SVG_RENDERER_HPP */
