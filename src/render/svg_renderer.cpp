#include "render/svg_renderer.hpp"

namespace Depixelize {

std::string SVGRenderer::get_extension() const
{
    return "svg";
}

void SVGRenderer::render(const std::string &filename, const std::vector<Shape> &shapes) const
{
}

} /* Depixelize */
