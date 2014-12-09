/*
 * Miscellaneous useful types
 */

#ifndef DEPIXELIZE_TYPES_HPP
#define DEPIXELIZE_TYPES_HPP

#include <deque>
#include <boost/polygon/voronoi.hpp>

namespace Depixelize {

class Point;

typedef boost::polygon::voronoi_diagram<double> vd_type;

} /* Depixelize */

#endif /* DEPIXELIZE_TYPES_HPP */
