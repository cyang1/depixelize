#ifndef DEPIXELIZE_SPLINE_OPTIMIZER_HPP
#define DEPIXELIZE_SPLINE_OPTIMIZER_HPP

#include <cstdint>
#include <vector>
#include <deque>
#include <unordered_map>

#include <opencv2/opencv.hpp>
#include <boost/polygon/voronoi.hpp>
#include <boost/polygon/point_data.hpp>
#include <boost/polygon/segment_data.hpp>

#include "geometry/types.hpp"
#include "geometry/point.hpp"
#include "geometry/edge.hpp"
#include "geometry/shape.hpp"

namespace Depixelize {

struct PointRef {
    uint32_t idx;
    bool can_optimize;
};

struct EdgeRef {
    uint32_t idx1, idx2;
    bool is_shading_edge;

    bool operator==(const EdgeRef &other) const
    {
        return (this->idx1 == other.idx1 && this->idx2 == other.idx2) ||
            (this->idx1 == other.idx2 && this->idx2 == other.idx1);
    }

    bool operator!=(const EdgeRef &other) const
    {
        return !(*this == other);
    }
};

typedef std::deque<PointRef> Path;

class SplineOptimizer
{
public:

    SplineOptimizer(const vd_type &vd,
                    const std::vector< boost::polygon::point_data<int> > &points,
                    const std::vector< boost::polygon::segment_data<int> > &segments,
                    const std::vector<uint32_t> &components,
                    uint32_t num_components,
                    const std::vector<cv::Vec3b> &colors)
        : vd(vd), point_data(points), segment_data(segments), components(components), num_components(num_components), colors(colors),
          component_paths(NULL), all_points(), point_map() { }

    virtual ~SplineOptimizer() {
        if (component_paths != NULL) {
            delete[] component_paths;
        }
    }

    void initialize();
    void optimize_splines();
    std::vector<Shape> make_shapes();

private:

    void join_edges(Path *path,
        const std::vector<EdgeRef> &edges,
        const std::unordered_multimap<uint32_t, EdgeRef> &adjacent_edges);

    bool should_optimize(EdgeRef e1, EdgeRef e2,
        const std::unordered_multimap<uint32_t, EdgeRef> &adjacent_edges);

    // These next 3 functions were borrowed from the voronoi_visualizer example
    void sample_curved_edge(const vd_type::edge_type& edge,
                            std::vector< boost::polygon::point_data<vd_type::coordinate_type> >* sampled_edge);
    boost::polygon::point_data<int> retrieve_point(const vd_type::cell_type& cell);
    boost::polygon::segment_data<int> retrieve_segment(const vd_type::cell_type& cell);

    const vd_type &vd;
    const std::vector< boost::polygon::point_data<int> > point_data;
    const std::vector< boost::polygon::segment_data<int> > segment_data;
    const std::vector<uint32_t> components;
    uint32_t num_components;
    const std::vector<cv::Vec3b> colors;

    Path* component_paths;
    std::vector<Point> all_points;
    std::unordered_map<Point, uint32_t> point_map;
};

} /* Depixelize */

#endif /* DEPIXELIZE_SPLINE_OPTIMIZER_HPP */
