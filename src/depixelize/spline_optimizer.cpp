#include "spline_optimizer.hpp"
#include "color_util.hpp"
#include "math_util.hpp"
#include "voronoi_visual_utils.hpp"

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

namespace Depixelize {

void SplineOptimizer::initialize()
{
    std::unordered_multimap<uint32_t, EdgeRef> adjacent_edges;
    std::unordered_set<Edge> seen_edges;

    std::vector<EdgeRef> *component_edges = new std::vector<EdgeRef>[this->num_components]();

    for (vd_type::const_edge_iterator it = this->vd.edges().begin(); it != this->vd.edges().end(); ++it) {
        // Cannot be visible if it isn't primary
        if (!it->is_primary()) {
            continue;
        }

        // Also assume that all primary edges are finite...
        if (!it->is_finite()) {
            std::cout << "WARNING: Primary, infinite edge..." << std::endl;
        }

        uint64_t cell_idx1 = it->cell()->source_index();
        uint64_t cell_idx2 = it->twin()->cell()->source_index();

        // The two sides of the edge have different components, its visible
        if (this->components[cell_idx1] != this->components[cell_idx2]) {
            Edge cur_edge = Edge(*it);
            // This may not work if we have a curved and a straight edge with
            // the same endpoints, but I don't think that will happen...
            if (seen_edges.count(cur_edge) > 0) {
                continue;
            }
            seen_edges.insert(cur_edge);

            bool is_shading_edge = shading_edge(this->colors[cell_idx1].val, this->colors[cell_idx2].val);

            using point_type = boost::polygon::point_data<vd_type::coordinate_type>;
            std::vector<point_type> samples;
            point_type vertex0(it->vertex0()->x(), it->vertex0()->y());
            point_type vertex1(it->vertex1()->x(), it->vertex1()->y());
            samples.push_back(vertex0);
            samples.push_back(vertex1);
            if (it->is_curved()) {
                this->sample_curved_edge(*it, &samples);
            }

            std::vector<uint32_t> path_points;
            std::vector<EdgeRef> path_edges;
            for (uint32_t i = 0; i < samples.size(); i++) {
                uint32_t pt_idx;
                Point cur_pt = Point(samples[i].x(), samples[i].y());
                auto pt_idx_it = this->point_map.find(cur_pt);
                if (pt_idx_it == this->point_map.end()) {
                    pt_idx = this->all_points.size();
                    this->all_points.push_back(cur_pt);
                    this->point_map.insert({ cur_pt, pt_idx });
                } else {
                    pt_idx = pt_idx_it->second;
                }

                path_points.push_back(pt_idx);
                if (i != 0) path_edges.push_back({ path_points[i - 1], path_points[i], is_shading_edge });
            }

            for (uint32_t i = 0; i < samples.size(); i++) {
                if (i != 0) adjacent_edges.insert({ path_points[i], path_edges[i - 1] });
                if (i != samples.size() - 1) {
                    adjacent_edges.insert({ path_points[i], path_edges[i] });
                    component_edges[this->components[cell_idx1]].push_back(path_edges[i]);
                    component_edges[this->components[cell_idx2]].push_back(path_edges[i]);
                }
            }
        }
    }

    int hist[] = { 0, 0, 0, 0, 0, 0 };
    for (uint32_t i = 0; i < this->all_points.size(); i++) {
        hist[adjacent_edges.count(i) - 1]++;
    }
    std::cout << "Number of edges with given valences" << std::endl;
    for (int i = 0; i < 6; i++) {
        std::cout << i + 1 << ": " << hist[i] << std::endl;
    }

    cv::Mat im1 = cv::Mat::zeros(39 * 20, 49 * 20, CV_8UC3);
    for (uint32_t i = 0; i < this->all_points.size(); i++) {
        auto range = adjacent_edges.equal_range(i);
        for_each(
            range.first,
            range.second,
            [&](std::unordered_multimap<uint32_t, EdgeRef>::value_type &p){
                Point pt1 = this->all_points[p.second.idx1];
                Point pt2 = this->all_points[p.second.idx2];
                cv::line(im1, cv::Point(10 * pt1.x, 10 * pt1.y), cv::Point(10 * pt2.x, 10 * pt2.y), cv::Scalar(0, 0, 255), 1, 8);
            }
        );
    }
    cv::imshow("Image", im1);
    cv::waitKey(0);

    this->component_paths = new Path[this->num_components]();

    for (uint32_t i = 0; i < this->num_components; i++) {
        join_edges(&this->component_paths[i], component_edges[i], adjacent_edges);
    }

    cv::Mat im2 = cv::Mat::zeros(39 * 20, 49 * 20, CV_8UC3);
    for (uint32_t i = 0; i < this->num_components; i++) {
        Path &cur_path = this->component_paths[i];

        for (uint32_t j = 0; j < cur_path.size(); j++) {
            PointRef prev_ptref = j == 0 ? cur_path.back() : cur_path[j - 1];
            PointRef cur_ptref = cur_path[j];
            Point pt1 = this->all_points[prev_ptref.idx];
            Point pt2 = this->all_points[cur_ptref.idx];

            cv::line(im2, cv::Point(10 * pt1.x, 10 * pt1.y), cv::Point(10 * pt2.x, 10 * pt2.y), cv::Scalar(255.0 * i / this->num_components, 255, 192), 1, 8);
        }
    }
    cv::cvtColor(im2, im2, CV_HSV2BGR);
    cv::imshow("Image", im2);
    cv::waitKey(0);

    delete[] component_edges;
}

void SplineOptimizer::optimize_splines()
{

}

std::vector<Shape> SplineOptimizer::make_shapes()
{
    std::vector<ColorPoint> *component_colors = new std::vector<ColorPoint>[this->num_components]();

    for (vd_type::const_cell_iterator it = this->vd.cells().begin(); it != this->vd.cells().end(); ++it) {
        Point centroid(0, 0);
        uint32_t num_points = 0;
        auto *edge = it->incident_edge();
        if (edge == NULL) {
            continue;
        }

        // Actually calculates the barycenter. Also leaves out infinite
        // edges when those should be clipped to the edge of the image
        do {
            if (edge->vertex0()) {
                centroid.x += edge->vertex0()->x();
                centroid.y += edge->vertex0()->y();
                num_points += 1;
            }

            if (edge->vertex1()) {
                centroid.x += edge->vertex1()->x();
                centroid.y += edge->vertex1()->y();
                num_points += 1;
            }

            edge = edge->next();
        } while (edge != it->incident_edge());

    }

    std::vector<Shape> ret;
    for (uint32_t i = 0; i < this->num_components; i++) {
        std::deque<Point> actual_path;
        for (auto &p : this->component_paths[i]) {
            actual_path.push_back(this->all_points[p.idx]);
        }

        ret.push_back(Shape(actual_path, component_colors[i]));
    }

    delete[] component_colors;

    std::sort(ret.begin(), ret.end());
    std::reverse(ret.begin(), ret.end());
    return ret;
}

void SplineOptimizer::join_edges(Path *path, const std::vector<EdgeRef> &edges,
    const std::unordered_multimap<uint32_t, EdgeRef> &adjacent_edges)
{
    std::unordered_multimap<uint32_t, EdgeRef> adj_list;
    for (auto &it : edges) {
        adj_list.insert({ it.idx1, it });
        adj_list.insert({ it.idx2, { it.idx2, it.idx1, it.is_shading_edge } });
    }

    path->push_back({ edges[0].idx1, true });
    std::unordered_set<uint32_t> used_pts;
    used_pts.insert(edges[0].idx1);

    EdgeRef cur_edge = edges[0];
    uint32_t next_pt = edges[0].idx2;
    bool next_should_optimize = true;
    do {
        path->push_back({ next_pt, next_should_optimize });
        used_pts.insert(next_pt);

        // There can be one, in that case this might break...
        auto range = adj_list.equal_range(next_pt);
        next_pt = all_points.size();
        for (auto it = range.first; it != range.second; ++it) {
            if (used_pts.count(it->second.idx2) == 0) {
                next_should_optimize = should_optimize(cur_edge, it->second, adjacent_edges);
                cur_edge = it->second;
                next_pt = cur_edge.idx2;
                break;
            }
        }
    } while (next_pt != all_points.size());

    // Have to do the last one separately after we know what the last point
    // in the path is
    auto last_range = adj_list.equal_range(path->back().idx);
    EdgeRef last_edge;
    for (auto it = last_range.first; it != last_range.second; ++it) {
        if (it->second.idx2 == path->front().idx) {
            last_edge = it->second;
            break;
        }
    }
    path->front().can_optimize = should_optimize(last_edge, edges[0], adjacent_edges);
}

// Invariant: e1.idx2 == e2.idx1
bool SplineOptimizer::should_optimize(EdgeRef e1, EdgeRef e2,
    const std::unordered_multimap<uint32_t, EdgeRef> &adjacent_edges)
{
    assert(e1.idx2 == e2.idx1);

    uint32_t common_pt = e1.idx2;
    uint32_t num_edges_around = adjacent_edges.count(common_pt);
    if (num_edges_around <= 2) {
        return true;
    } else if (num_edges_around >= 4) {
        // If there are at least 4 edges, then this point is fixed.
        return false;
    }

    EdgeRef last_edge;
    auto range = adjacent_edges.equal_range(common_pt);
    for (auto it = range.first; it != range.second; ++it) {
        if (it->second != e1 && it->second != e2) {
            last_edge = it->second;
            break;
        }
    }
    uint32_t last_pt = last_edge.idx1 == common_pt ? last_edge.idx2 : last_edge.idx1;

    // If we have 2 contour edges and 1 shading edge, then we can resolve the
    // ambiguity easily: true iff both e1 and e2 are contour edges
    int num_shading_edges = last_edge.is_shading_edge + e1.is_shading_edge + e2.is_shading_edge;
    if (num_shading_edges == 1) {
        return last_edge.is_shading_edge;
    }

    // Everything else failed, we have to measure the angles
    Point e1_vec(this->all_points[e1.idx1].x - this->all_points[common_pt].x,
        this->all_points[e1.idx1].y - this->all_points[common_pt].y);
    Point e2_vec(this->all_points[e2.idx2].x - this->all_points[common_pt].x,
        this->all_points[e2.idx2].y - this->all_points[common_pt].y);
    Point last_vec(this->all_points[last_pt].x - this->all_points[common_pt].x,
        this->all_points[last_pt].y - this->all_points[common_pt].y);

    // Technically, we want the pair closest to 180 degrees. However,
    // vector_angle will always return the principle angle (<180), so the
    // largest angle will also be the one closest to 180
    double e1_e2_angle = vector_angle(e1_vec, e2_vec);
    return e1_e2_angle > vector_angle(e2_vec, last_vec) && e1_e2_angle > vector_angle(e1_vec, last_vec);
}

void SplineOptimizer::sample_curved_edge(const vd_type::edge_type& edge,
    std::vector< boost::polygon::point_data<vd_type::coordinate_type> >* sampled_edge)
{
  boost::polygon::point_data<int> point = edge.cell()->contains_point() ?
      retrieve_point(*edge.cell()) :
      retrieve_point(*edge.twin()->cell());
  boost::polygon::segment_data<int> segment = edge.cell()->contains_point() ?
      retrieve_segment(*edge.twin()->cell()) :
      retrieve_segment(*edge.cell());
  boost::polygon::voronoi_visual_utils<vd_type::coordinate_type>::discretize(
      point, segment, 0.05, sampled_edge);
}

boost::polygon::point_data<int> SplineOptimizer::retrieve_point(const vd_type::cell_type& cell)
{
  vd_type::cell_type::source_index_type index = cell.source_index();
  vd_type::cell_type::source_category_type category = cell.source_category();
  if (category == boost::polygon::SOURCE_CATEGORY_SINGLE_POINT) {
    return this->point_data[index];
  }
  index -= this->point_data.size();
  if (category == boost::polygon::SOURCE_CATEGORY_SEGMENT_START_POINT) {
    return low(this->segment_data[index]);
  } else {
    return high(this->segment_data[index]);
  }
}

boost::polygon::segment_data<int> SplineOptimizer::retrieve_segment(const vd_type::cell_type& cell)
{
  vd_type::cell_type::source_index_type index = cell.source_index() - this->point_data.size();
  return this->segment_data[index];
}

} /* Depixelize */
