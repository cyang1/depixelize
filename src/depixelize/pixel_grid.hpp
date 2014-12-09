#ifndef DEPIXELIZE_PIXEL_GRID_HPP
#define DEPIXELIZE_PIXEL_GRID_HPP

#include <cstdint>
#include <vector>

#include <opencv2/opencv.hpp>
#include <boost/polygon/point_data.hpp>
#include <boost/polygon/segment_data.hpp>

namespace Depixelize {

using NeighborEdgeList = boost::polygon::segment_data<int> (*)[8];

class PixelGrid
{
public:

    PixelGrid() : img(), neighbors(NULL), neighbor_edges(NULL) { };

    virtual ~PixelGrid()
    {
        if (neighbors != NULL) {
            delete[] neighbors;
            delete[] neighbor_edges;
        }
    }

    void initialize(const cv::Mat img);

    void make_planar();

    int get_data(std::vector< boost::polygon::point_data<int> >* points,
                 std::vector< boost::polygon::segment_data<int> >* edges,
                 std::vector<cv::Vec3b>* colors,
                 std::vector<uint32_t>* components);

private:

    // Curves, sparse pixels, and islands heuristic from paper
    int rate_diagonal(boost::polygon::segment_data<int> d);

    int curves_heuristic(boost::polygon::segment_data<int> d);
    int sparse_pixels_heuristic(boost::polygon::segment_data<int> d);
    int islands_heuristic(boost::polygon::segment_data<int> d);

    inline bool in_bounds(boost::polygon::point_data<int> p, int rmin, int rmax, int cmin, int cmax);

    inline void split_edge(boost::polygon::segment_data<int> e,
                           boost::polygon::segment_data<int>* out_1,
                           boost::polygon::segment_data<int>* out_2);

    inline boost::polygon::segment_data<int> reverse_edge(boost::polygon::segment_data<int> e);

    // edge_num must correspond to the number of the edge for the source point
    void remove_edge(boost::polygon::segment_data<int> e, int edge_num);

    uint32_t get_components(uint32_t *components);

    cv::Mat img;
    uint8_t *neighbors;
    NeighborEdgeList neighbor_edges;
};

} /* Depixelize */

#endif /* DEPIXELIZE_PIXEL_GRID_HPP */
