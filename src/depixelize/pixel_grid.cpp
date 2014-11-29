#include "pixel_grid.hpp"
#include "color_util.hpp"

#include <climits>
#include <set>
#include <queue>
#include <algorithm>

#define NUM_DIRS 8
#define DIM 2

#define NUM_FORWARD_DIRS 4
#define IDX_R  0
#define IDX_BR 1
#define IDX_B  2
#define IDX_BL 3
#define IDX_L  4
#define IDX_TL 5
#define IDX_T  6
#define IDX_TR 7

#define DIR_R  (1 << IDX_R)
#define DIR_BR (1 << IDX_BR)
#define DIR_B  (1 << IDX_B)
#define DIR_BL (1 << IDX_BL)
#define DIR_L  (1 << IDX_L)
#define DIR_TL (1 << IDX_TL)
#define DIR_T  (1 << IDX_T)
#define DIR_TR (1 << IDX_TR)

#define DIR_TO_IDX(dir) (CHAR_BIT * sizeof(int) - __builtin_clz(dir) - 1)
#define IDX_TO_DIR(idx) (1 << idx)

#define DIAGONALS_EXIST(tl, tr) (((tl) & DIR_BR) && ((tr) & DIR_BL))
#define BOX_EXISTS(tl, tr, bl) (((tl) & DIR_R) && ((tl) & DIR_B) && ((tr) & DIR_B) && ((bl) & DIR_R))

#define MIRROR_EDGE(num) (((num) + NUM_DIRS / 2) % NUM_DIRS)
#define VALENCE(neighbor_mask) __builtin_popcount(neighbor_mask)

#define SPARSE_PIXEL_WINDOW 8

namespace Depixelize {

/*
 * These are row, col. They are also ordered explicitly as follows:
 *  5  6  7
 *   \ | /
 * 4 - x - 0
 *   / | \
 *  3  2  1
 */
int DIRECTIONS[NUM_DIRS][DIM] = { { 0, 1 }, { 1, 1 }, { 1, 0 }, { 1, -1 }, { 0, -1 }, { -1, -1 }, { -1, 0 }, { -1, 1 } };

void PixelGrid::initialize(const cv::Mat img)
{
    this->img = img;
    this->neighbors = new unsigned char[img.rows * img.cols];
    this->neighbor_edges = new boost::polygon::segment_data<int>[img.rows * img.cols][8];

    // Initialize the neighbors bitmasks
    std::fill(this->neighbors, this->neighbors + img.rows * img.cols, 0xFF);
    if (img.rows == 1) {
        std::fill(this->neighbors, this->neighbors + img.cols, 0x11); // 0b00010001, no top or bottom neighbors
    } else {
        std::fill(this->neighbors, this->neighbors + img.cols, 0x8F); // 0b00011111, no top neighbors
        std::fill(this->neighbors + (img.rows - 1) * img.cols, this->neighbors + img.rows * img.cols, 0xF1); // 0b11110001, no bottom neighbors
    }
    for (int i = 0; i < img.rows; i++) {
        this->neighbors[i * img.cols] &= 0xC7;// 0b11000111, no left neighbors
        this->neighbors[i * img.cols + (img.cols - 1)] &= 0x7C;// 0b01111100, no right neighbors
    }

    // Initialize the neighbor_edges properly
    // NOTE: This data structure is actually entirely static, it is only created
    // for convenience
    for (int r = 0; r < img.rows; r++) {
        for (int c = 0; c < img.cols; c++) {
            boost::polygon::point_data<int> cur_pt(r, c);
            for (int i = 0; i < NUM_DIRS; i++) {
                boost::polygon::point_data<int> this_pt(r + DIRECTIONS[i][0], c + DIRECTIONS[i][1]);
                this->neighbor_edges[r * img.cols + c][i] = boost::polygon::segment_data<int>(cur_pt, this_pt);
            }
        }
    }
}

void PixelGrid::make_planar()
{
    const int ch = this->img.channels();

    for (int r = 0; r < this->img.rows; r++) {
        unsigned char* row = this->img.ptr<unsigned char>(r);
        for (int c = 0; c < this->img.cols; c++) {
            unsigned char cur_bitmask = this->neighbors[r * this->img.cols + c];
            auto cur_edges = this->neighbor_edges[r * this->img.cols + c];

            // Take advantage of symmetry by only going through forward edges
            for (int i = 0; i < NUM_FORWARD_DIRS; i++) {
                if (cur_bitmask & IDX_TO_DIR(i)) {
                    boost::polygon::segment_data<int> e = cur_edges[i];
                    boost::polygon::point_data<int> other = e.high();
                    if (dissimilar_colors(&row[ch * c], this->img.at<cv::Vec3b>(other.x(), other.y()).val)) {
                        remove_edge(e, i);
                    }
                }
            }
        }
    }

    for (int r = 0; r < this->img.rows - 1; r++) {
        for (int c = 0; c < this->img.cols - 1; c++) {
            unsigned char tl_bitmask = this->neighbors[r * this->img.cols + c];
            unsigned char tr_bitmask = this->neighbors[r * this->img.cols + c + 1];
            unsigned char bl_bitmask = this->neighbors[(r + 1) * this->img.cols + c];

            boost::polygon::segment_data<int> d1 = this->neighbor_edges[r * this->img.cols + c][IDX_BR];
            boost::polygon::segment_data<int> d2 = this->neighbor_edges[r * this->img.cols + c + 1][IDX_BL];

            // Both diagonals exist, need to resolve to make planar
            if (DIAGONALS_EXIST(tl_bitmask, tr_bitmask)) {
                if (BOX_EXISTS(tl_bitmask, tr_bitmask, bl_bitmask)) {
                    remove_edge(d1, IDX_BR);
                    remove_edge(d2, IDX_BL);
                } else {
                    int r1 = rate_diagonal(d1);
                    int r2 = rate_diagonal(d2);
                    if (r1 <= r2) remove_edge(d1, IDX_BR);
                    if (r2 <= r1) remove_edge(d2, IDX_BL);
                }
            }
        }
    }
}

void PixelGrid::get_data(std::vector< boost::polygon::point_data<int> > &points,
                         std::vector< boost::polygon::segment_data<int> > &edges)
{
    points.clear();
    edges.clear();

    for (int r = 0; r < this->img.rows; r++) {
        for (int c = 0; c < this->img.cols; c++) {
            unsigned char cur_bitmask = this->neighbors[r * this->img.cols + c];
            if (cur_bitmask == 0) {
                points.push_back(boost::polygon::point_data<int>(r, c));
                continue;
            }
            for (int i = 0; i < NUM_FORWARD_DIRS; i++) {
                if (cur_bitmask & IDX_TO_DIR(i)) {
                    edges.push_back(this->neighbor_edges[r * this->img.cols + c][i]);
                }
            }
        }
    }
}

int PixelGrid::rate_diagonal(boost::polygon::segment_data<int> d)
{
    return curves_heuristic(d) + sparse_pixels_heuristic(d) + islands_heuristic(d);
}

int PixelGrid::curves_heuristic(boost::polygon::segment_data<int> d)
{
    const int w = this->img.cols;

    std::set< boost::polygon::segment_data<int> > seen_edges;
    seen_edges.insert(d);
    seen_edges.insert(reverse_edge(d));
    std::queue< boost::polygon::point_data<int> > nodes;
    nodes.push(d.low());
    nodes.push(d.high());

    while (!nodes.empty()) {
        boost::polygon::point_data<int> cur_pt = nodes.front();
        nodes.pop();

        // x and y are r and c in our system
        int r = cur_pt.x();
        int c = cur_pt.y();

        unsigned char bitmask = this->neighbors[r * w + c];
        auto edges = this->neighbor_edges[r * w + c];
        if (VALENCE(bitmask) != 2) {
            continue;
        }
        for (int i = 0; i < NUM_DIRS; i++) {
            if (bitmask & IDX_TO_DIR(i)) {
                boost::polygon::segment_data<int> cur_edge = edges[i];
                if (seen_edges.count(cur_edge) == 0) {
                    seen_edges.insert(cur_edge);
                    seen_edges.insert(reverse_edge(cur_edge));

                    // Low is always the current node
                    nodes.push(cur_edge.high());
                }
            }
        }
    }

    // The curves heuristic votes with the number of edges in the sequence
    return seen_edges.size() / 2;
}

int PixelGrid::sparse_pixels_heuristic(boost::polygon::segment_data<int> d)
{
    const int w = this->img.cols;

    int rmin = std::min(d.low().x(), d.high().x());
    int rmax = std::max(d.low().x(), d.high().x());
    int cmin = std::min(d.low().y(), d.high().y());
    int cmax = std::min(d.low().y(), d.high().y());

    int bb_rmin = std::max(rmin - SPARSE_PIXEL_WINDOW / 2 - 1, 0);
    int bb_rmax = std::min(rmax + SPARSE_PIXEL_WINDOW / 2 - 1, this->img.rows - 1);
    int bb_cmin = std::max(cmin - SPARSE_PIXEL_WINDOW / 2 - 1, 0);
    int bb_cmax = std::min(cmax + SPARSE_PIXEL_WINDOW / 2 - 1, this->img.cols - 1);

    std::set< boost::polygon::point_data<int> > seen_nodes;
    seen_nodes.insert(d.low());
    seen_nodes.insert(d.high());
    std::queue< boost::polygon::point_data<int> > frontier;
    frontier.push(d.low());
    frontier.push(d.high());

    while (!frontier.empty()) {
        boost::polygon::point_data<int> cur_pt = frontier.front();
        frontier.pop();

        // x and y are r and c in our system
        int r = cur_pt.x();
        int c = cur_pt.y();

        unsigned char bitmask = this->neighbors[r * w + c];
        auto edges = this->neighbor_edges[r * w + c];
        for (int i = 0; i < NUM_DIRS; i++) {
            if (bitmask & IDX_TO_DIR(i)) {
                boost::polygon::point_data<int> next_pt = edges[i].high();

                if (seen_nodes.count(next_pt) == 0 && in_bounds(next_pt, bb_rmin, bb_rmax, bb_cmin, bb_cmax)) {
                    seen_nodes.insert(next_pt);
                    frontier.push(next_pt);
                }
            }
        }
    }

    return -seen_nodes.size();
}

int PixelGrid::islands_heuristic(boost::polygon::segment_data<int> d)
{
    const int w = this->img.cols;

    unsigned char low_mask = this->neighbors[d.low().x() * w + d.low().y()];
    unsigned char high_mask = this->neighbors[d.high().x() * w + d.high().y()];

    if (VALENCE(low_mask) == 1 || VALENCE(high_mask) == 1) {
        return 5;
    }
    return 0;
}

inline bool PixelGrid::in_bounds(boost::polygon::point_data<int> p, int rmin, int rmax, int cmin, int cmax)
{
    int r = p.x();
    int c = p.y();
    return rmin <= r && r <= rmax && cmin <= c && c <= cmax;
}

inline boost::polygon::segment_data<int> PixelGrid::reverse_edge(boost::polygon::segment_data<int> e)
{
    return boost::polygon::segment_data<int>(e.high(), e.low());
}

void PixelGrid::remove_edge(boost::polygon::segment_data<int> e, int edge_num)
{
    // No need to actually remove edges from the neighbor_edges structure
    boost::polygon::point_data<int> pt1 = e.low();
    boost::polygon::point_data<int> pt2 = e.high();

    this->neighbors[pt1.x() * this->img.cols + pt1.y()] &= ~(1 << edge_num);
    this->neighbors[pt2.x() * this->img.cols + pt2.y()] &= ~(1 << MIRROR_EDGE(edge_num));
}

} /* Depixelize */
