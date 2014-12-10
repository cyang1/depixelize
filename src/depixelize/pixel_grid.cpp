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

#define DIR_TO_IDX(dir) (CHAR_BIT * sizeof(uint8_t) - __builtin_clz(dir) - 1)
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
const int DIRECTIONS[NUM_DIRS][DIM] = { { 0, 1 }, { 1, 1 }, { 1, 0 }, { 1, -1 }, { 0, -1 }, { -1, -1 }, { -1, 0 }, { -1, 1 } };

void PixelGrid::initialize(const cv::Mat img)
{
    this->img = img;
    this->neighbors = new uint8_t[img.rows * img.cols];
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
        uint8_t* row = this->img.ptr<uint8_t>(r);
        for (int c = 0; c < this->img.cols; c++) {
            uint8_t cur_bitmask = this->neighbors[r * this->img.cols + c];
            auto cur_edges = this->neighbor_edges[r * this->img.cols + c];

            // Take advantage of symmetry by only going through forward edges
            for (int i = 0; i < NUM_FORWARD_DIRS; i++) {
                if (cur_bitmask & IDX_TO_DIR(i)) {
                    if (dissimilar_colors(&row[ch * c], this->img.at<cv::Vec3b>(r + DIRECTIONS[i][0], c + DIRECTIONS[i][1]).val)) {
                        remove_edge(cur_edges[i], i);
                    }
                }
            }
        }
    }

    std::vector< std::pair<boost::polygon::segment_data<int>, uint8_t> > edges_to_remove;
    for (int r = 0; r < this->img.rows - 1; r++) {
        for (int c = 0; c < this->img.cols - 1; c++) {
            uint8_t tl_bitmask = this->neighbors[r * this->img.cols + c];
            uint8_t tr_bitmask = this->neighbors[r * this->img.cols + c + 1];
            uint8_t bl_bitmask = this->neighbors[(r + 1) * this->img.cols + c];

            boost::polygon::segment_data<int> d1 = this->neighbor_edges[r * this->img.cols + c][IDX_BR];
            boost::polygon::segment_data<int> d2 = this->neighbor_edges[r * this->img.cols + c + 1][IDX_BL];

            // Both diagonals exist, need to resolve to make planar
            if (DIAGONALS_EXIST(tl_bitmask, tr_bitmask)) {
                if (BOX_EXISTS(tl_bitmask, tr_bitmask, bl_bitmask)) {
                    edges_to_remove.push_back({ d1, IDX_BR });
                    edges_to_remove.push_back({ d2, IDX_BL });
                } else {
                    int r1 = rate_diagonal(d1);
                    int r2 = rate_diagonal(d2);
                    if (r1 <= r2) edges_to_remove.push_back({ d1, IDX_BR });
                    if (r2 <= r1) edges_to_remove.push_back({ d2, IDX_BL });
                }
            }
        }
    }

    // Remove these edges in a second pass because edge removal heuristics
    // depend on the original graph
    for (auto &p : edges_to_remove) {
        remove_edge(p.first, p.second);
    }
}

int PixelGrid::get_data(std::vector< boost::polygon::point_data<int> >* points,
                         std::vector< boost::polygon::segment_data<int> >* edges,
                         std::vector<cv::Vec3b>* colors,
                         std::vector<uint32_t>* components)
{
    std::vector<cv::Vec3b> pt_colors;
    std::vector<cv::Vec3b> edge_colors;
    std::vector<uint64_t> pt_components;
    std::vector<uint64_t> edge_components;

    points->clear();
    edges->clear();

    uint32_t *img_components = new uint32_t[this->img.rows * this->img.cols];
    uint32_t num_components = this->get_components(img_components);

    // for (int r = 0; r < this->img.rows; r++) {
    //     for (int c = 0; c < this->img.cols; c++) {
    //         printf("%2d ", img_components[r * this->img.cols + c]);
    //     }
    //     std::cout << std::endl;
    // }

    // Yes, accessing colors one at a time is inefficient. Probably not much
    // less efficient than constructing a cv::Vec3b each time though.
    for (int r = 0; r < this->img.rows; r++) {
        for (int c = 0; c < this->img.cols; c++) {
            uint8_t cur_bitmask = this->neighbors[r * this->img.cols + c];
            if (cur_bitmask == 0) {
                points->push_back(boost::polygon::point_data<int>(2 * c, 2 * r));
                pt_colors.push_back(this->img.at<cv::Vec3b>(r, c));
                pt_components.push_back(img_components[r * this->img.cols + c]);
                continue;
            }
            auto cur_edges = this->neighbor_edges[r * this->img.cols + c];
            for (int i = 0; i < NUM_FORWARD_DIRS; i++) {
                if (cur_bitmask & IDX_TO_DIR(i)) {
                    int next_r = r + DIRECTIONS[i][0];
                    int next_c = c + DIRECTIONS[i][1];

                    boost::polygon::segment_data<int> e = cur_edges[i];
                    boost::polygon::point_data<int> p1(2 * e.low().y(), 2 * e.low().x()),
                        p2(2 * e.high().y(), 2 * e.high().x());
                    boost::polygon::segment_data<int> e1, e2;
                    this->split_edge(boost::polygon::segment_data<int>(p1, p2), &e1, &e2);

                    edges->push_back(e1);
                    edges->push_back(e2);
                    edge_colors.push_back(this->img.at<cv::Vec3b>(r, c));
                    edge_colors.push_back(this->img.at<cv::Vec3b>(next_r, next_c));
                    edge_components.push_back(img_components[r * this->img.cols + c]);
                    edge_components.push_back(img_components[next_r * this->img.cols + next_c]);
                }
            }
        }
    }

    colors->insert(colors->end(), pt_colors.begin(), pt_colors.end());
    colors->insert(colors->end(), edge_colors.begin(), edge_colors.end());
    components->insert(components->end(), pt_components.begin(), pt_components.end());
    components->insert(components->end(), edge_components.begin(), edge_components.end());

    return num_components;
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

        uint8_t bitmask = this->neighbors[r * w + c];
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

        uint8_t bitmask = this->neighbors[r * w + c];
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

    uint8_t low_mask = this->neighbors[d.low().x() * w + d.low().y()];
    uint8_t high_mask = this->neighbors[d.high().x() * w + d.high().y()];

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

inline void PixelGrid::split_edge(boost::polygon::segment_data<int> e,
                                  boost::polygon::segment_data<int>* out_1,
                                  boost::polygon::segment_data<int>* out_2)
{
    boost::polygon::point_data<int> pt1(e.low().x(), e.low().y());
    boost::polygon::point_data<int> pt2(e.high().x(), e.high().y());
    boost::polygon::point_data<int> mid_pt((pt1.x() + pt2.x()) / 2, (pt1.y() + pt2.y()) / 2);

    *out_1 = boost::polygon::segment_data<int>(pt1, mid_pt);
    *out_2 = boost::polygon::segment_data<int>(mid_pt, pt2);
}

inline boost::polygon::segment_data<int> PixelGrid::reverse_edge(boost::polygon::segment_data<int> e)
{
    return boost::polygon::segment_data<int>(e.high(), e.low());
}

void PixelGrid::remove_edge(boost::polygon::segment_data<int> e, int edge_idx)
{
    // No need to actually remove edges from the neighbor_edges structure
    boost::polygon::point_data<int> pt1 = e.low();
    boost::polygon::point_data<int> pt2 = e.high();

    this->neighbors[pt1.x() * this->img.cols + pt1.y()] &= ~(IDX_TO_DIR(edge_idx));
    this->neighbors[pt2.x() * this->img.cols + pt2.y()] &= ~(IDX_TO_DIR(MIRROR_EDGE(edge_idx)));
}

uint32_t PixelGrid::get_components(uint32_t *components)
{
    bool *seen = new bool[this->img.rows * this->img.cols]();
    uint32_t component_num = 0;

    for (int r = 0; r < this->img.rows; r++) {
        for (int c = 0; c < this->img.cols; c++) {
            int pos = r * this->img.cols + c;
            if (!seen[pos]) {
                std::queue<int> q;
                q.push(pos);
                components[pos] = component_num;
                seen[pos] = true;

                while (!q.empty()) {
                    int cur_pos = q.front();
                    q.pop();

                    uint8_t bitmask = this->neighbors[cur_pos];
                    for (int i = 0; i < NUM_DIRS; i++) {
                        if (bitmask & IDX_TO_DIR(i)) {
                            int new_pos = cur_pos + (DIRECTIONS[i][0] * this->img.cols + DIRECTIONS[i][1]);

                            if (!seen[new_pos]) {
                                components[new_pos] = component_num;
                                seen[new_pos] = true;
                                q.push(new_pos);
                            }
                        }
                    }
                }

                component_num++;
            }
        }
    }

    delete[] seen;
    return component_num;
}

} /* Depixelize */
