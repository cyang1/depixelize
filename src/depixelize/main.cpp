#include <cstdlib>
#include <cstdint>
#include <getopt.h>
#include <iostream>
#include <vector>

#include <boost/filesystem/path.hpp>
#include <boost/polygon/voronoi.hpp>
#include <opencv2/opencv.hpp>

#include "pixel_grid.hpp"
#include "spline_optimizer.hpp"
#include "render/svg_renderer.hpp"

using namespace Depixelize;

/**
 * Struct of the program options.
 */
struct Options
{
    // not allocated, pointed it to something static
    const char* input_filename;
    // not allocated, pointed it to something static
    const char* output_filename;
    // output file scale (defaults to 4)
    int output_scale;
};

/**
 * Prints out a nice usage message for this program
 */
static void print_usage(const char* progname)
{
    printf("Usage: %s [options] input_filename\n", progname);
    printf("Program Options:\n");
    printf("  -o  --output <FILENAME>    Puts the resulting image into <FILENAME> (default: [input_filename]_out)\n");
    printf("  -s  --scale  <INT>         Scale the original image by this factor (default: 4)\n");
    printf("  -h  --help                 This message\n");
}

/**
 * Parses args into an Options struct. Returns true on success, false on failure.
 */
static bool parse_args(Options* options, int argc, char* argv[])
{
    if (argc < 2) {
        print_usage(argv[0]);
        return false;
    }

    options->input_filename = NULL;
    options->output_filename = NULL;
    options->output_scale = 4;

    // parse commandline options ////////////////////////////////////////////
    int opt;
    static struct option long_options[] = {
        {"help",     0, 0,  'h'},
        {"output",   1, 0,  'o'},
        {"scale",    1, 0,  's'},
        {0 ,0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "o:s:h", long_options, NULL)) != EOF) {
        switch (opt) {
        case 'o':
            options->output_filename = optarg;
            break;
        case 's':
            options->output_scale = atoi(optarg);
            break;
        case 'h':
        default:
            print_usage(argv[0]);
            return false;
        }
    }
    // end parsing of commandline options //////////////////////////////////////

    if (optind + 1 > argc) {
        fprintf(stderr, "Error: missing input filename\n");
        print_usage(argv[0]);
        return false;
    }

    options->input_filename = argv[optind];

    if (options->output_filename == NULL) {
        boost::filesystem::path filename = boost::filesystem::path(options->input_filename).stem();
        filename += "_out";
        options->output_filename = filename.c_str();
    }

    return true;
}

int main(int argc, char* argv[])
{
    Options opt;
    if (!parse_args(&opt, argc, argv)) {
        return 1;
    }

    cv::Mat img = cv::imread(opt.input_filename);
    cv::cvtColor(img, img, CV_BGR2YCrCb);

    PixelGrid pg;
    pg.initialize(img);
    pg.make_planar();

    std::vector< boost::polygon::point_data<int> > points;
    std::vector< boost::polygon::segment_data<int> > edges;
    std::vector<cv::Vec3b> colors;
    std::vector<uint32_t> components;
    uint32_t num_components = pg.get_data(&points, &edges, &colors, &components);

    // std::cout << points.size() << std::endl;
    // for (auto it = points.begin(); it != points.end(); ++it) {
    //     std::cout << it->y() << " " << it->x() << std::endl;
    // }

    // std::cout << edges.size() << std::endl;
    // for (auto it = edges.begin(); it != edges.end(); ++it) {
    //     std::cout << it->low().y() << " " << it->low().x() << " " << it->high().y() << " " << it->high().x() << std::endl;
    // }

    boost::polygon::voronoi_diagram<double> vd;
    construct_voronoi(points.begin(), points.end(), edges.begin(), edges.end(), &vd);

    SplineOptimizer sp_opt(vd, points, edges, components, num_components, colors);
    sp_opt.initialize();
    sp_opt.optimize_splines();
    std::vector<Shape> shapes = sp_opt.make_shapes();

    // pg.get_data above doubles the coordinate locations so that they stay integers
    SVGRenderer svg_renderer(img.cols * 2, img.rows * 2, opt.output_scale / 2.0);
    svg_renderer.render(svg_renderer.make_full_filename(opt.output_filename), shapes);

    return 0;
}
