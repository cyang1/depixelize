#ifndef DEPIXELIZE_COLOR_UTIL_HPP
#define DEPIXELIZE_COLOR_UTIL_HPP

#include <cstdint>
#include <opencv2/opencv.hpp>

namespace Depixelize {

inline cv::Vec3b yuv2rgb(const cv::Vec3b &yuv)
{
    // Hackity hack hack hack
    cv::Mat M(1, 1, CV_8UC3, cv::Scalar::all(0));
    M.at<cv::Vec3b>(0, 0) = yuv;
    cv::cvtColor(M, M, CV_YCrCb2RGB);
    return M.at<cv::Vec3b>(0, 0);
}

inline bool dissimilar_colors(const uint8_t yuv1[3], const uint8_t yuv2[3])
{
    // From the hqx algorithm, Stephin 2003
    return abs(yuv1[0] - yuv2[0]) > 48
        || abs(yuv1[1] - yuv2[1]) > 7
        || abs(yuv1[2] - yuv2[2]) > 6;
}

inline bool shading_edge(const uint8_t yuv1[], const uint8_t yuv2[])
{
    // Magic numbers taken from Kopf-Lischinski algorithm
    // Used for resolving ambiguities at t-junctions
    return abs(yuv1[0] - yuv2[0]) <= 100
        && abs(yuv1[1] - yuv2[1]) <= 100
        && abs(yuv1[2] - yuv2[2]) <= 100;
}

inline bool contour_edge(const uint8_t yuv1[], const uint8_t yuv2[])
{
    return !shading_edge(yuv1, yuv2);
}

} /* Depixelize */

#endif /* DEPIXELIZE_COLOR_UTIL_HPP */
