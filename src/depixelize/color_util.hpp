#ifndef DEPIXELIZE_COLOR_UTIL_HPP
#define DEPIXELIZE_COLOR_UTIL_HPP

namespace Depixelize {

inline bool dissimilar_colors(const unsigned char yuv1[3], const unsigned char yuv2[3])
{
    // From the hqx algorithm, Stephin 2003
    return abs(yuv1[0] - yuv2[0]) > 48
        || abs(yuv1[1] - yuv2[1]) > 7
        || abs(yuv1[2] - yuv2[2]) > 6;
}


} /* Depixelize */

#endif /* DEPIXELIZE_COLOR_UTIL_HPP */
