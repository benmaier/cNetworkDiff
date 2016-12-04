#ifndef __DIFFUSION_H__
#define __DIFFUSION_H__

#include "Utilities.h"
#include <random>
#include <tuple>

using namespace std;

tuple < double, double > mmfpt_and_mean_coverage_time(
        size_t N,
        vector < pair <size_t, size_t > > edge_list,
        double coverage_ratio = 1.0,
        size_t seed = 0
        );

#endif
