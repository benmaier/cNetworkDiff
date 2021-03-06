#ifndef __DIFFUSION_H__
#define __DIFFUSION_H__

#include "Utilities.h"
#include <random>
#include <tuple>

using namespace std;

size_t saturation_time(
        size_t N,
        vector < pair <size_t, size_t > > edge_list,
        size_t source,
        size_t target,
        size_t N_walker,
        size_t N_saturation,
        size_t seed
        );

tuple < vector < double >, double > gmfpts_and_mean_coverage_time(
        size_t N,
        vector < pair <size_t, size_t > > edge_list,
        double coverage_ratio = 1.0,
        size_t seed = 0
        );

tuple < double, double > mgmfpt_and_mean_coverage_time(
        size_t N,
        vector < pair <size_t, size_t > > edge_list,
        double coverage_ratio = 1.0,
        size_t seed = 0
        );

tuple < double, double > mmfpt_and_mean_coverage_time_power_law(
        size_t N,
        double alpha,
        double coverage_ratio = 1.0,
        size_t seed = 0
        );

tuple < double, double > mmfpt_and_mean_coverage_time_power_law_k(
        size_t N,
        double alpha,
        double k,
        double coverage_ratio,
        size_t seed
        );

tuple < double, double > mmfpt_and_mean_coverage_time_meanfield_MHRN(
        size_t B,
        size_t L,
        double xi,
        double coverage_ratio,
        size_t seed
        );

vector < double > cover_times(
        size_t N,
        vector < pair <size_t, size_t > > edge_list,
        double coverage_ratio,
        size_t seed
        );
#endif
