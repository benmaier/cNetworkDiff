/* 
 * The MIT License (MIT)
 * Copyright (c) 2016, Benjamin Maier
 *
 * Permission is hereby granted, free of charge, to any person 
 * obtaining a copy of this software and associated documentation 
 * files (the "Software"), to deal in the Software without 
 * restriction, including without limitation the rights to use, 
 * copy, modify, merge, publish, distribute, sublicense, and/or 
 * sell copies of the Software, and to permit persons to whom the 
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall 
 * be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-
 * INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN 
 * AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF 
 * OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
 * IN THE SOFTWARE.
 */

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <set>
#include <utility>
#include <random>
#include <cmath>
#include <numeric>
#include <random>
#include <ctime>
#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Utilities.h"
#include "diffusion.h"
//#include "ResultClasses.h"
//#include "test.h"

using namespace std;
namespace py = pybind11;

PYBIND11_MODULE(cNetworkDiffusion, m) {
    m.doc() = "For fast simulations of random walks on networks.";
    
    m.def("saturation_time", &saturation_time, R"pbdoc(Simulates N_walker random walkers starting on a single source node on the network (given as edge list) and returns the time it takes until at least N_saturation of themarrive at the target node. RNG is non-deterministically initialized if seed = 0.)pbdoc",
            py::arg("N"),
            py::arg("edge_list"),
            py::arg("source"),
            py::arg("target"),
            py::arg("N_walker"),
            py::arg("N_saturation"),
            py::arg("seed") = 0
            );

    m.def("mgmfpt_and_mean_coverage_time", &mgmfpt_and_mean_coverage_time, "Simulates N_walker random walks on the network given as edge list and returns the mean global mean first passage time (MGMFPT) and the mean coverage time. RNG is non-deterministically initialized if seed = 0.",
            py::arg("N"),
            py::arg("edge_list"),
            py::arg("coverage_ratio") = 1.0,
            py::arg("seed") = 0
            );

    m.def("mmfpt_and_mean_coverage_time", &mgmfpt_and_mean_coverage_time, "FOR BACKWARDS COMPATIBILITY ONLY. Simulates N_walker random walks on the network given as edge list and returns the mean global mean first passage time (MGMFPT) and the mean coverage time. RNG is non-deterministically initialized if seed = 0.",
            py::arg("N"),
            py::arg("edge_list"),
            py::arg("coverage_ratio") = 1.0,
            py::arg("seed") = 0
            );

    m.def("gmfpts_and_mean_coverage_time", &gmfpts_and_mean_coverage_time, "Simulates N_walker random walks on the network given as edge list and returns the global mean first passage time (MMFPT) per target node and the overall mean coverage time. RNG is non-deterministically initialized if seed = 0.",
            py::arg("N"),
            py::arg("edge_list"),
            py::arg("coverage_ratio") = 1.0,
            py::arg("seed") = 0
            );

    m.def("mmfpts_and_mean_coverage_time", &gmfpts_and_mean_coverage_time, "FOR BACKWARDS COMPATIBILITY ONLY. Simulates N_walker random walks on the network given as edge list and returns the global mean first passage time (MMFPT) per target node and the overall mean coverage time. RNG is non-deterministically initialized if seed = 0.",
            py::arg("N"),
            py::arg("edge_list"),
            py::arg("coverage_ratio") = 1.0,
            py::arg("seed") = 0
            );
    m.def("mmfpt_and_mean_coverage_time_power_law", &mmfpt_and_mean_coverage_time_power_law, "Simulates N_walker random walks on a ring network with power-law jump distribution and returns the mean mean first passage time (MMFPT) and the mean coverage time. RNG is non-deterministically initialized if seed = 0.",
            py::arg("N"),
            py::arg("alpha"),
            py::arg("coverage_ratio") = 1.0,
            py::arg("seed") = 0
            );

    m.def("mmfpt_and_mean_coverage_time_power_law_k", &mmfpt_and_mean_coverage_time_power_law_k, "Simulates N_walker random walks on a ring network with power-law jump distribution where the distribution sums to k. Returns the mean mean first passage time (MMFPT) and the mean coverage time. RNG is non-deterministically initialized if seed = 0.",
            py::arg("N"),
            py::arg("alpha"),
            py::arg("k"),
            py::arg("coverage_ratio") = 1.0,
            py::arg("seed") = 0
            );

    m.def("mmfpt_and_mean_coverage_time_meanfield_MHRN", &mmfpt_and_mean_coverage_time_meanfield_MHRN, "Simulates N_walker random walks on a ring network with MHRN jump distribution where the distribution sums to k. Returns the mean mean first passage time (MMFPT) and the mean coverage time. RNG is non-deterministically initialized if seed = 0.",
            py::arg("B"),
            py::arg("L"),
            py::arg("xi"),
            py::arg("coverage_ratio") = 1.0,
            py::arg("seed") = 0
            );

    m.def("cover_times", &cover_times,
          "Simulates N random walkers, each starting on a different node. Returns a list of the number of time steps each took to visit every node at least once.",
            py::arg("N"),
            py::arg("edge_list"),
            py::arg("coverage_ratio") = 1.0,
            py::arg("seed") = 0
            );

}
