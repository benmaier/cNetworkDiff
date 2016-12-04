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

#include "Utilities.h"
#include "diffusion.h"

#include <iostream>
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
#include <assert.h>

using namespace std;

tuple < double, double > mmfpt_and_mean_coverage_time(
        size_t N,
        vector < pair <size_t, size_t > > edge_list,
        double coverage_ratio,
        size_t seed
        )
{
    assert( coverage_ratio>0 && coverage_ratio<=1.0);

    vector < vector < size_t > * > G = get_neighbor_list(N,edge_list);

    double mfpt = 0.; // mean first passage time
    double covt = 0.; // coverage time

    size_t nmax = coverage_ratio * N;

    size_t N_w = N;
    vector < size_t > current_nodes(N_w);
    vector < set < size_t > * > already_visited;

    vector < size_t > remaining_walkers;

    for(size_t node = 0; node < N; node++)
    {
        size_t walker = node;
        current_nodes[walker] = node;
        already_visited.push_back(new set < size_t >);
        already_visited[walker]->insert(node);
        remaining_walkers.push_back(walker);
    }

    //initialize random generators
    default_random_engine generator(seed);
    uniform_real_distribution<double> uni_distribution(0.,1.);
    
    size_t t = 1;

    while (remaining_walkers.size()>0)
    {
        vector < size_t > next_to_pop;
        for(auto const& walker: remaining_walkers)
        {
            size_t u = current_nodes[walker];
            size_t k = G[u]->size();
            size_t neigh = G[u]->at(uni_distribution(generator) * k);

            current_nodes[walker] = neigh;

            if (already_visited[walker]->find(neigh) == already_visited[walker]->end())
            {
                already_visited[walker]->insert(neigh);
                mfpt += t;
            }

            if (already_visited[walker]->size()==nmax)
            {
                covt += t;
                next_to_pop.push_back(walker);
            }

        }

        vector<size_t>::iterator n_it;
        for (n_it = next_to_pop.begin(); n_it != next_to_pop.end(); ++n_it)
        {
            vector<size_t>::iterator to_delete = find(remaining_walkers.begin(),remaining_walkers.end(),*n_it);
            *to_delete = remaining_walkers.back();
            remaining_walkers.pop_back();
        }

        t++;
    }

    //free memory
    for(size_t node = 0; node < N; node++)
    {
        delete already_visited[node];
        delete G[node];
    }

    covt /= N*(N-1);
    mfpt /= N*(N-1);

    return make_pair(mfpt, covt);

}
