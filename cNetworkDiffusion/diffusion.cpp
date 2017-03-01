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
    size_t number_of_pairs = 0;

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
                number_of_pairs++;
            }

            if (already_visited[walker]->size()==nmax)
            {
                covt += t;
                next_to_pop.push_back(walker);
            }

        }

        for (size_t walker_id=0; walker_id<next_to_pop.size(); walker_id++)
        {
            vector<size_t>::iterator to_delete = find(remaining_walkers.begin(),
                                                      remaining_walkers.end(),
                                                      next_to_pop[walker_id]
                                                     );
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

    covt /= N;
    mfpt /= number_of_pairs;

    return make_pair(mfpt, covt);

}

tuple < double, double > mmfpt_and_mean_coverage_time_meanfield_MHRN(
        size_t B,
        size_t L,
        double k,
        double xi,
        double coverage_ratio,
        size_t seed
        )
{
    assert( coverage_ratio>0 && coverage_ratio<=1.0);

    size_t N = pow(B,L);

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

    vector < double > layer_pmf = get_p_MHRN(B,L,k,xi);
    for(size_t l = 1; l<=L; l++)
        layer_pmf[l-1] *= ( pow(B,l) - pow(B,l-1) );


    
    size_t t = 1;
    size_t number_of_pairs = 0;

    while (remaining_walkers.size()>0)
    {
        vector < size_t > next_to_pop;
        for(auto const& walker: remaining_walkers)
        {
            size_t u = current_nodes[walker];
            size_t neigh = get_random_neighbor_MHRN(u,B,layer_pmf,generator,uni_distribution);

            current_nodes[walker] = neigh;

            if (already_visited[walker]->find(neigh) == already_visited[walker]->end())
            {
                already_visited[walker]->insert(neigh);
                mfpt += t;
                number_of_pairs++;
            }

            if (already_visited[walker]->size()==nmax)
            {
                covt += t;
                next_to_pop.push_back(walker);
            }

        }

        for (size_t walker_id=0; walker_id<next_to_pop.size(); walker_id++)
        {
            vector<size_t>::iterator to_delete = find(remaining_walkers.begin(),
                                                      remaining_walkers.end(),
                                                      next_to_pop[walker_id]
                                                     );
            *to_delete = remaining_walkers.back();
            remaining_walkers.pop_back();
        }

        t++;
    }

    //free memory
    for(size_t node = 0; node < N; node++)
    {
        delete already_visited[node];
    }

    covt /= N;
    mfpt /= number_of_pairs;

    return make_pair(mfpt, covt);

}
tuple < double, double > mmfpt_and_mean_coverage_time_power_law(
        size_t N,
        double alpha,
        double coverage_ratio,
        size_t seed
        )
{
    assert( coverage_ratio>0 && coverage_ratio<=1.0);


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

    size_t N_half = N / 2;
    vector < double > pmf(N-1);
    double a0 = 0.0;

    for(size_t neigh = 0; neigh< N_half; neigh++)
    {
        pmf[neigh] = pow(neigh+1.0,alpha);
        pmf[N-2-neigh] = pmf[neigh];
        if (neigh == N-2-neigh)
            a0 += pmf[neigh];
        else
            a0 += 2*pmf[neigh];
    }


    
    size_t t = 1;
    size_t number_of_pairs = 0;

    while (remaining_walkers.size()>0)
    {
        vector < size_t > next_to_pop;
        for(auto const& walker: remaining_walkers)
        {
            size_t u = current_nodes[walker];
            size_t neigh = arg_choose_from_vector(pmf,generator,uni_distribution,a0) + 1;
            neigh = (u + neigh) % N;

            current_nodes[walker] = neigh;

            if (already_visited[walker]->find(neigh) == already_visited[walker]->end())
            {
                already_visited[walker]->insert(neigh);
                mfpt += t;
                number_of_pairs++;
            }

            if (already_visited[walker]->size()==nmax)
            {
                covt += t;
                next_to_pop.push_back(walker);
            }

        }

        for (size_t walker_id=0; walker_id<next_to_pop.size(); walker_id++)
        {
            vector<size_t>::iterator to_delete = find(remaining_walkers.begin(),
                                                      remaining_walkers.end(),
                                                      next_to_pop[walker_id]
                                                     );
            *to_delete = remaining_walkers.back();
            remaining_walkers.pop_back();
        }

        t++;
    }

    //free memory
    for(size_t node = 0; node < N; node++)
    {
        delete already_visited[node];
    }

    covt /= N;
    mfpt /= number_of_pairs;

    return make_pair(mfpt, covt);

}

tuple < double, double > mmfpt_and_mean_coverage_time_power_law_k(
        size_t N,
        double alpha,
        double k,
        double coverage_ratio,
        size_t seed
        )
{
    assert( coverage_ratio>0 && coverage_ratio<=1.0);


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

    size_t N_half = N / 2;
    vector < double > pmf(N-1);
    double a0 = 0.0;

    for(size_t neigh = 0; neigh< N_half; neigh++)
    {
        pmf[neigh] = pow(neigh+1.0,alpha);
        pmf[N-2-neigh] = pmf[neigh];
        if (neigh == N-2-neigh)
            a0 += pmf[neigh];
        else
            a0 += 2*pmf[neigh];
    }

    double overflow_probability = 0.0;

    for(size_t neigh = 0; neigh< N_half; neigh++)
    {
        pmf[neigh] *= k/a0;
        pmf[N-2-neigh] = pmf[neigh];

        if (pmf[neigh]>=1.)
        {
            overflow_probability += pmf[neigh] - 1.0;
            pmf[neigh] = 1.0;
            pmf[N-2-neigh] = 1.0;
        }
        else if (overflow_probability>0.)
        {
            double dp = 1.0 - pmf[neigh];
            if (dp>overflow_probability)
            {
                pmf[neigh] += overflow_probability;
                if (neigh != N-2-neigh)
                    pmf[N-2-neigh] += overflow_probability;
                else if (1.0-pmf[neigh] > overflow_probability)
                    pmf[neigh] += overflow_probability;
                else
                    pmf[neigh] = 1.;

                overflow_probability = 0.0;
            }
            else
            {
                pmf[neigh] = 1.0;
                pmf[N-2-neigh] = 1.0;
                overflow_probability -= dp;
            }
        }
    }

    a0 = k;


    
    size_t t = 1;
    size_t number_of_pairs = 0;

    while (remaining_walkers.size()>0)
    {
        vector < size_t > next_to_pop;
        for(auto const& walker: remaining_walkers)
        {
            size_t u = current_nodes[walker];
            size_t neigh = arg_choose_from_vector(pmf,generator,uni_distribution,a0) + 1;
            neigh = (u + neigh) % N;

            current_nodes[walker] = neigh;

            if (already_visited[walker]->find(neigh) == already_visited[walker]->end())
            {
                already_visited[walker]->insert(neigh);
                mfpt += t;
                number_of_pairs++;
            }

            if (already_visited[walker]->size()==nmax)
            {
                covt += t;
                next_to_pop.push_back(walker);
            }

        }

        for (size_t walker_id=0; walker_id<next_to_pop.size(); walker_id++)
        {
            vector<size_t>::iterator to_delete = find(remaining_walkers.begin(),
                                                      remaining_walkers.end(),
                                                      next_to_pop[walker_id]
                                                     );
            *to_delete = remaining_walkers.back();
            remaining_walkers.pop_back();
        }

        t++;
    }

    //free memory
    for(size_t node = 0; node < N; node++)
    {
        delete already_visited[node];
    }

    covt /= N;
    mfpt /= number_of_pairs;

    return make_pair(mfpt, covt);

}
