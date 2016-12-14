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


using namespace std;

size_t arg_choose_from_vector(
        vector < double > const & weights, 
        default_random_engine & generator, 
        uniform_real_distribution<double> & distribution,
        double a0
      )
{
    if (a0 == 0.0)
        double a0 = accumulate(weights.begin(), weights.end(), 0.0);

    double rProduct = distribution(generator) * a0;

    int event = 0;
    int N = weights.size();

    if (N==0)
    {
        throw length_error( "Rate list is empty." );
    }

    double sum_event = 0.0;
    while ( (event<N) and not ( (sum_event < rProduct) and (rProduct <= sum_event+weights[event]) ) )
    {
        sum_event += weights[event];
        event++;
    }

    return event;
}


vector < set < size_t > * > get_neighbor_set(
        size_t N,
        vector < pair < size_t, size_t > > &edge_list
        )
{
    vector < set < size_t > * > G;

    for(size_t node = 0; node< N; node++)
        G.push_back(new set <size_t>);

    for(auto edge: edge_list)
    {
        G[edge.first]->insert(edge.second);
        G[edge.second]->insert(edge.first);
    }

    return G;
}

vector < vector < size_t > * > get_neighbor_list(
        size_t N,
        vector < pair < size_t, size_t > > &edge_list
        )
{
    vector < vector < size_t > * > G;

    for(size_t node = 0; node< N; node++)
        G.push_back(new vector <size_t>);

    for(auto const& edge: edge_list)
    {
        G[edge.first]->push_back(edge.second);
        G[edge.second]->push_back(edge.first);
    }

    return G;
}

void rm_element(
        vector < size_t > &vec,
        size_t idx
        )
{
    vec[idx] = vec.back();
    vec.pop_back();
}

void add_nodes_belonging_to_this_component(
        size_t start_node,
        const vector < set < size_t > * > &G,
        set < size_t > * comp,
        vector < bool > &already_visited
       )
{
    comp->insert(start_node);
    already_visited[start_node] = true;
    for(auto const& neigh: *G[start_node])
    {
        if ( not already_visited[neigh] )
        {
            add_nodes_belonging_to_this_component(neigh,G,comp,already_visited);
        }
    }
}

// returns a list of sets, each set contains the nodes belonging to the
// component
vector < set < size_t > * > get_components(
        const vector < set < size_t > * > &G
        )
{
    vector < bool > already_visited;
    for(size_t node=0; node<G.size(); node++)
        already_visited.push_back(false);

    vector < set <size_t> * > components;

    for(size_t node=0; node<G.size(); node++)
    {
        if ( not already_visited[node] )
        {
            components.push_back(new set <size_t>);
            add_nodes_belonging_to_this_component(node,G,components.back(),already_visited);
        }
    }

    return components;
}

vector < set < size_t > * > get_components_from_edgelist(
        size_t N,
        vector < pair < size_t,size_t > > &edge_list
        )
{
    vector < set < size_t > * > G;
    for(size_t node=0; node<N; node++)
        G.push_back( new set <size_t> );

    for(auto edge: edge_list)
    {
        G[edge.first]->insert(edge.second);
        G[edge.second]->insert(edge.first);
    }

    return get_components(G);
}

// replace the graph in-place with the giant component
void get_giant_component(
        vector < set < size_t > * > &G
        )
{
    vector < set < size_t > * > components = get_components(G);
    size_t max = components[0]->size();
    size_t max_comp = 0;
    for(size_t comp = 1; comp < components.size(); comp++)
        if (components[comp]->size()>max)
        {
            max = components[comp]->size();
            max_comp = comp;
        }

    for(size_t node = 0; node<G.size(); node++)
    {
        if ( components[max_comp]->find(node) == components[max_comp]->end() )
        {
            // if current node not in giant component, empty the neighbor set
            G[node]->clear();
        }
    }

    //free memory
    for(size_t comp = 0; comp<components.size(); comp++)
        delete components[comp];
}



