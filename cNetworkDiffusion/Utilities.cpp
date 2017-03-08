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
        a0 = accumulate(weights.begin(), weights.end(), 0.0);

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


size_t get_random_neighbor_MHRN(
        size_t const & u,
        size_t const & B,
        vector < double > const & p_node_in_layers,
        default_random_engine & generator, 
        uniform_real_distribution<double> & distribution
     )
{
    //size_t L = p_node_in_layers.size();
    /*
    vector < double > p_node_in_layers;
    
    for(size_t l = 1; l <= L; l++)
    {
        size_t nodes_in_layer = pow(B,l)-pow(B,l-1);
        p_node_in_layers.push_back( p_l[l-1] * nodes_in_layer );
    }
    */

    size_t neighbor_layer = arg_choose_from_vector(
                                p_node_in_layers,
                                generator,
                                distribution
                               );

    //because we start counting the layers at 1
    neighbor_layer++;

    size_t nodes_in_layer_module = pow(B,neighbor_layer);
    size_t layer_module_id = u / nodes_in_layer_module;
    size_t neigh_min = layer_module_id * nodes_in_layer_module;
    size_t neigh_max = (layer_module_id+1) * nodes_in_layer_module - 1;

    size_t nodes_in_lower_layer_module = pow(B,neighbor_layer-1);
    size_t lower_layer_module_id = u / nodes_in_lower_layer_module;
    size_t forbidden_min = lower_layer_module_id * nodes_in_lower_layer_module;
    size_t forbidden_max = (lower_layer_module_id+1) * nodes_in_lower_layer_module - 1;

    uniform_int_distribution<size_t> random_integer(neigh_min,neigh_max); // inclusive interval
    
    size_t v = random_integer(generator);

    // draw from this layer as long as the drawn neighbor v is in the same module
    // as the orginal node one layer lower
    while ( (v>=forbidden_min) and (v<=forbidden_max) )
        v = random_integer(generator);

    return v;
}

vector < double > get_p_MHRN(
        size_t const & B,
        size_t const & L,
        double const & k,
        double const & xi
        )
{
    double p1;
    if ( xi == 1.0 ) 
        p1 = k / double(B-1) / double(L);
    else
        p1 = k / double(B-1) * (1.0-xi) / (1.0-pow(xi,L));

    if (p1>1.0)
        throw domain_error("The lowest layer connection probability is >1.0, meaning that either xi is too small or k is too large.");

    vector < double > p(L);
    p[0] = p1;
    for (size_t l=2; l<=L; l++)
        p[l-1] = p1 * pow(xi/double(B),l-1);

    return p;
}
