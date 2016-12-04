#include "CastResult.h"
#include "diffusion.h"
#include "math.h"
#include "matrix.h"
#include "mex.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t N, seed;
    double coverage_ratio;
    
    if (nrhs!=4)
    {
        mexPrintf("Got %d input arguments. ", nrhs);
        throw length_error("Invalid number of input arguments. Has to be 4.");
    }

    if (nlhs!=2)
    {
        mexPrintf("Got %d output arguments. ", nrhs);
        throw length_error("Invalid number of output arguments. Has to be 2.");
    }

    vector < pair <size_t,size_t> > edgelist = get_edgelist(prhs[0]);

    read_single_value(prhs[1],N);
    read_single_value(prhs[2],coverage_ratio);
    read_single_value(prhs[3],seed);

    pair < double, double > > result = mmfpt_and_mean_coverage_time(N,edge_list,coverage_ratio,seed); 

    plhs[0] = result.first;
    plhs[1] = result.second;
}
