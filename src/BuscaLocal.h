#ifndef BUSCALOCAL_H
#define BUSCALOCAL_H

#define EPSILON 0.0000001

using namespace std;

#include "Parameters.h"
#include "Construction.h"

class BuscaLocal{
    public:
        bool redundancySwapAdjacent(Parameters *parameters, unsigned op1, unsigned op2);
        bool redundancySwap(Parameters *parameters, const vector<int>& U, int index1, int index2);
        bool swap(Parameters *parameters, vector<int>& U, int makespan);
        bool relocate(Parameters *parameters, vector<int>& U, int makespan);
        bool relocate2(Parameters *parameters, vector<int>& U, int makespan);

};

#endif