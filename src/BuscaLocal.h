#ifndef BUSCALOCAL_H
#define BUSCALOCAL_H

#define EPSILON 0.0000001

using namespace std;

#include "Parameters.h"
#include "Construction.h"

class BuscaLocal{
    private:
    public:
        Parameters *parameters;
        vector<int> U;
        vector<int> endTimeOperations;
        int makespan;

        BuscaLocal(Parameters *parameters, vector<int> U, vector<int> endTimeOperations, int makespan);
        bool redundancySwapAdjacent(unsigned op1, unsigned op2);
        bool redundancySwap(const vector<int>& U, int index1, int index2);
        bool swap();
        bool relocate();
        bool relocateBlock(int blockSize);
        bool relocateBlock2();
        bool relocateBlock3();
        bool relocateBlock4();
        bool searchNeighborhood(unsigned int i);
        void runSearchTotal();

};

#endif