#ifndef Genetic_H
#define Genetic_H

#include "Parameters.h"
#include "Construction.h"
#include "Population.h"
#include "LocalSearch.h"
#include "time.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <list>
#include <math.h>

using namespace std;

class Genetic
{
    private:
        Parameters *parameters;

        LocalSearch *BL;

        int nbIterWithoutImprov; ///< number of iterations without improvement

        int nbIter; ///< number of iterations

    public:
        Population *population;

        Individual* crossoverOX(Individual *parent1, Individual *parent2);
        Genetic(Parameters *parameters, Population *population, LocalSearch *BL);
        void evolve(int maxIterWithoutImprov);
        Individual* RR(Individual * parent1);
        Individual* swapGenetic(Individual * parent1);

};

#endif //Genetic_H