#ifndef GENETICO_H
#define GENETICO_H

#include "Parameters.h"
#include "Construction.h"
#include "Populacao.h"
#include "BuscaLocal.h"
#include "time.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <list>
#include <math.h>

using namespace std;

class Genetico
{
    private:
        Parameters *parameters;

        BuscaLocal *BL;

        int nbIterWithoutImprov; ///< number of iterations without improvement

        int nbIter; ///< number of iterations

    public:
        Populacao *population;

        Individuo* crossoverOX(Individuo *parent1, Individuo *parent2);
        Genetico(Parameters *parameters, Populacao *population, BuscaLocal *BL);
        void evolve(int maxIterWithoutImprov);

};

#endif //GENETICO_H