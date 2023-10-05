#ifndef GENETICO_H
#define GENETICO_H

#include "Parameters.h"
#include "Construction.h"
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

    public:
        int crossoverOX(const vector<int>& parent1,const vector<int>& parent2, vector<int>& child);
        Genetico(Parameters *parameters);
        void evolve();

};

#endif //GENETICO_H