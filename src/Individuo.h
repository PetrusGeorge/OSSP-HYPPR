#ifndef INDIVIDUO_H
#define INDIVIDUO_H

#include <vector>
#include <list>
#include <iostream>
#include "Parameters.h"
#include "Construction.h"

using namespace std;

class Individuo;

typedef struct{

    Individuo *individual;
    double distance;

} Proximity;

class Individuo{
    private:
        Parameters *parameters;
        vector<int> endTimeOperations;
    public:
        Individuo(Parameters *parameters);
        Individuo();

        vector<int> chromosome;
        unsigned int makespan;

        float fitRank;
        list<Proximity> closest;

        unsigned int calcMakespan();

        int calculateDistance(Individuo *individuo);
        void addClose(Individuo *indiv);
        void removeClose(Individuo *indiv);
        double distToClosests(int n);

        void recopyIndividual(Individuo *indiv);
};

#endif