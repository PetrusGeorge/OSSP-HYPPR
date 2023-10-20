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
       
    public:
        Parameters *parameters;
        Individuo(Parameters *parameters);
        Individuo(Parameters *parameters, int c);
        Individuo();

        vector<int> chromosome;
        unsigned int makespan;
        vector<int> endTimeOperations;

        unsigned age;

        float fitRank;
        float fitnessExt;
        float divRank;

        list<Proximity> closest;

        unsigned int calcMakespan();

        int calculateDistance(Individuo *individuo);
        void addClose(Individuo *indiv);
        void removeClose(Individuo *indiv);
        double distToClosests(int n);

        void recopyIndividual(Individuo *indiv);
        bool verifySequence();
};

#endif