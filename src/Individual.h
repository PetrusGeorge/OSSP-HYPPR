#ifndef Individual_H
#define Individual_H

#include <vector>
#include <list>
#include <iostream>
#include "Parameters.h"
#include "Construction.h"

using namespace std;

class Individual;

typedef struct{

    Individual *individual;
    double distance;

} Proximity;

class Individual{
    private:
       
    public:
        Parameters *parameters;
        Individual(Parameters *parameters);
        Individual(Parameters *parameters, int c);
        Individual();

        vector<int> chromosome;
        unsigned int makespan;
        vector<int> endTimeOperations;

        unsigned age;

        float fitRank;
        float fitnessExt;
        float divRank;

        list<Proximity> closest;

        unsigned int calcMakespan();

        int calculateDistance(Individual *Individual);
        void addClose(Individual *indiv);
        void removeClose(Individual *indiv);
        double distToClosests(int n);

        void recopyIndividual(Individual *indiv);
        bool verifySequence();
};

#endif