#ifndef INDIVIDUO_H
#define INDIVIDUO_H

#include <vector>
#include <list>
#include <iostream>
#include "Parameters.h"
#include "Construction.h"

using namespace std;

typedef struct{

    Individuo *individual;

    double distance;
}Proximity;

class Individuo{
    private:
        Parameters *parameters;
        vector<int> endTimeOperations;
    public:
        vector<int> chromosome;
        unsigned int makespan;
        float fitRank;
        unsigned int calcMakespan();
        Individuo(Parameters *parameters);


    list<Proximity> closest;

    int calculateDistance(Individuo *individuo);
};

#endif