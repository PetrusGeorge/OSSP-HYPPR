#ifndef Population_H
#define Population_H

#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include "Individual.h"
#include "LocalSearch.h"

using namespace std;

/// Structure to manage the population
struct SubPopulation {
    vector<Individual*> individuals; ///< individuals in the population
    int numberIndividuals; ///< number of individuals in the population
};

class Population {

private:
    
    // busca_local
    // void education(Individual *indiv);

    // insere Individual na população
    int placeIndividual(SubPopulation *subPop, Individual *indiv);
    void education(Individual *indiv);

    Parameters *parameters;
    LocalSearch *BL;
public:
    // Individual *trainer; //< structure to instantiate only one individual for local search procedures
    SubPopulation *subPopulation; ///< structure to handle the population
    Population(Parameters *parameters, LocalSearch *BS);
    ~Population();
    int addIndividual(Individual *indiv);
    void updateProximity(SubPopulation *subPop, Individual *indiv);
    void removeIndividual(SubPopulation *subPop, int p);

    int selectToRemove(SubPopulation *subPop);
    void updateAge();
    void evalExtFit(SubPopulation *subPop);

    Individual* getIndividualBinT();
    Individual* getBestIndividual();
    void diversify();

};

#endif