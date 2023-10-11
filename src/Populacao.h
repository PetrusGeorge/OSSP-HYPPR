#ifndef POPULACAO_H
#define POPULACAO_H

#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include "Individuo.h"

using namespace std;

/// Structure to manage the population
struct SubPopulacao {
    vector<Individuo*> individuals; ///< individuals in the population
    int numberIndividuals; ///< number of individuals in the population
};

class Populacao {

private:
    
    // busca_local
    // void education(Individual *indiv);

    // insere individuo na população
    int placeIndividual(SubPopulacao *subPop, Individuo *indiv);

    Parameters *parameters;
public:
    // Individual *trainer; //< structure to instantiate only one individual for local search procedures
    SubPopulacao *subPopulation; ///< structure to handle the population
    Populacao(Parameters *parameters);
    ~Populacao();
    int addIndividual(Individuo *indiv);
    void updateProximity(SubPopulacao *subPop, Individuo *indiv);
    void removeIndividual(SubPopulacao *subPop, int p);

    int selectToRemove(SubPopulacao *subPop);
    void updateAge();
    void evalExtFit(SubPopulacao *subPop);
    bool fitExist(SubPopulacao *subPop, Individuo *indiv);

    Individuo* getIndividualBinT();
    void diversify();

};

#endif