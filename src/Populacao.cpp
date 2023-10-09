#include "Populacao.h"

Populacao::Populacao(Parameters *parameters){

    this->parameters = parameters;
    Individuo *randomIndiv;
    subPopulation = new SubPopulacao();
    subPopulation->numberIndividuals = 0;

    // trainer = new Individual(parameters);
    // trainer->localSearch = new LocalSearch(parameters, trainer);


    // Create the initial population
    for (unsigned int i = 0; i < parameters->populationSize; i++) {
        randomIndiv = new Individuo(parameters);
        // education(randomIndiv);
        addIndividuo(randomIndiv);
        delete randomIndiv;
    }
}

Populacao::~Populacao() {
    int size;
    if (subPopulation != nullptr) {
        size = (int)subPopulation->individuals.size();
        for (int i = 0; i < size; i++)
            delete subPopulation->individuals[i];
        delete subPopulation;
    }
    // delete trainer->localSearch;
    // delete trainer;
}

void Populacao::updateProximity(SubPopulacao *subPop, Individuo *indiv) {
    for (int k = 0; k < subPop->numberIndividuals; k++) {
        if (subPop->individuals[k] != indiv) {
            subPop->individuals[k]->addClose(indiv);
            indiv->addClose(subPop->individuals[k]);
        }
    }
}

int Populacao::placeIndividual(SubPopulacao *subPop, Individuo *indiv) {
    Individuo *myIndiv = new Individuo(parameters);


    myIndiv->recopyIndividual(indiv);

    int i = (int)subPop->individuals.size() - 1;
    subPop->individuals.push_back(myIndiv);

    while (i >= 0) {
        if (subPop->individuals[i]->makespan >= indiv->makespan) {
            subPop->individuals[i + 1] = subPop->individuals[i];
            i--;
        }
        else {
            subPop->individuals[i + 1] = myIndiv;
            subPop->numberIndividuals++;
            updateProximity(subPop, subPop->individuals[i + 1]);

            return i + 1; // success
        }
    }

    subPop->individuals[0] = myIndiv;
    subPop->numberIndividuals++;
    updateProximity(subPop, subPop->individuals[0]);

    return 0;
}

int Populacao::addIndividuo(Individuo *indiv) {
    SubPopulacao *subPop;
    int k, result;
    bool firstIt = true;

    subPop = subPopulation;
    result = placeIndividual(subPop, indiv);

    // Keep only the survivors if the maximum size of the population has been reached
    // if (result != -1 && subPop->numberIndividuals > parameters->populationSize + parameters->maxPopulationSize) {
    //     while (subPop->numberIndividuals > parameters->populationSize) {
    //         k = selectToRemove(subPop);
    //         removeIndividual(subPop, k);
    //         if (firstIt) {
    //             firstIt = false;
    //         }
    //     }
    // }
    return result;
}

void Populacao::removeIndividuo(SubPopulacao *subPop, int p) {
    Individuo *remIndiv = subPop->individuals[p];

    // Place individual at the end
    for (int i = p + 1; i < (int)subPop->individuals.size(); i++)
        subPop->individuals[i - 1] = subPop->individuals[i];

    // Remove it from the population
    subPop->individuals.pop_back();
    subPop->numberIndividuals--;

    // Remove it from the proximity structures
    for (int i = 0; i < subPop->numberIndividuals; i++)
        subPop->individuals[i]->removeClose(remIndiv);

    delete remIndiv;
}

