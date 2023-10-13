#include "Genetico.h"

Genetico::Genetico(Parameters *parameters, Populacao *population, BuscaLocal *BL) {
    this->parameters = parameters;
    this->population = population;
    this->BL = BL;
}

void Genetico::evolve(int maxIterWithoutImprov){
    // Individuals used for crossover
    Individuo *parent1;
    Individuo *parent2;

    nbIterWithoutImprov = 1;
    int nbIterWithoutImprovDiv = 1;
    nbIter = 1;

    string temp;
    int place;
    clock_t debut = clock();

    // Child reference
    Individuo *offspring = new Individuo(parameters);

    while (nbIterWithoutImprov < maxIterWithoutImprov) {
        //cout << "makespan: " << population->getBestIndividual()->makespan << " nbIterWithoutImprov: "<< nbIterWithoutImprov << endl;
        // CROSSOVER
        parent1 = population->getIndividualBinT(); // Pick individual by binary tournament
        parent2 = population->getIndividualBinT(); // Pick individual by binary tournament

        offspring = crossoverOX(parent1, parent2); // OX crossover
        cout << offspring->makespan << endl;

        // LOCAL SEARCH
        for(int i = 0; i < 10; i++){
            BL->runSearchTotal(offspring);
        }
        

        cout << offspring->makespan << endl;

        // Tries to add child to population
        place = population->addIndividual(offspring);
        if (place == -2) {
            return;
        }

        if (place == 0) { // A new best solution has been found
            nbIterWithoutImprov = 1;
            nbIterWithoutImprovDiv = 1;
        }
        else
            nbIterWithoutImprov++;

        nbIterWithoutImprovDiv++;
        nbIter++;

        // DIVERSIFICATION
        // Max iterations without improvement resulting in diversification reached
        if (nbIterWithoutImprovDiv == parameters->maxDiversify) {
            population->diversify();
            if (parameters->terminate) {
                return;
            }
            nbIterWithoutImprovDiv = 1;
        }
    }
    parameters->nbIter = (unsigned int) nbIter;
}

Individuo* Genetico::crossoverOX(Individuo *parent1, Individuo *parent2) {
    
    // Beginning and end of the crossover zone
    unsigned int begin = rand() % parameters->numJobs * parameters->numJobs;
    unsigned int end = rand() % parameters->numJobs * parameters->numJobs;

    Individuo *child = new Individuo(parameters, 0);
    child->recopyIndividual(parent1);

    while (end == begin && parameters->numJobs > 1)
        end = rand() % parameters->numJobs;

    if (begin > end) {
        unsigned int temp = begin;
        begin = end;
        end = temp;
    }

    for (unsigned int i = 0; i < parameters->positionsOffspring.size(); i++) {
        parameters->positionsOffspring[i] = false;
    }
    for (unsigned int i = begin; i <= end; i++) {
        parameters->positionsOffspring[parent1->chromosome[i]-1] = true;
    }
   
    // Copy unused values of parent2 to child sequentially
    //2 4 7 5 3 6 11 12 1 15 13 10 8 14 9 8 
    unsigned int pos = end + 1, i = end + 1;

    while (pos < parameters->positionsOffspring.size()) {
        if (!parameters->positionsOffspring[parent2->chromosome[i]-1]) {
            child->chromosome[pos] = parent2->chromosome[i];
            parameters->positionsOffspring[parent2->chromosome[i]-1] = true;
            pos++;
        }
        i++;
        if (i == parameters->positionsOffspring.size()) {
            i = 0;
        }
    }
    if (i == parameters->positionsOffspring.size()) {
        i = 0;
    }
    pos = 0;
    while (pos < begin) {
        if (!parameters->positionsOffspring[parent2->chromosome[i]-1]) {
            child->chromosome[pos] = parent2->chromosome[i];
            pos++;
        }
        i++;
        if (i == parameters->positionsOffspring.size()) {
            i = 0;
        }

    }

    child->makespan = calculateMakespan(parameters, child->chromosome, child->endTimeOperations);

    return child;
    
}

