#include "Genetico.h"

Genetico::Genetico(Parameters *parameters) {
    this->parameters = parameters;
}

void Genetico::evolve(){

}

int Genetico::crossoverOX(const vector<int>& parent1,const vector<int>& parent2, vector<int>& child) {
    
    // Beginning and end of the crossover zone
    unsigned int begin = rand() % parameters->numJobs;
    unsigned int end = rand() % parameters->numJobs;

    while (end == begin && parameters->numJobs > 1)
        end = rand() % parameters->numJobs;

    if (begin > end) {
        unsigned int temp = begin;
        begin = end;
        end = temp;
    }

    // Copy part of parent1 to child
    child = parent1;

    for (unsigned int i = 0; i < parameters->positionsOffspring.size(); i++) {
        parameters->positionsOffspring[i] = false;
    }
    for (unsigned int i = begin; i <= end; i++) {
        parameters->positionsOffspring[parent1[i]] = true;
    }
   
    // Copy unused values of parent2 to child sequentially
    unsigned int pos = end + 1, i = end + 1;

    while (pos < parameters->positionsOffspring.size()) {
        if (!parameters->positionsOffspring[parent2[i]]) {
            child[pos] = parent2[i];
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
        if (!parameters->positionsOffspring[parent2[i]]) {
            child[pos] = parent2[i];
            pos++;
        }
        i++;
        if (i == parameters->positionsOffspring.size()) {
            i = 0;
        }
    }

    return calculateMakespan(parameters, child);
    
}

