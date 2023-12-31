#include "Population.h"

Population::Population(Parameters *parameters, LocalSearch *BL){

    this->parameters = parameters;
    this->BL = BL;

    Individual *randomIndiv;
    subPopulation = new SubPopulation();
    subPopulation->numberIndividuals = 0;

    // Create the initial population
    for (unsigned int i = 0; i < parameters->populationSize; i++) {
        randomIndiv = new Individual(parameters);
        //education(randomIndiv);
        addIndividual(randomIndiv);
        delete randomIndiv;
    }
}

void Population::education(Individual *indiv) {
    Individual *trainer = new Individual(parameters, 0);
    trainer->recopyIndividual(indiv);
    BL->runSearchTotal(indiv);
    indiv->recopyIndividual(trainer);
    delete trainer;
}

Population::~Population() {
    int size;
    if (subPopulation != nullptr) {
        size = (int)subPopulation->individuals.size();
        for (int i = 0; i < size; i++)
            delete subPopulation->individuals[i];
        delete subPopulation;
    }
}

void Population::updateProximity(SubPopulation *subPop, Individual *indiv) {
    for (int k = 0; k < subPop->numberIndividuals; k++) {
        if (subPop->individuals[k] != indiv) {
            subPop->individuals[k]->addClose(indiv);
            indiv->addClose(subPop->individuals[k]);
        }
    }
}

int Population::placeIndividual(SubPopulation *subPop, Individual *indiv) {
    Individual *myIndiv = new Individual(parameters, 0);

    myIndiv->recopyIndividual(indiv);

    int i = (int)subPop->individuals.size() - 1;
    subPop->individuals.push_back(myIndiv);

    while (i >= 0) {
        if(subPop->individuals[i]->makespan == indiv->makespan){
            subPop->individuals.erase(subPop->individuals.begin()+i+1);
            return -1;
        }
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

    // if(subPop->numberIndividuals >= parameters->populationSize && indiv->makespan == subPop->individuals[0]->makespan){
    //     subPop->individuals.erase(subPop->individuals.begin());
    //     return -1;
    // }

    subPop->individuals[0] = myIndiv;
    subPop->numberIndividuals++;
    updateProximity(subPop, subPop->individuals[0]);

    return 0;
}

int Population::addIndividual(Individual *indiv) {
    SubPopulation *subPop;
    int k, result;
    bool firstIt = true;

    subPop = subPopulation;
    // cout << "pré-place" << endl;
    result = placeIndividual(subPop, indiv);
    // cout << "pos-place" << endl;

    // Keep only the survivors if the maximum size of the population has been reached
    if (result != -1 && subPop->numberIndividuals > parameters->populationSize + parameters->maxPopulationSize) {
        while (subPop->numberIndividuals > parameters->populationSize) {
            k = selectToRemove(subPop);
            removeIndividual(subPop, k);
            if (firstIt) {
                firstIt = false;
            }
        }
    }
    return result;
}

void Population::removeIndividual(SubPopulation *subPop, int p) {
    Individual *remIndiv = subPop->individuals[p];

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

int Population::selectToRemove(SubPopulation *subPop) {
    // Select one individual to be eliminated from the population
    vector<int> position;
    int temp, toRemove;

    updateAge();
    evalExtFit(subPop);

    for (int i = 0; i < subPop->numberIndividuals; i++)
        position.push_back(i);

    // Add a penalty in case of clone
    for (int i = 1; i < subPop->numberIndividuals; i++) {
        if (subPop->individuals[i]->distToClosests(1) <= 0.001) { // in solution space
            subPop->individuals[i]->fitnessExt += 5;
        }
        // if (fitExist(subPop, subPop->individuals[i])) { // in objective space
        //     subPop->individuals[i]->fitnessExt += 5;
        // }
    }

    // Rank elements by extended fitness and select out the worst
    for (int n = 0; n < subPop->numberIndividuals; n++) {
        for (int i = 0; i < subPop->numberIndividuals - n - 1; i++) {
            if (subPop->individuals[position[i]]->fitnessExt > subPop->individuals[position[i + 1]]->fitnessExt) {
                temp = position[i + 1];
                position[i + 1] = position[i];
                position[i] = temp;
            }
        }
    }
    toRemove = position[subPop->numberIndividuals - 1];
    
    return toRemove;
}

void Population::updateAge() {
    for (int i = 0; i < subPopulation->numberIndividuals; i++)
        subPopulation->individuals[i]->age++;
}


void Population::evalExtFit(SubPopulation *subPop) {
    int temp;
    vector<int> position;
    vector<double> distances;

    for (int i = 0; i < subPop->numberIndividuals; i++) {
        position.push_back(i);
        distances.push_back(subPop->individuals[i]->distToClosests(parameters->numberCloseIndividuals));
    }

    // Rank individuals in terms of contribution to diversity
    for (int n = 0; n < subPop->numberIndividuals; n++) {
        for (int i = 0; i < subPop->numberIndividuals - n - 1; i++) {
            if (distances[position[i]] < distances[position[i + 1]] - 0.000001) {
                temp = position[i + 1];
                position[i + 1] = position[i];
                position[i] = temp;
            }
        }
    }

    // Compute the biased fitness
    for (int i = 0; i < subPop->numberIndividuals; i++) {
        subPop->individuals[position[i]]->divRank = (float)i / (float)(subPop->numberIndividuals - 1);
        subPop->individuals[position[i]]->fitRank = (float)subPop->individuals[position[i]]->makespan / (float)(subPop->numberIndividuals - 1);
        subPop->individuals[position[i]]->fitnessExt = subPop->individuals[position[i]]->fitRank + ((float) 1.0 - (float)parameters->numberElite / (float)subPop->numberIndividuals) * subPop->individuals[position[i]]->divRank;
    }
}

Individual* Population::getIndividualBinT() {
    Individual *individual1;
    Individual *individual2;
    int place1, place2;

    // Pick the first individual
    place1 = rand() % (subPopulation->numberIndividuals);
    individual1 = subPopulation->individuals[place1];

    // Pick the second individual
    place2 = rand() % (subPopulation->numberIndividuals);
    individual2 = subPopulation->individuals[place2];

    evalExtFit(subPopulation);

    // Keep the best one
    if (individual1->fitnessExt < individual2->fitnessExt){
        return individual1;
    }
    else{
        return individual2;
    }
        
}

Individual *Population::getBestIndividual() {
    if (subPopulation->numberIndividuals != 0)
        return subPopulation->individuals[0];
    else
        return nullptr;
}

void Population::diversify() {
    Individual *randomIndiv;

    // Remove 70% of the population
    while (subPopulation->numberIndividuals > (int)(0.3 * (double)parameters->populationSize)) {
        delete subPopulation->individuals[subPopulation->numberIndividuals - 1];
        subPopulation->individuals.pop_back();
        subPopulation->numberIndividuals--;
    }

    // Create new individuals until minimum population size is reached
    for (unsigned int i = 0; i < parameters->populationSize; i++) {
        randomIndiv = new Individual(parameters);
        //education(randomIndiv);
        addIndividual(randomIndiv);
        delete randomIndiv;
    }
}

// bool Population::fitExist(SubPopulation *subPop, Individual *indiv) {
//     unsigned int count = 0;
//     unsigned int distance = indiv->makespan;

//     for (int i = 0; i < subPop->numberIndividuals; i++) {
//         if (subPop->individuals[i]->makespan == distance) {
//             count++;
//         }
//     }
//     if (count <= 1)
//         return false;
//     else
//         return true;
// }
