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
    int i =1;

    while (nbIterWithoutImprov < maxIterWithoutImprov) {
        cout << i << " " << nbIterWithoutImprov << endl;
        i++;

        parent1 = population->getIndividualBinT(); // Pick individual by binary tournament
        //parent2 = population->getIndividualBinT(); // Pick individual by binary tournament

        //offspring = crossoverOX(parent1, parent2); // OX crossover
        offspring = RR(parent1);
        BL->runSearchTotal(offspring); // Local search
        //if(!offspring->verifySequence()) continue;

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

Individuo* Genetico::RR(Individuo * parent1) {
    Individuo *offspring = new Individuo(parameters, 0);

    offspring->recopyIndividual(parent1);

    //cout << "INICIAL: " << offspring->makespan << endl;

    vector <int> :: iterator it;

    int m = rand() % parameters->numTools;

    vector<int> opsOfMachines;
    vector<int> indexOfMachines;

    //cout << "machine m: " << m+1 << endl;

    for(int i =0; i < offspring->chromosome.size(); i++){
        int op = offspring->chromosome[i];
        int machineOp = decryptJobMachineIndex(op, parameters->numJobs).second;

        if(machineOp == m) {
            opsOfMachines.push_back(op);
            indexOfMachines.push_back(i);
        }
    }

    // show opInseridos
    /*
    for(int i =0; i < opsOfMachines.size(); i++){
        cout << opsOfMachines[i] << " ";
    } cout << endl;

    //show offspring chromosome
    for(int i =0; i < offspring->chromosome.size(); i++){
        cout << offspring->chromosome[i] << " ";
    } cout << endl;

    // show indexMachine

    for(int i =0; i < indexOfMachines.size(); i++){
        cout << indexOfMachines[i] << " ";
    } cout << endl;
    */


    while(indexOfMachines.size() > 0){
        int indexRandIndex = rand() % indexOfMachines.size();
        int randIndex = indexOfMachines[indexRandIndex];

        int indexRandOp = rand() % opsOfMachines.size();
        int randOp = opsOfMachines[indexRandOp];

        //cout << "randIndex: " << randIndex << " randOp: " << randOp << endl;

        indexOfMachines.erase(indexOfMachines.begin() + indexRandIndex);
        opsOfMachines.erase(opsOfMachines.begin() + indexRandOp);

        offspring->chromosome[randIndex] = randOp;
        /*

        for(int i =0; i < offspring->chromosome.size(); i++){
            cout << offspring->chromosome[i] << " ";
        } cout << endl;
        */
    }
    
    //cout <<  offspring->calcMakespan() << endl;

    /*
    offspring->calcMakespan();
    while(itr <= parameters->jobsMovidos) {
        
        int op = opInseridos[itr];
        vector<int> M, J;

        double bestMk = 9999999;
        vector<int> times, bestTimes;
        int p;

        for(int i =0; i < offspring->chromosome.size()-1; i++){

            times = offspring->endTimeOperations;
            offspring->chromosome.insert(offspring->chromosome.begin() + i, op);

            calculateJMbyIndex(parameters, times, i-1 >= 0 ? i-1 : 0, offspring->chromosome, M, J);
            int mk = updateMakespan(parameters, times, i-1 >= 0 ? i-1 : 0, offspring->chromosome, M, J);

            if(mk < bestMk){
                bestMk = mk;
                bestTimes = times;
                p = i;
            }

            offspring->chromosome.erase(offspring->chromosome.begin() + i);
        }

        offspring->chromosome.insert(offspring->chromosome.begin() + p, op);
        offspring->endTimeOperations = bestTimes;

        itr++;
    }

    
*/

    offspring->makespan = calculateMakespan(parameters, offspring->chromosome, offspring->endTimeOperations);
    return offspring;

}