#include "Genetic.h"

Genetic::Genetic(Parameters *parameters, Population *population, LocalSearch *BL) {
    this->parameters = parameters;
    this->population = population;
    this->BL = BL;
}

Individual * Genetic::swapGenetic(Individual *parent1){
    Individual* offspring = new Individual(parameters, 0);

    offspring->recopyIndividual(parent1);
    vector<int> M, J;
    vector<int> bestSequence, bestTimes;

    for(int i =0; i < parameters->numJobs; i++){
        for(int j =0; j < parameters->numTools; j++){
            bestSequence = offspring->chromosome;
            bestTimes = offspring->endTimeOperations;
            int bestMakespan = offspring->makespan;

            for(int k = j+1; k < parameters->numTools; k++){
                int indexOp = j+i*parameters->numTools;
                int indexOp2 = k+i*parameters->numTools;
                
                if(decryptJobMachineIndex(offspring->chromosome[indexOp], parameters->numJobs).first != decryptJobMachineIndex(offspring->chromosome[indexOp2], parameters->numJobs).first){
                    continue;
                }

                vector<int> sequence = offspring->chromosome;
                vector<int> times = offspring->endTimeOperations;

                swap(sequence[indexOp], sequence[indexOp2]);

                int mk = calculateMakespan(parameters, sequence, times);

                // if(indexOp>indexOp2) {
                //     calculateJMbyIndex(parameters, times, indexOp2-1 >= 0 ? indexOp2-1 : 0, sequence, M, J);
                //     mk = updateMakespan(parameters, times, indexOp2-1 >= 0 ? indexOp2-1 : 0, sequence, M, J);
                // }
                // else {
                //     calculateJMbyIndex(parameters, times, indexOp-1 >= 0 ? indexOp-1 : 0, sequence, M, J);
                //     mk = updateMakespan(parameters, times, indexOp-1 >= 0 ? indexOp-1 : 0, sequence, M, J);
                // }

                if(mk < bestMakespan){
                   // cout << "MELHORA" << endl;
                    bestTimes = times;
                    bestSequence = sequence;
                    bestMakespan = mk;
                }

            }

            offspring->makespan = bestMakespan;
            offspring->chromosome = bestSequence;
            offspring->endTimeOperations = bestTimes;
        }
    }

    return offspring;

}

void Genetic::evolve(int maxIterWithoutImprov){
    // Individuals used for crossover
    Individual *parent1;
    Individual *parent2;

    nbIterWithoutImprov = 1;
    int nbIterWithoutImprovDiv = 1;
    nbIter = 1;

    string temp;
    int place;
    clock_t debut = clock();

    // Child reference
    Individual *offspring = new Individual(parameters);
    int i =1;

    while (nbIterWithoutImprov < maxIterWithoutImprov) {
        cout << i << " " << nbIterWithoutImprov << endl;
        i++;

        parent1 = population->getIndividualBinT(); // Pick individual by binary tournament
        //parent2 = population->getIndividualBinT(); // Pick individual by binary tournament

        //offspring = crossoverOX(parent1, parent2); // OX crossover
        offspring = RR(parent1);
        offspring = swapGenetic(offspring);
        //if(!offspring->verifySequence()) continue;

        // Tries to add child to population
        place = population->addIndividual(offspring);

        if (place == -2) {
            return;
        }

        if (place == 0) { // A new best solution has been found
            nbIterWithoutImprov = 1;
            nbIterWithoutImprovDiv = 1;

            cout << population->getBestIndividual()->makespan << endl;
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

Individual* Genetic::crossoverOX(Individual *parent1, Individual *parent2) {
    
    // Beginning and end of the crossover zone
    unsigned int begin = rand() % parameters->numJobs * parameters->numJobs;
    unsigned int end = rand() % parameters->numJobs * parameters->numJobs;

    Individual *child = new Individual(parameters, 0);
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

Individual* Genetic::RR(Individual * parent1) {
    Individual *offspring = new Individual(parameters, 0);

    offspring->recopyIndividual(parent1);

    vector <int> :: iterator it;

    int m = rand() % parameters->numTools;
    //cout << "machine: " << m+1 << endl;

    vector<int> opsOfMachines;

    int erase = 0;

    for(int i =0; i < offspring->chromosome.size(); i++){
        int op = offspring->chromosome[i];
        int machineOp = decryptJobMachineIndex(op, parameters->numJobs).second;

        if(machineOp == m) {
            opsOfMachines.push_back(op);

            int index = find(offspring->chromosome.begin(), offspring->chromosome.end(), op) - offspring->chromosome.begin();
            offspring->chromosome.erase(offspring->chromosome.begin() + index);
        }
    }

    // show offspring

    // for(int i =0; i < offspring->chromosome.size(); i++){
    //     cout << offspring->chromosome[i] << " ";
    // } cout << endl;


    while(opsOfMachines.size() > 0){
        int op = opsOfMachines[0];
        int machinesPassed = 0;

        vector<int> sequence, times;
        vector<int> M, J, machinesBlock;

        //cout << "operação: " << op << endl;

        int bestMk = 99999;
        int bestIndex;
        vector<int> bestTimes = offspring->endTimeOperations;

        for(int i =0; i <= offspring->chromosome.size(); i++){
            sequence = offspring->chromosome;
            times = offspring->endTimeOperations;

            //cout << "machines passed: " << machinesPassed << " i: " << i << endl;

            if(machinesPassed == parameters->numTools-1 && decryptJobMachineIndex(offspring->chromosome[i], parameters->numJobs).second == m){
                machinesPassed = 0;
                continue;
            }

            if(machinesPassed == parameters->numTools-1){        
                //cout << "dentro do if: " << i << endl;
                
                sequence.insert(sequence.begin() + i, op);
                //show sequence

                // for(int i =0; i < sequence.size(); i++){
                //     cout << sequence[i] << " ";
                // } cout << endl;

                calculateJMbyIndex(parameters, times, i-1 > 0 ? i-1 : 0, sequence, M, J, machinesBlock);
                int mk = updateMakespan(parameters, times, i-1 > 0 ? i-1 : 0, sequence, M, J, machinesBlock);

                // calculateMakespan(parameters, sequence, times);

                if(mk < bestMk){
                    bestMk = mk;
                    bestIndex = i;
                    bestTimes = times;
                }

                machinesPassed = 1;
                continue;
            }

            //cout << "pos-if" << endl;

            int machineOp = decryptJobMachineIndex(offspring->chromosome[i], parameters->numJobs).second;
            if(machineOp != m) {
                machinesPassed++;
            }
        }

        offspring->chromosome.insert(offspring->chromosome.begin() + bestIndex, op);
        offspring->endTimeOperations = bestTimes;

        opsOfMachines.erase(opsOfMachines.begin());
    }

    // while(indexOfMachines.size() > 0){
    //     int indexRandIndex = rand() % indexOfMachines.size();
    //     int randIndex = indexOfMachines[indexRandIndex];

    //     int indexRandOp = rand() % opsOfMachines.size();
    //     int randOp = opsOfMachines[indexRandOp];

    //     indexOfMachines.erase(indexOfMachines.begin() + indexRandIndex);
    //     opsOfMachines.erase(opsOfMachines.begin() + indexRandOp);

    //     offspring->chromosome[randIndex] = randOp;
    // } 

    
    offspring->makespan = calculateMakespan(parameters, offspring->chromosome, offspring->endTimeOperations);
    return offspring;

}
