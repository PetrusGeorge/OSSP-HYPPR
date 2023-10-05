#include "BuscaLocal.h"
#include <math.h>


BuscaLocal::BuscaLocal(Parameters *parameters, vector<int> U, vector<int> endTimeOperations, int makespan){
    this->makespan = makespan;
    this->parameters = parameters;
    this->U = U;
    this->endTimeOperations = endTimeOperations;
}

bool BuscaLocal::redundancySwapAdjacent(unsigned op1, unsigned op2){
    pair<int, int> o1 = decryptJobMachineIndex(op1, parameters->numJobs);
    pair<int, int> o2 = decryptJobMachineIndex(op2, parameters->numJobs);

    if(o1.first != o2.first && o1.second != o2.second){
        return true;
    }

    return false;
}

bool BuscaLocal::redundancySwap(const vector<int>& U, int index1, int index2){
    pair<int, int> o1 = decryptJobMachineIndex(U[index1], parameters->numJobs);
    pair<int, int> o2 = decryptJobMachineIndex(U[index2], parameters->numJobs);

    if(o1.first != o2.first && o1.second != o2.second){
        return true;
    }

    for(int i = index1+1; i < index2; i++){
        pair<int, int> o = decryptJobMachineIndex(U[i], parameters->numJobs);
        if(o.first == o1.first || o.first == o2.first)
            return true;
        if(o.second == o1.second || o.second == o2.second)
            return true;
    }

    return false;

}


bool BuscaLocal::swap(){
    vector<int> sequence;
    vector<int> bestSequence = U;
    int bestMakespan = makespan;

    vector<int> J, M;
    vector<int> times;
    vector<int> bestTimes  = endTimeOperations;

    for(int i = 0; i < U.size()-1; i++){
        for(int j = i+1; j < U.size()-1; j++){

            sequence = U;
            times = endTimeOperations;

            if(j == i+1 && redundancySwapAdjacent(U[i], U[j]))
                continue;
            else if(j != i+1 && redundancySwap(U, i, j))
                continue;
            
            int aux = sequence[i];
            sequence[i] = sequence[j];
            sequence[j] = aux;

            calculateJMbyIndex(parameters, times, j-1, sequence, M, J);
            int mk = updateMakespan(parameters, times, j-1, sequence, M, J);

            if(mk < bestMakespan){
                bestMakespan = mk;
                bestSequence = sequence;
                bestTimes = times;
            }

        }
    }

    if(bestMakespan < makespan){
        U = bestSequence;
        makespan = bestMakespan;
        endTimeOperations = bestTimes;
        return true;
    }

    return false;
}

bool BuscaLocal::relocate(){
    vector<int> sequence;
    vector<int> bestSequence = U;
    int bestMakespan = makespan;

    vector<int> M, J;
    vector<int> times;
    vector<int> bestTimes = endTimeOperations;

    for(int i = 0; i < U.size(); i++){
        for(int j = 0; j < U.size(); j++){
            if(j == i)
                continue;

            sequence = U;
            times = endTimeOperations;
            int op = sequence[i];

            if(j == i+1 && redundancySwapAdjacent(U[i], U[j]))
                continue;

            sequence.erase(sequence.begin()+i);

            int mk;
            if(i>j) {
                sequence.insert(sequence.begin()+j, op);
                calculateJMbyIndex(parameters, times, j-1 >= 0 ? j-1 : 0, sequence, M, J);
                mk = updateMakespan(parameters, times, j-1 >= 0 ? j-1 : 0, sequence, M, J);
            }
            else {
                sequence.insert(sequence.begin()+j-1, op);
                calculateJMbyIndex(parameters, times, i-1 >= 0 ? i-1 : 0, sequence, M, J);
                mk = updateMakespan(parameters, times, i-1 >= 0 ? i-1 : 0, sequence, M, J);
            }

            if(mk < bestMakespan){
                bestMakespan = mk;
                bestSequence = sequence;
                bestTimes = times;
            }
        }
    }

    if(bestMakespan < makespan){
        U = bestSequence;
        makespan = bestMakespan;
        endTimeOperations = bestTimes;
        return true;
    }

    return false;

}

bool BuscaLocal::relocateBlock(int blockSize){
    vector<int> sequence;
    vector<int> bestSequence = U;
    int bestMakespan = makespan;

    vector<int> M, J;
    vector<int> times;
    vector<int> bestTimes = endTimeOperations;

    for(int i = 0; i < U.size()-blockSize; i++){
        for(int j = 0; j < U.size(); j++){
            int equal =0;

            for(int k = i; k < i+blockSize; k++){
                if(j == k){
                    equal = 1;
                    break;
                }
            }

            if(equal) continue;
            
            sequence = U;
            times = endTimeOperations;

            for(int w =0; w < blockSize; w++){
                if(i > j){
                    sequence.insert(sequence.begin()+j+w, sequence[i+w]);
                    sequence.erase(sequence.begin()+i+1+w);
                }
                else {
                    sequence.insert(sequence.begin()+j+1, sequence[i]);
                    sequence.erase(sequence.begin()+i);
                }
            }

            int mk;
            if(i>j) {
                calculateJMbyIndex(parameters, times, j-1 >= 0 ? j-1 : 0, sequence, M, J);
                mk = updateMakespan(parameters, times, j-1 >= 0 ? j-1 : 0, sequence, M, J);
            }
            else {
                calculateJMbyIndex(parameters, times, i-1 >= 0 ? i-1 : 0, sequence, M, J);
                mk = updateMakespan(parameters, times, i-1 >= 0 ? i-1 : 0, sequence, M, J);
            }

            if(mk < bestMakespan){
                bestMakespan = mk;
                bestSequence = sequence;
                bestTimes = times;
            }
        }
    }

    if(bestMakespan < makespan){
        U = bestSequence;
        makespan = bestMakespan;
        endTimeOperations = bestTimes;
        return true;
    }

    return false;

}

bool BuscaLocal::searchNeighborhood(unsigned int i) {
       
    if (i == 1) {
        return relocate();
    } else if (i == 2) {
        return relocateBlock(i);
    } else if (i == 3) {    
        return relocateBlock(i);
    } else if (i == 4) {
        return relocateBlock(i);
    } else { 
        return swap();
    }
    
}

void BuscaLocal::runSearchTotal() {

    bool foundBetter;
    double r;

    while (true) {
        foundBetter = false;
        for (unsigned int i = 1; i <= parameters->nViz && !foundBetter; i++) {
            foundBetter = searchNeighborhood(i);
        }
        if (!foundBetter) break;
    }

}