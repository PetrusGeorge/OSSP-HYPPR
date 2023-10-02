#include "BuscaLocal.h"
#include <math.h>

bool BuscaLocal::redundancySwapAdjacent(Parameters *parameters, unsigned op1, unsigned op2){
    pair<int, int> o1 = decryptJobMachineIndex(op1, parameters->numJobs);
    pair<int, int> o2 = decryptJobMachineIndex(op2, parameters->numJobs);

    if(o1.first != o2.first && o1.second != o2.second){
        return true;
    }

    return false;
}

bool BuscaLocal::redundancySwap(Parameters *parameters, const vector<int>& U, int index1, int index2){
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


bool BuscaLocal::swap(Parameters *parameters, vector<int>& U, int makespan){
    vector<int> sequence;
    vector<int> bestSequence = U;
    int bestMakespan = makespan;

    for(int i = 0; i < U.size()-1; i++){
        for(int j = i+1; j < U.size()-1; j++){
            sequence = U;
            if(j == i+1 && redundancySwapAdjacent(parameters, U[i], U[j]))
                continue;
            else if(j != i+1 && redundancySwap(parameters, U, i, j))
                continue;
            
            int aux = sequence[i];
            sequence[i] = sequence[j];
            sequence[j] = aux;

            int mk = calculateMakespan(parameters, sequence);
            if(mk < bestMakespan){
                bestMakespan = mk;
                bestSequence = sequence;
            }

        }
    }

    if(bestMakespan < makespan){
        U = bestSequence;
        return true;
    }

    return false;
}

bool BuscaLocal::relocate(Parameters *parameters, vector<int>& U, int makespan){
    vector<int> sequence;
    vector<int> bestSequence = U;
    int bestMakespan = makespan;

    for(int i = 0; i < U.size()-1; i++){
        for(int j = 0; j < U.size()-1; j++){
            if(j == i)
                continue;
            sequence = U;
            int op = sequence[i];

            if(j == i+1 && redundancySwapAdjacent(parameters, U[i], U[j]))
                continue;

            sequence.erase(sequence.begin()+i);
            if(i>j) sequence.insert(sequence.begin()+j, op);
            else sequence.insert(sequence.begin()+j-1, op);

            int mk = calculateMakespan(parameters, sequence);
            if(mk < bestMakespan){
                bestMakespan = mk;
                bestSequence = sequence;
            }
        }
    }

    if(bestMakespan < makespan){
        U = bestSequence;
        return true;
    }

    return false;

}

bool BuscaLocal::relocate2(Parameters *parameters, vector<int>& U, int makespan){
    vector<int> sequence;
    vector<int> bestSequence = U;
    int bestMakespan = makespan;

    for(int i = 0; i < U.size()-2; i++){
        for(int j = 0; j < U.size()-1; j++){
            if(j == i || j == i+1)
                continue;
            
            sequence = U;
            int op1 = sequence[i];
            int op2 = sequence[i+1];

            for(int w = 0; w <2; w++){
                sequence.erase(sequence.begin()+i);
            }

            if(i > j){
                sequence.insert(sequence.begin()+j, op2);
                sequence.insert(sequence.begin()+j, op1);
            }
            else {
                sequence.insert(sequence.begin()+j-1, op2);
                sequence.insert(sequence.begin()+j-1, op1);
            }
            

            int mk = calculateMakespan(parameters, sequence);

            if(mk < bestMakespan){
                bestMakespan = mk;
                bestSequence = sequence;
            }
        }
    }

    if(bestMakespan < makespan){
        U = bestSequence;
        return true;
    }

    return false;

}