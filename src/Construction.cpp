#include "Construction.h"
#include <math.h>

pair<unsigned, unsigned> decryptJobMachineIndex(unsigned op, unsigned numJobs){
    unsigned machineNumber = ceil(op * 1.0/ numJobs * 1.0);
    unsigned jobNumber = op  - (machineNumber-1) * numJobs;

    return make_pair(jobNumber-1, machineNumber-1);
}

int calculateMakespan(Parameters *parameters, vector <int> U, vector<int>& endTimeOperations) {
    vector<int> M(parameters->numTools, 0);//Machines acumulated time
    vector<int> J(parameters->numJobs, 0);//Jobs acumulated time
    //U is the sequence of operations
    while(!U.empty()){
     
        int op = U[0];
        U.erase(U.begin());
        
        pair<unsigned, unsigned> indexJM = decryptJobMachineIndex(op, parameters->numJobs);

        int acumulatedJ = (M[indexJM.second] - J[indexJM.first]) > EPSILON ? (M[indexJM.second] - J[indexJM.first]) : 0;
        int acumulatedM = (J[indexJM.first] - M[indexJM.second]) > EPSILON ? J[indexJM.first] - M[indexJM.second] : 0;

        J[indexJM.first] += acumulatedJ + parameters->jobsToolsMatrix[indexJM.first][indexJM.second];
        M[indexJM.second] += acumulatedM + parameters->jobsToolsMatrix[indexJM.first][indexJM.second];

        endTimeOperations[op-1] = M[indexJM.second];

        //cout << "J[" << jobNumber-1 << "]: " << J[jobNumber-1] << " M[" << machineNumber-1 << "]: " << M[machineNumber-1] << endl << endl;        
    }

    return *max_element(M.begin(), M.end()) ;
}

/* int calculateMakespan(vector<int> U, Parameters *parameters, vector<int>& M, vector<int>& J) {

    M = vector<int>(parameters->numTools, 0);
    J = vector<int>(parameters->numJobs, 0);

    while(!U.empty()){
     
        int op = U[0];
        U.erase(U.begin());
        
        pair<unsigned, unsigned> indexJM = decryptJobMachineIndex(op, parameters->numJobs);

        int acumulatedJ = (M[indexJM.second] - J[indexJM.first]) > EPSILON ? (M[indexJM.second] - J[indexJM.first]) : 0;
        int acumulatedM = (J[indexJM.first] - M[indexJM.second]) > EPSILON ? J[indexJM.first] - M[indexJM.second] : 0;

        J[indexJM.first] += acumulatedJ + parameters->jobsToolsMatrix[indexJM.first][indexJM.second];
        M[indexJM.second] += acumulatedM + parameters->jobsToolsMatrix[indexJM.first][indexJM.second];

        //cout << "J[" << jobNumber-1 << "]: " << J[jobNumber-1] << " M[" << machineNumber-1 << "]: " << M[machineNumber-1] << endl << endl;        
    }

    return *max_element(M.begin(), M.end()) ;
} */

inline void updateMJ(Parameters *parameters, vector<int>& M, vector<int>& J, unsigned op, vector<int>& endTimeOperations){
    pair<unsigned, unsigned> indexJM = decryptJobMachineIndex(op, parameters->numJobs);

    int acumulatedJ = (M[indexJM.second] - J[indexJM.first]) > EPSILON ? (M[indexJM.second] - J[indexJM.first]) : 0;
    int acumulatedM = (J[indexJM.first] - M[indexJM.second]) > EPSILON ? J[indexJM.first] - M[indexJM.second] : 0;

    J[indexJM.first] += acumulatedJ + parameters->jobsToolsMatrix[indexJM.first][indexJM.second];
    M[indexJM.second] += acumulatedM + parameters->jobsToolsMatrix[indexJM.first][indexJM.second];

    endTimeOperations[op-1] = M[indexJM.second];
}

int updateMakespan(Parameters *parameters, vector<int>& endTimeOperations, int index, vector<int> U, vector<int>& M, vector<int>& J){
    for(int i = index+1; i < U.size(); i++){
        unsigned op = U[i];
        updateMJ(parameters, M, J, op, endTimeOperations);
    }

    return *max_element(M.begin(), M.end());
}

void calculateJMbyIndex(Parameters *parameters, const vector<int>& endTimeOperations, int index, vector<int> U, vector<int>& M, vector<int>& J){
    M = vector<int>(parameters->numTools, 0);
    J = vector<int>(parameters->numJobs, 0);
    int j = 0, m = 0;
    for(int i = index; i >= 0; i--){
        unsigned op = U[i];
        pair<unsigned, unsigned> indexJM = decryptJobMachineIndex(op, parameters->numJobs);

        if(j == parameters->numJobs && m == parameters->numTools){
            break;
        }

        if(J[indexJM.first] == 0){
            J[indexJM.first] = endTimeOperations[op-1];
            j++;
        }
        if(M[indexJM.second] == 0){
            M[indexJM.second] = endTimeOperations[op-1];
            m++;
        }
    }
}

inline int calculateMakespanOp(Parameters *parameters, vector<int> M, const vector<int> J, unsigned op){
    pair<unsigned, unsigned> indexJM = decryptJobMachineIndex(op, parameters->numJobs);

    int acumulatedM = (J[indexJM.first] - M[indexJM.second]) > EPSILON ? J[indexJM.first] - M[indexJM.second] : 0;

    M[indexJM.second] += acumulatedM + parameters->jobsToolsMatrix[indexJM.first][indexJM.second];

    return *max_element(M.begin(), M.end()) ;
}

double EMC(Parameters *parameters,const vector<vector<unsigned>>& P, unsigned op){

    double lbJobs = 0;
    double lbMachines = 0;

    unsigned machineNumber = ceil(op * 1.0/ parameters->numJobs * 1.0);
    unsigned jobNumber = op  - (machineNumber-1) * parameters->numJobs;

    for(int i = 0; i < parameters->numJobs; i++){
        double lb = 0;
        for(int j = 0; j < parameters->numTools; j++){

            if(i == jobNumber-1 && j == machineNumber-1){
                j++;
                continue;
            }
          lb += P[i][j];
        }
        lbJobs = max(lbJobs, lb);
    }

    for(int i = 0; i < parameters->numTools; i++){
        double lb = 0;
        for(int j = 0; j < parameters->numJobs; j++){
            
            if(i == jobNumber-1 && j == machineNumber-1){
                j++;
                continue;
            }
          lb += P[i][j];
        }
        lbMachines = max(lbMachines, lb);
    }

    return max(lbJobs, lbMachines);
}

unsigned int BICH_MIH(Parameters *parameters, vector<int>& endTimeOperations, vector<int>& chromose){
    // parameters->numJobs = 3;
    // parameters->numTools = 3;
    // parameters->positionsOffspring.resize(3);
    vector<int> M(parameters->numTools, 0);
    vector<int> J(parameters->numJobs, 0);
    endTimeOperations = vector<int>(parameters->numJobs * parameters->numTools, 0);
    vector<vector<unsigned>> P = parameters->jobsToolsMatrix;
    vector<int> currentSequence;

    unsigned operationsNum = parameters->numJobs * parameters->numTools;

    while(currentSequence.size() < operationsNum){
        vector<pair<double, int>> ops;

        for(int i = 0; i < parameters->numJobs; i++){
            for(int j = 0; j < parameters->numTools; j++){
                if(P[i][j] != 0){
                    unsigned op = j * parameters->numJobs + i + 1;
                    double bich_mih = (EMC(parameters, P, op) + calculateMakespanOp(parameters, M, J, op))*(1 - parameters->alpha);

                    bich_mih += (parameters->alpha * (J[i] - M[j] > EPSILON ? J[i] - M[j] : 0));

                    ops.push_back(make_pair(bich_mih, op));
                }
            }
        }

        sort(ops.begin(), ops.end());

        double beta = ((double) rand() / RAND_MAX) + 0.00001;//0.00001 to avoid 0
        int choseIndex = rand() % ((int) ceil(beta * ops.size()));//random function that has more chance to be a lower number
        int choseOp = ops[choseIndex].second;
        currentSequence.push_back(choseOp);
        updateMJ(parameters, M, J, choseOp, endTimeOperations);

        pair<unsigned, unsigned> indexJM = decryptJobMachineIndex(choseOp, parameters->numJobs);
        
        P[indexJM.first][indexJM.second] = 0;
    }

    chromose = currentSequence;
    return *max_element(M.begin(), M.end());
}