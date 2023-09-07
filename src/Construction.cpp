#include "Construction.h"
#include <math.h>
#include <queue>

int calculateMakespan(vector<int> U, Parameters *parameters) {

    vector<int> M(parameters->numTools, 0);
    vector<int> J(parameters->numJobs, 0);

    while(!U.empty()){
     
        int op = U[0];
        U.erase(U.begin());
        
        unsigned machineNumber = ceil(op * 1.0/ parameters->numJobs * 1.0);
        unsigned jobNumber = op  - (machineNumber-1) * parameters->numJobs;

        int acumulatedJ = (M[machineNumber-1] - J[jobNumber-1]) > 0 ? (M[machineNumber-1] - J[jobNumber-1]) : 0;
        int acumulatedM = (J[jobNumber-1] - M[machineNumber-1]) > 0 ? J[jobNumber-1] - M[machineNumber-1] : 0;

        J[jobNumber-1] += acumulatedJ + parameters->jobsToolsMatrix[jobNumber-1][machineNumber-1];
        M[machineNumber-1] += acumulatedM + parameters->jobsToolsMatrix[jobNumber-1][machineNumber-1];

        //cout << "J[" << jobNumber-1 << "]: " << J[jobNumber-1] << " M[" << machineNumber-1 << "]: " << M[machineNumber-1] << endl << endl;        
    }

    return *max_element(M.begin(), M.end()) ;
}

int calculateMakespan(vector<int> U, Parameters *parameters, vector<int>& M, vector<int>& J) {

    M = vector<int>(parameters->numTools, 0);
    J = vector<int>(parameters->numJobs, 0);

    while(!U.empty()){
     
        int op = U[0];
        U.erase(U.begin());
        
        unsigned machineNumber = ceil(op * 1.0/ parameters->numJobs * 1.0);
        unsigned jobNumber = op  - (machineNumber-1) * parameters->numJobs;

        int acumulatedJ = (M[machineNumber-1] - J[jobNumber-1]) > 0 ? (M[machineNumber-1] - J[jobNumber-1]) : 0;
        int acumulatedM = (J[jobNumber-1] - M[machineNumber-1]) > 0 ? J[jobNumber-1] - M[machineNumber-1] : 0;

        J[jobNumber-1] += acumulatedJ + parameters->jobsToolsMatrix[jobNumber-1][machineNumber-1];
        M[machineNumber-1] += acumulatedM + parameters->jobsToolsMatrix[jobNumber-1][machineNumber-1];

        //cout << "J[" << jobNumber-1 << "]: " << J[jobNumber-1] << " M[" << machineNumber-1 << "]: " << M[machineNumber-1] << endl << endl;        
    }

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

vector<int> BICH_MIH(Parameters *parameters){
    vector<int> M;
    vector<int> J;
    vector<vector<unsigned>> P = parameters->jobsToolsMatrix;
    vector<int> currentSequence;

    unsigned operationsNum = parameters->numJobs * parameters->numTools;

    while(currentSequence.size() < operationsNum){
        priority_queue <pair<double, int>> ops;

        for(int i = 0; i < parameters->numJobs; i++){
            for(int j = 0; j < parameters->numTools; j++){
                if(P[i][j] != 0){
                    unsigned op = j * parameters->numJobs + i + 1;
                    double bich_mih = (EMC(parameters, P, op) + calculateMakespan(currentSequence, parameters, M, J))*(1 - parameters->alpha);

                    bich_mih += (parameters->alpha * (J[i] - M[j] > 0 ? J[i] - M[j] : 0));

                    ops.push(make_pair(bich_mih, op));
                }
            }
        }

        unsigned bestOp = ops.top().second;
        currentSequence.push_back(bestOp);

        unsigned machineNumber = ceil(bestOp * 1.0/ parameters->numJobs * 1.0);
        unsigned jobNumber = bestOp  - (machineNumber-1) * parameters->numJobs;
        
        P[jobNumber-1][machineNumber-1] = 0;
    }

    return currentSequence;
}