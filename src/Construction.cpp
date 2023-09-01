#include "Construction.h"
#include <math.h>

unsigned calculateMakespan(vector<unsigned> U, Parameters *parameters) {

    vector<int> M(parameters->numTools);//Sum of time in each machine
    vector<int> J(parameters->numJobs);//Sum of time in each job

    while(!U.empty()){
     
        int op = U[0];
        U.erase(U.begin());
        
        unsigned machineNumber = ceil(op * 1.0/ parameters->numTools * 1.0);
        unsigned jobNumber = op  - (machineNumber-1) * parameters->numTools;

        int acumulatedJ = (M[machineNumber-1] - J[jobNumber-1]) > 0 ? (M[machineNumber-1] - J[jobNumber-1]) : 0;
        int acumulatedM = (J[jobNumber-1] - M[machineNumber-1]) > 0 ? J[jobNumber-1] - M[machineNumber-1] : 0;

        J[jobNumber-1] += acumulatedJ + parameters->jobsToolsMatrix[jobNumber-1][machineNumber-1];
        M[machineNumber-1] += acumulatedM + parameters->jobsToolsMatrix[jobNumber-1][machineNumber-1];

        cout << "J[" << jobNumber-1 << "]: " << J[jobNumber-1] << " M[" << machineNumber-1 << "]: " << M[machineNumber-1] << endl << endl;        
    }

    return *max_element(M.begin(), M.end()) ;
}

void Construction(Parameters *parameters) {
    // ...
}