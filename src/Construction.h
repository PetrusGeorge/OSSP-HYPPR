#ifndef CONSTRUCTION_H
#define CONSTRUCTION_H

#define EPSILON 0.0000001

using namespace std;

#include "Parameters.h"

vector<int> BICH_MIH(Parameters *parameters, vector<int>& endTimeOperations);
int calculateMakespan(Parameters *parameters, vector<int> U);
pair<unsigned, unsigned> decryptJobMachineIndex(unsigned op, unsigned numJobs);

#endif