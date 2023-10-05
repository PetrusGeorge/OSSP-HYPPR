#ifndef CONSTRUCTION_H
#define CONSTRUCTION_H

#define EPSILON 0.0000001

using namespace std;

#include "Parameters.h"

unsigned int BICH_MIH(Parameters *parameters, vector<int>& endTimeOperations, vector<int>& chromosome);
int calculateMakespan(Parameters *parameters, vector<int> U);
pair<unsigned, unsigned> decryptJobMachineIndex(unsigned op, unsigned numJobs);
int updateMakespan(Parameters *parameters, vector<int>& endTimeOperations, int index, vector<int> U, vector<int>& M, vector<int>& J);
void calculateJMbyIndex(Parameters *parameters, const vector<int>& endTimeOperations, int index, vector<int> U, vector<int>& M, vector<int>& J);

#endif