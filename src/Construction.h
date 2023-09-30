#ifndef CONSTRUCTION_H
#define CONSTRUCTION_H

#define EPSILON 0.0000001

using namespace std;

#include "Parameters.h"

vector<int> BICH_MIH(Parameters *parameters, vector<int>& endTimeOperations);
int calculateJMbyIndex(Parameters *parameters, const vector<int>& endTimeOperations, int index, vector<int> U, vector<int>& M, vector<int>& J);
int calculateMakespan(Parameters *parameters, vector<int>& endTimeOperations, vector<int> U);
void updateMakespan(Parameters *parameters, vector<int>& M, vector<int>& J, int index, vector<int> U, vector<int>& endTimeOperations);
#endif