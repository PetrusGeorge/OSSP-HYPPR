#ifndef CONSTRUCTION_H
#define CONSTRUCTION_H

using namespace std;

#include "Parameters.h"

inline unsigned calculateMakespan(vector<unsigned> U, Parameters *parameters);
inline unsigned calculateMakespan(vector<unsigned> U, Parameters *parameters, vector<int>& M, vector<int>& J);
vector<int> BICH_MIH(Parameters *parameters);

#endif