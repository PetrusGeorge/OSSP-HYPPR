#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <time.h>

using namespace std;

class Parameters {

    private:
        void setMethodParams();

    public:

        bool terminate;

        // Number of iterations on genetic algorithm
        unsigned int nbIter;
        
        // Problem Parameters
        string instancePath;
        
        // Instance variables
        unsigned int numJobs;
        unsigned int numTools;
        unsigned int maxCapacity;
        vector<vector<unsigned int> > jobsToolsMatrix;

        unsigned int nIterNeighborhood;

        // Alpha for BICH_MIH
        double alpha;

        // Method Parameters
        unsigned int populationSize;
        unsigned int maxPopulationSize;
        unsigned int numberCloseIndividuals;
        unsigned int maxDiversify;
        unsigned int numberElite;

        vector<bool> positionsOffspring; ///< temporary structure for crossover operation



        Parameters(string instancePath);
};

#endif //PARAMETERS_H