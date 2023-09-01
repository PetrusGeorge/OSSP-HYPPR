/**
 * @class Parameters
 *
 * @brief Store main parameters
 *
 * Read input instances and store parameters
 * that are used throughout the algorithm.
 *
 * @author Jordana Mecler
 *
 * Contact: jmecler@inf.puc-rio.br
 *
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <iostream>
#include "math.h"
#include <time.h>
#include <algorithm>

using namespace std;

class Parameters {
public:
    // Random Seed
    int seed;

    bool terminate;
    double cpuTime;

    // Parameters to be displayed in solution file
    unsigned int nbIter;
    unsigned int improvesPrimary;
    unsigned int improvesSecondary;
    unsigned int ties;

    unsigned int nViz;

    // Problem Parameters
    string instancesPaths;
    string instancesNames;
    string solutionPath;

    vector<string> files;
    unsigned int numJobs;
    unsigned int numTools;
    unsigned int maxCapacity;
    vector<vector<unsigned int> > jobsToolsMatrix;

    // Method Parameters
    unsigned int populationSize;
    unsigned int maxPopulationSize;
    unsigned int numberElite;
    unsigned int numberCloseIndividuals;
    unsigned int maxDiversify;
    unsigned int jobsMovidos;

    vector<vector<unsigned int> > L; ///< auxiliary matrix for KTNS
    vector<unsigned int> W_n; ///< auxiliary vector for KTNS
    vector<bool> used; ///< maps if tool is loaded in KTNS
    vector<unsigned int> modifiedSolution; ///< modified solution structure for local search
    vector<unsigned int> tempSolution; ///< temporary structure for local search
    vector<bool> positionsOffspring; ///< temporary structure for crossover operation
    // vector<vector<unsigned int> > orderLSi; ///< order of search in local search
    // vector<vector<unsigned int> > orderLSii; 
    // vector<vector<unsigned int> > orderLSs;

    vector<unsigned int> listaSwap;

    vector<unsigned int> moves;
    vector<unsigned int> ktnsCache;
//    vector<bool> improvedRemoving;
//    vector<vector<unsigned int>> jobsDistance;
    vector<vector<unsigned int> > loadedMatrix; ///< instance matrix after KTNS

    /// Initialize the method params.
    void setMethodParams();

    /// Initialize some method auxiliary parameters.
    void setAuxiliaryParams();

    /// Get the set of all files with a certain path and prefix.
    void getFiles();

    /**
     * Read instance file and assign its values to variables.
     *
     * @param file path of the file to be read
     */
    void readFile(const string &file);

    /**
     * Constructor assign values from command line and instances read.
     *
     * @param seed random seed
     * @param instancePath path to the set of instances to be processed
     * @param instancesNames prefix of the instances in the set
     * @param solutionPath path to solution file
     * @param populationSize minimum population size
     * @param maxPopulationSize maximum subpopulation (beyond population size)
     * @param numberElite number of elite individuals considered for fitness
     * @param numberCloseIndividuals number of close individuals considered for distance calculation
     * @param maxDiversify maximum iterations without improvement resulting in diversification
     */
    Parameters(int seed, string instancePath, string instancesNames, string solutionPath, unsigned int populationSize, unsigned int maxPopulationSize, unsigned int numberElite, unsigned int numberCloseIndividuals, unsigned int maxDiversify);

    int ** C;
    int ** F;
    int ** P;

    double tempoIdle;
    pair<double, int> viz1;
    pair<double, int> viz2;
    pair<double, int> viz3;

    /// Destructor
    ~Parameters();
};

#endif //PARAMETERS_H