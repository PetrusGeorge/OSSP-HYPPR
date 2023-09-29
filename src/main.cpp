/**
 * @file main.cpp
 *
 * @author Jordana Mecler
 *
 * Contact: jmecler@inf.puc-rio.br
 *
 */

//#include <Windows.h>
//#include <stdlib.h>
#include "CommandLine.h"
#include "Parameters.h"
#include "Population.h"
#include "Genetic.h"
#include "Util.h"
#include <chrono>

#include "Construction.h"

#define min(a, b) ((a < b) ? a : b)

using namespace std;

int main(int argc, char* argv[]) {

    double bestCost;
    double averageCost, averageTime;

    unsigned int runs = 10;
    Parameters *parameters;
    Population *population;
    clock_t nb_ticks_allowed;

    try {
        // Read arguments from command line and set problem parameters
        CommandLine commandLine(argc, argv);
        if (!commandLine.isValid())
            throw string(
                    "(ERROR) Wrong command line. Format: [instance_path] [instance_names(-)] [solution_path] [seed] [population_size] [max_population_size] [number_elite] [number_close_individuals] [max_diversify]");
        nb_ticks_allowed = commandLine.getCpuTime() * CLOCKS_PER_SEC;
        parameters = new Parameters(commandLine.getSeed(), commandLine.getInstancesPaths(), commandLine.getInstancesNames(),
                            commandLine.getSolutionPath(), commandLine.getPopulationSize(), commandLine.getMaxPopulationSize(),
                            commandLine.getNumberElite(), commandLine.getNumberCloseIndividuals(), commandLine.getMaxDiversify());
//        runs = commandLine.getRuns();

        //cout << "Starting " << commandLine.getInstancesPaths() << "/" << commandLine.getInstancesNames() <<  endl;
        //cout << endl;

        double cost, totalAvgCost = 0, totalBestCost = 0;
        double totalTime = 0;
        double totalImprovPrim = 0, totalImprovSec = 0;

        // Run algorithm for every file in the set of instances
        for (const string &item : parameters->files) {
//            const clock_t totalStartTime = clock();
//            bestCost = numeric_limits<double>::infinity();
//            averageCost = 0;
//            averageTime = 0;
//            parameters->improvesPrimary = 0;
//            parameters->improvesSecondary = 0;

            // Run the algorithm on the same file a number of 'runs' times
            // to generate average results
//            for (unsigned int i = 0; i < runs; i++) {
                const clock_t startTime = clock();
                
                parameters->terminate = false;
                parameters->readFile(item);

                parameters->setAuxiliaryParams();
                                
                parameters->cpuTime = cpuTime();
                // Start population
                population = new Population(parameters);
                // cout << population->getBestIndividual()->solutionCost.evaluation << endl;
                // cout << cpuTime() - parameters->cpuTime << endl << endl;

                // for (int i = 0; i < population->getBestIndividual()->chromosome.size(); i++) {
                //     cout << population->getBestIndividual()->chromosome[i] << " ";
                // }
                
                // cout << endl << endl;

                // getchar();

                /*
                vector<int> teste;
                vector<int> endTimeOperations;

                teste = BICH_MIH(parameters, endTimeOperations);
                */

               vector<int> teste = {9,3,5,6,4,8,7,2,1};
               vector<int> endTimeOperations(9, 0);
               calculateMakespan(parameters, endTimeOperations, teste);

                for(auto a : teste){
                    cout << a << " ";
                }cout << endl;

                for(int i =0; i < endTimeOperations.size(); i++){
                    cout << endTimeOperations[i] << " ";
                }

                vector<int> M, J;
                cout << endl;

                cout << updateMakespan(parameters, endTimeOperations, 3, teste, M, J) << endl;

                // show M e J vectors

                for(auto a : M){
                    cout << a << " ";
                }cout << endl;

                for(auto a : J){
                    cout << a << " ";
                }cout << endl;
                // Run genetic algorithm
                /* Genetic solver(parameters, population, nb_ticks_allowed, false);
                solver.evolve(min(1001, 1000));
                //solver.evolve(min(parameters->numJobs * 20, 1000));

                cost = population->getBestIndividual()->solutionCost.evaluation;
                totalTime = (float(clock() - startTime) / CLOCKS_PER_SEC);
                double finalTime = cpuTime() - parameters->cpuTime; */
                
                // string instance_name = parameters->instancesNames.substr(parameters->instancesNames.rfind("/") + 1);
                // ofstream outClientFile("Resultados/Teste20010.csv", ios::app);
                // outClientFile << instance_name << "," << cost << "," << finalTime << "," << 4 << ",";
                // for (int i = 0; i < parameters->numJobs; i++) {
                //     outClientFile << population->getBestIndividual()->chromosome[i] << " ";
                // } 
                // outClientFile << endl;
                
                delete population;

                // Write result to solution file
//                outputFile.open(parameters->solutionPath, ofstream::out | ofstream::app);
//                outputFile << cost << ","  << totalTime << endl;
//                outputFile.close();

//                averageCost += cost;
//                if (cost < bestCost) {
//                    bestCost = cost;
//                }
//            cout << "Run " << i << endl;
//            }
//            averageCost /= runs;
//            averageTime /= runs;

//            totalAvgCost += averageCost;
//            totalBestCost += bestCost;
//            totalTime += averageTime;

//            parameters->improvesPrimary /= runs;
//            parameters->improvesSecondary /= runs;

//            totalImprovPrim += parameters->improvesPrimary;
//            totalImprovSec += parameters->improvesSecondary;

        }

//        totalAvgCost /= parameters->files.size();
//        totalBestCost /= parameters->files.size();
//        totalTime /= parameters->files.size();

//        totalImprovPrim /= parameters->files.size();
//        totalImprovSec /= parameters->files.size();

        // Write average result to solution file
//        outputFile.open(parameters->solutionPath, ofstream::out | ofstream::app);
//        outputFile << totalBestCost << "," << totalAvgCost << ","  << totalTime << endl;
//        outputFile.close();

        delete parameters;

        //cout << "Finishing " << commandLine.getInstancesPaths() << "/" <<  commandLine.getInstancesNames() << " without errors" << endl;
        //cout << endl;

        return 0;
    }
    catch (const string &error) {
        cout << error << endl;
        return 0;
    }
}