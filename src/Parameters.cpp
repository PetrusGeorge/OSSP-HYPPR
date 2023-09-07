/**
 * @file Parameters.cpp
 *
 * Implements Parameters class methods
 *
 * @author Jordana Mecler
 *
 * Contact: jmecler@inf.puc-rio.br
 *
 */

//#include <Windows.h>
#include <iterator>
#include <glob.h>
#include "Parameters.h"

void Parameters::setMethodParams() {
    populationSize = 20;
    maxPopulationSize = 40;
    numberElite = 10;
    numberCloseIndividuals = 3;
    maxDiversify = 1000;
    terminate = false;
    jobsMovidos = 4;

    // improvesPrimary = 0;
    // improvesSecondary = 0;
    // ties = 0;
}

void Parameters::setAuxiliaryParams() {
    /*vector<vector<unsigned int> > L(numTools, vector<unsigned int>(numJobs));
    vector<unsigned int> W_n(numTools);
    vector<bool> used(numTools);
    vector<bool> positionsOffspring(numJobs);
    vector<bool> improvedRemoving(numJobs);
//    vector<vector<unsigned int>> jobsDistance(numJobs, vector<unsigned int>(numJobs));
    vector<vector<unsigned int> > loadedMatrix(numJobs, vector<unsigned int>(numTools));
    this->L = L;
    this->W_n = W_n;
    this->used = used;
    this->positionsOffspring = positionsOffspring;
//    this->improvedRemoving = improvedRemoving;
    this->loadedMatrix = loadedMatrix;*/

    // cout << numJobs << endl;
    vector<bool> positionsOffspring(numJobs);
    this->positionsOffspring = positionsOffspring;

    if (numJobs != 500) {
        nViz = 5;
    } else {
        nViz = 4;
    }

    this->moves = {1, 2, 3, 4, 5};
    // if (this->orderLSi.empty()) {
    //     for (unsigned int i = 0; i < numJobs; i++) {
    //         for (unsigned int j = 0; j < numJobs; j++) {
    //             if (!((int)i - 1 <= (int)j && (int)j < (int)i + 1)) {
    //                 this->orderLSi.push_back({i,j});
    //             }
    //         }
    //     }
    // }

    // if (this->orderLSs.empty()) {
    //     for (unsigned int i = 0; i < numJobs; i++) {
    //         for (unsigned int j = i + 1; j < numJobs; j++) {
    //             this->orderLSs.push_back({i,j});
    //         }
    //     }
    // }

    // if (this->orderLSii.empty()) {
    //     for (unsigned int i = 0; i < numJobs + 1 - 2; i++) {
    //         for (unsigned int j = 0; j < numJobs; j++) {
    //             if (!((int)i - 1 <= (int)j && (int)j < (int)i + 2)) {
    //                 this->orderLSii.push_back({i,j});
    //             }
    //         }
    //     }
    // }

    if (this->listaSwap.empty()) {
        for (unsigned int i = 0; i < numJobs; i++) {
            this->listaSwap.push_back(i);
        }
    }

}

void Parameters::getFiles() {
    // Get the set of instances given the path and prefix determined

    vector<string> prefixList;
    if (instancesNames.find('L') != string::npos) {
        string val;
        stringstream prefixStream;
        prefixStream.str(instancesNames);
        while (getline(prefixStream, val, '-')) {
            val.append("-");
            prefixList.push_back(val);
        }
    }
    else {
        prefixList.push_back(instancesNames);
    }

    for (const string &item : prefixList) {
        string pattern(instancesPaths);
        pattern.append("/").append(item).append("*");
//        WIN32_FIND_DATAA data;
//        HANDLE hFind;
//        cout << pattern << endl;
//        if ((hFind = FindFirstFileA(pattern.c_str(), &data)) != INVALID_HANDLE_VALUE) {
//            do {
//                string file = instancesPaths;
//                file.append("\\").append(data.cFileName);
//                files.push_back(file);
//            } while (FindNextFileA(hFind, &data) != 0);
//            FindClose(hFind);
//        }
        glob_t globResult;
        glob(pattern.c_str(), GLOB_TILDE, nullptr, &globResult);
        for (unsigned int i = 0; i < globResult.gl_pathc; ++i) {
            files.push_back(string(globResult.gl_pathv[i]));
        }
        globfree(&globResult);
    }

    if (files.empty()) {
        cout << "No files found matching the specified prefix. Exiting." << endl;
    }
}

double calculateAlpha(unsigned numTools) {
    // Calculate alpha parameter
    double a;
    switch (numTools){
        case 3:
        case 4:
        case 15:
            a = 0.89;
            break;
        case 5:
            a = 0.59;
            break;
        case 6:
            a = 0.11;
            break;
        case 7:
        case 8:
        case 20:
            a = 0.22;
            break;
        case 9:
        case 10:
            a = 0.67;
            break;
        default:
            a = 0.1;
    }

    return a;
}

void Parameters::readFile(const string &file) {
    // Read instance file

    ifstream inFile(file);

    if (inFile.good()) {
        string sLine;
        getline(inFile, sLine);
        numJobs = (unsigned int)stoi(sLine);
        getline(inFile, sLine);
        numTools = (unsigned int)stoi(sLine);

        alpha = calculateAlpha(numTools);
        
        unsigned int j = 0;
        vector<vector<unsigned int> > lines(numJobs, vector<unsigned int>(numTools, 0));
        while (getline(inFile, sLine)) {
            cout << "Instancia: " << sLine << endl;
            istringstream iss(sLine);
            vector<string> tokens {(istream_iterator<string>(iss)), istream_iterator<string>() };;
            for (unsigned int i = 0; i < tokens.size(); i++) {                
                cout << "TOKEN: " << tokens[i] << endl;
                lines[j][i] = (unsigned int)stoi(tokens[i]);
            }
            j++;
        }
        jobsToolsMatrix = lines;
    }

    cout << "MATRIZ DE CUSTOS\n";
    for(int  i = 0; i < jobsToolsMatrix.size(); i++){
        for(int j = 0; j < jobsToolsMatrix[i].size(); j++){
            cout << jobsToolsMatrix[i][j] << " "; //PROCESSAMENTO DO JOB I NA MAQUINA J
        }
        cout << endl;
    }
    cout << endl;
    inFile.close();

    C = new int * [numJobs + 1];
    F = new int * [numJobs + 1];
    P = new int * [numJobs + 1];    

    for (int i = 0; i < numJobs + 1; i++) {
        C[i] = new int [numTools + 1];
        F[i] = new int [numTools + 1];
        P[i] = new int [numTools + 1];
    }

    for(int j = 0; j < numTools + 1; j++) {
        F[0][j] = 0;
        P[0][j] = 0;
        C[0][j] = 0;
    }

    //C = TEMPO EM QUE A I-ÉSIMA TAREFA TERMINA NA MAQUINA J
    //DESCOBRIR O QUE É O F E O P
    //PROCURAR UMA CONSTRUÇAO
    
}

Parameters::Parameters(int seed, string instancesPaths, string instancesNames, string solutionPath, unsigned int populationSize, unsigned int maxPopulationSize, unsigned int numberElite, unsigned int numberCloseIndividuals, unsigned int maxDiversify) {
    //this->instancesPaths = instancesPaths;
    this->instancesNames = instancesNames;
    //this->solutionPath = solutionPath;
    files.push_back(instancesNames);
//    this->populationSize = populationSize;
//    this->maxPopulationSize = maxPopulationSize;
//    this->numberElite = numberElite;
//    this->numberCloseIndividuals = numberCloseIndividuals;
//    this->maxDiversify = maxDiversify;
    /*this->seed = seed;
    if (seed == 0) // using the time to generate a seed when seed = 0
        srand((unsigned int)time(nullptr));
    else
        srand((unsigned int) this->seed);*/

    unsigned seed1 = time(0);
    // cout << seed1 << endl;
    srand(seed1);

    //getFiles();

    // Setting the method parameters
    setMethodParams();

    // setAuxiliaryParams();

}

Parameters::~Parameters() {

    for (int i = 0; i < numJobs + 1; i++) {
        delete[] C[i];
        delete[] F[i];
        delete[] P[i];
    }

    delete[] C;
    delete[] F;
    delete[] P;

}