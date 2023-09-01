/**
 * @file Genetic.cpp
 *
 * Implements Genetic class methods
 *
 * @author Jordana Mecler
 *
 * Contact: jmecler@inf.puc-rio.br
 *
 */

#include "Genetic.h"
#include "Util.h"

// Main code of the HGA
void Genetic::evolve(int maxIterWithoutImprov) {

    // Individuals used for crossover
    Individual *parent1;
    Individual* parent2;

    nbIterWithoutImprov = 1;
    int nbIterWithoutImprovDiv = 1;
    nbIter = 1;

    string temp;
    int place;
    clock_t debut = clock();

    // Child reference
    offspring = new Individual(parameters);
    // Individual used for local search
    trainer = new Individual(parameters);
    trainer->localSearch = new LocalSearch(parameters, trainer);

    int cont = 1;
    bool controle = true;
    // while (nbIterWithoutImprov < maxIterWithoutImprov) {
    while (cpuTime() - parameters->cpuTime < 100 * parameters->numJobs * parameters->numTools / 1000) {     
        cont++;
        parent1 = population->getIndividualBinT(); // Pick individual by binary tournament

        RR(parent1); //ruin-and-recreate
        
        // LOCAL SEARCH
        trainer->recopyIndividual(trainer, offspring);
        trainer->localSearch->runSearchTotal();
        offspring->recopyIndividual(offspring, trainer);

        cout << population->getBestIndividual()->solutionCost.evaluation << endl;
        cout << cpuTime() - parameters->cpuTime << endl << endl;

        // Tries to add child to population
        place = population->addIndividual(offspring);
        if (place == -2) {
            return;
        }

        // if (place == 0) { // A new best solution has been found
        //     nbIterWithoutImprov = 1;
        //     nbIterWithoutImprovDiv = 1;
        // }
        // else {
        //     nbIterWithoutImprov++;
        // }

        // nbIterWithoutImprovDiv++;
        // nbIter++;

        // DIVERSIFICATION
        // Max iterations without improvement resulting in diversification reached
        // if (nbIterWithoutImprovDiv == parameters->maxDiversify) {
        //     population->diversify();
        //     if (parameters->terminate) {
        //         return;
        //     }
        //     nbIterWithoutImprovDiv = 1;
        // }

    }
    
    cout << cont << endl;
    getchar();
    // cout << endl << endl << cont << endl << endl;
    parameters->nbIter = (unsigned int) nbIter;
    // cout << nbIter << endl;
    
}

void Genetic::crossover3(Individual * parent1) {

    offspring->chromosome = parent1->chromosome;
    offspring->E = parent1->E;
    offspring->Q = parent1->Q;

    int jobsMovidos = 8;
    //const int minJobs = 8;
    //const int maxJobs = parameters->numJobs / 4;
    //int jobsMovidos = minJobs + rand () % (maxJobs - 1);
    
    int menorJobMovido, menorAux;
    int maiorJobMovido, maiorAux;
    vector <int> LJ;
    vector <int> auxiliar(parameters->numJobs);
    
    for (int i = 0; i < offspring->chromosome.size(); i++) {
        auxiliar[offspring->chromosome[i]] = i;
    }

    menorJobMovido = 9999;
    maiorJobMovido = -9999;

    for (int i = 1; i <= jobsMovidos; i++) {
        int seleciona = rand() % offspring->chromosome.size();

        LJ.push_back(offspring->chromosome[seleciona]);

        menorAux = maiorAux = auxiliar[offspring->chromosome[seleciona]];

        if (menorAux < menorJobMovido) {
            menorJobMovido = menorAux;
        }

        if (maiorAux > maiorJobMovido) {
            maiorJobMovido = maiorAux;
        }
       
        offspring->chromosome.erase(offspring->chromosome.begin() + seleciona);
    }

    int p;
    bool improved = false;

    for(int i = menorJobMovido + 1;i < offspring->chromosome.size() + 1; i++) {
        offspring->E[i][0] = offspring->E[i-1][1];
        for(int j = 1;j < parameters->numTools+1; j++) {
            if (j != parameters->numTools) {
                offspring->E[i][j] = max(offspring->E[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1], offspring->E[i-1][j+1]);
            } else {
                offspring->E[i][j] = offspring->E[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1];
            }
        }
        // sol.jobsInfo[offspring->chromosome[i-1]].posicao = i - 1;
    }
    
    int jobInicio = offspring->chromosome.size() - 1 - (parameters->numJobs - maiorJobMovido - 1);
    int maquinaInicio, aux;
    
    aux = jobInicio;
    for(int i = offspring->chromosome.size() - aux;i < offspring->chromosome.size() + 1; i++) {
        offspring->Q[i][0] = offspring->Q[i-1][1];
        maquinaInicio = parameters->numTools-1;
        for(int j = 1;j < parameters->numTools+1; j++) {
            if (j != parameters->numTools) {
                offspring->Q[i][j] = max(offspring->Q[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[jobInicio]][maquinaInicio], offspring->Q[i-1][j+1]);
                maquinaInicio--;
            } else {
                offspring->Q[i][j] = offspring->Q[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[jobInicio]][maquinaInicio];
                maquinaInicio--;
            }
        }
        jobInicio--;  
    }

    int ** C = new int * [parameters->numJobs + 1];
    for (int i = 0; i < parameters->numJobs + 1; i++) {
         C[i] = new int [parameters->numTools + 1];
    }

    int maximo = -9999;
    int vAux, guardaValor;
    while (!LJ.empty()) {
        int mkspPtb = 99999;
        int k = rand() % LJ.size();
        
        for (int i = 0; i < offspring->chromosome.size() + 1; i++) {
                vAux = parameters->numTools;

                for (int j = 1; j <= 1; j++) {
                    guardaValor = max(offspring->Q[offspring->chromosome.size() - i][j] + parameters->jobsToolsMatrix[LJ[k]][parameters->numTools - j], offspring->Q[offspring->chromosome.size()- i][j + 1]);
                    C[i][j] = offspring->E[i][vAux] + max(offspring->Q[offspring->chromosome.size() - i][j] + parameters->jobsToolsMatrix[LJ[k]][parameters->numTools - j], offspring->Q[offspring->chromosome.size() - i][j + 1]);
                
                    if (C[i][j] > maximo) {
                        maximo = C[i][j];
                    }

                    vAux--;
                
                }

                for (int j = 2; j < parameters->numTools; j++) {
                    C[i][j] = offspring->E[i][vAux] + max(guardaValor + parameters->jobsToolsMatrix[LJ[k]][parameters->numTools - j], offspring->Q[offspring->chromosome.size() - i][j + 1]);
                    guardaValor = max(guardaValor + parameters->jobsToolsMatrix[LJ[k]][parameters->numTools - j], offspring->Q[offspring->chromosome.size() - i][j + 1]);

                    if (C[i][j] > maximo) {
                        maximo = C[i][j];
                    }

                    vAux--;
                
                }

                for (int j = parameters->numTools; j <= parameters->numTools; j++) {
                    C[i][j] = offspring->E[i][vAux] + guardaValor + parameters->jobsToolsMatrix[LJ[k]][parameters->numTools - j];
                
                    if (C[i][j] > maximo) {
                        maximo = C[i][j];
                    }

                    vAux--;

                }
                
                if (maximo < mkspPtb) {
                    improved = true;
                    mkspPtb = maximo;
                    p = i;
                }
                maximo = -9999;
        }
    
        offspring->chromosome.insert(offspring->chromosome.begin() + p, LJ[k]);
        LJ.erase(LJ.begin() + k);

        for(int i = p + 1; i < offspring->chromosome.size() + 1; i++) {
            offspring->E[i][0] = offspring->E[i-1][1];
            // offspring->solutionCost.zeroBlocks = offspring->solutionCost.zeroBlocks - offspring->idleOrBlocked[i-1];
            // offspring->idleOrBlocked[i-1] = 0;
            for(int j = 1;j < parameters->numTools+1; j++) {
                if (j != parameters->numTools) {
                    offspring->E[i][j] = max(offspring->E[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1], offspring->E[i-1][j+1]);
                } else {
                    offspring->E[i][j] = offspring->E[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1];
                }

                // offspring->idleOrBlocked[i-1] = offspring->idleOrBlocked[i-1] + offspring->E[i][j] - offspring->E[i-1][j] - parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1];
            }
            // sol.jobsInfo[offspring->chromosome[i-1]].posicao = i-1;
            // offspring->solutionCost.zeroBlocks = offspring->solutionCost.zeroBlocks + offspring->idleOrBlocked[i-1];
            // sol.jobsInfo[i-1].tempo_sistema = offspring->E[i][parameters->numTools] - offspring->E[i][1];
        }
        
        jobInicio = p;
        for(int i = offspring->chromosome.size() - p;i < offspring->chromosome.size() + 1; i++) {
            offspring->Q[i][0] = offspring->Q[i-1][1];
            maquinaInicio = parameters->numTools-1;
            for(int j = 1;j < parameters->numTools+1; j++) {
                if (j != parameters->numTools) {
                    offspring->Q[i][j] = max(offspring->Q[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[jobInicio]][maquinaInicio], offspring->Q[i-1][j+1]);
                    maquinaInicio--;
                } else {
                    offspring->Q[i][j] = offspring->Q[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[jobInicio]][maquinaInicio];
                    maquinaInicio--;
                }
            }
            jobInicio--;  
        }

    }

    offspring->solutionCost.evaluation = offspring->E[parameters->numJobs][parameters->numTools];

    for (int i = 0; i < parameters->numJobs+1; i++) {
        delete[] C[i];
    }

    delete[] C;

    if (offspring->solutionCost.evaluation != offspring->E[parameters->numJobs][parameters->numTools] || offspring->solutionCost.evaluation != offspring->Q[parameters->numJobs][parameters->numTools]) {
        cout << "Problema na perturbação!!";
        exit(1);
    }

    offspring->caminhoCritico();
    offspring->jobsInfo.resize(parameters->numJobs);

}

void Genetic::RR(Individual * parent1) {

    offspring->chromosome = parent1->chromosome;
    offspring->E = parent1->E;
    offspring->Q = parent1->Q;

    //const int minJobs = 8;
    //const int maxJobs = parameters->numJobs / 4;
    // int jobsMovidos = 5 + rand() % (12 - 4);
    // parameters->jobsMovidos = 8 + rand() % (12 - 7);

    vector <int> jobsInseridos;
    vector <int> jobsPos;
    vector <int> :: iterator it;

    // int ** C = new int * [parameters->numJobs + 1];
    // for (int i = 0; i < parameters->numJobs + 1; i++) {
    //      C[i] = new int [parameters->numTools + 1];
    // }

    int cont = 0;

    while (cont <= parameters->jobsMovidos) {
        int escolha = rand() % offspring->chromosome.size();
        int job = offspring->chromosome[escolha];

        it = find(jobsInseridos.begin(), jobsInseridos.end(), job);

        while(it != jobsInseridos.end()) {
            escolha = rand() % offspring->chromosome.size();
            job = offspring->chromosome[escolha];

            it = find(jobsInseridos.begin(), jobsInseridos.end(), job);
        }

        jobsInseridos.push_back(job);
        jobsPos.push_back(escolha);

        cont++;

    }

    int itr = 0;
    
    while(itr < parameters->jobsMovidos) {

        int escolha = jobsPos[itr];
        int job = jobsInseridos[itr];

        int elimina = find(offspring->chromosome.begin(), offspring->chromosome.end(), job) - offspring->chromosome.begin();

        offspring->chromosome.erase(offspring->chromosome.begin() + elimina);

        int p;
        bool improved = false;

        for(int i = elimina + 1;i < offspring->chromosome.size() + 1; i++) {
            offspring->solutionCost.zeroBlocks = offspring->solutionCost.zeroBlocks - offspring->idleOrBlocked[i-1];
            offspring->idleOrBlocked[i-1] = 0;        
        
            offspring->E[i][0] = offspring->E[i-1][1];
            for(int j = 1;j < parameters->numTools+1; j++) {
                if (j != parameters->numTools) {
                    offspring->E[i][j] = max(offspring->E[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1], offspring->E[i-1][j+1]);
                } else {
                    offspring->E[i][j] = offspring->E[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1];
                }

                offspring->idleOrBlocked[i-1] = offspring->idleOrBlocked[i-1] + offspring->E[i][j] - offspring->E[i-1][j] - parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1];
            }
            offspring->solutionCost.zeroBlocks = offspring->solutionCost.zeroBlocks + offspring->idleOrBlocked[i-1];
            // sol.jobsInfo[i-1].posicao = i - 1;
        }
    
        int jobInicio = offspring->chromosome.size() - 1 - (parameters->numJobs - elimina - 1);
        int maquinaInicio, aux;
    
        aux = jobInicio;
        for(int i = offspring->chromosome.size() - aux;i < offspring->chromosome.size() + 1; i++) {
            offspring->Q[i][0] = offspring->Q[i-1][1];
            maquinaInicio = parameters->numTools-1;
            for(int j = 1;j < parameters->numTools+1; j++) {
                if (j != parameters->numTools) {
                    offspring->Q[i][j] = max(offspring->Q[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[jobInicio]][maquinaInicio], offspring->Q[i-1][j+1]);
                    maquinaInicio--;
                } else {
                    offspring->Q[i][j] = offspring->Q[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[jobInicio]][maquinaInicio];
                    maquinaInicio--;
                }
            }
            jobInicio--;  
        }

        int mkspPtb = 99999;
        
        int maximo = -9999;
        int vAux, guardaValor;
        for (int i = 0; i < offspring->chromosome.size() + 1; i++) {
            if (i == escolha) {
                continue;
            }

            vAux = parameters->numTools;

            for (int j = 1; j <= 1; j++) {
                guardaValor = max(offspring->Q[offspring->chromosome.size() - i][j] + parameters->jobsToolsMatrix[job][parameters->numTools - j], offspring->Q[offspring->chromosome.size() - i][j + 1]);
                parameters->C[i][j] = offspring->E[i][vAux] + guardaValor;
                
                if (parameters->C[i][j] > maximo) {
                    maximo = parameters->C[i][j];
                }

                vAux--;
                
            }

            for (int j = 2; j < parameters->numTools; j++) {
                guardaValor = max(guardaValor + parameters->jobsToolsMatrix[job][parameters->numTools - j], offspring->Q[offspring->chromosome.size() - i][j + 1]);
                parameters->C[i][j] = offspring->E[i][vAux] + guardaValor;
                    
                if (parameters->C[i][j] > maximo) {
                    maximo = parameters->C[i][j];
                }

                vAux--;
                
            }

            for (int j = parameters->numTools; j <= parameters->numTools; j++) {
                parameters->C[i][j] = offspring->E[i][vAux] + guardaValor + parameters->jobsToolsMatrix[job][parameters->numTools - j];
                
                if (parameters->C[i][j] > maximo) {
                    maximo = parameters->C[i][j];
                }

                vAux--;

            }
                
            if (maximo < mkspPtb) {
                improved = true;
                mkspPtb = maximo;
                p = i;
            }
            maximo = -9999;
        }

        // cout << job << " " << p << " " << endl;
        offspring->chromosome.insert(offspring->chromosome.begin() + p, job);

        for(int i = p + 1; i < offspring->chromosome.size() + 1; i++) {
            offspring->E[i][0] = offspring->E[i-1][1];
            offspring->solutionCost.zeroBlocks = offspring->solutionCost.zeroBlocks - offspring->idleOrBlocked[i-1];
            offspring->idleOrBlocked[i-1] = 0;
            for(int j = 1;j < parameters->numTools+1; j++) {
                if (j != parameters->numTools) {
                    offspring->E[i][j] = max(offspring->E[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1], offspring->E[i-1][j+1]);
                } else {
                    offspring->E[i][j] = offspring->E[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1];
                }

                offspring->idleOrBlocked[i-1] = offspring->idleOrBlocked[i-1] + offspring->E[i][j] - offspring->E[i-1][j] - parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1];
            }
            // sol.jobsInfo[i-1].posicao = i-1;
            offspring->solutionCost.zeroBlocks = offspring->solutionCost.zeroBlocks + offspring->idleOrBlocked[i-1];
            // sol.jobsInfo[i-1].tempo_sistema = offspring->E[i][parameters->numTools] - offspring->E[i][1];
        }
        
        jobInicio = p;
        for(int i = offspring->chromosome.size() - p;i < offspring->chromosome.size() + 1; i++) {
            offspring->Q[i][0] = offspring->Q[i-1][1];
            maquinaInicio = parameters->numTools-1;
            for(int j = 1;j < parameters->numTools+1; j++) {
                if (j != parameters->numTools) {
                    offspring->Q[i][j] = max(offspring->Q[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[jobInicio]][maquinaInicio], offspring->Q[i-1][j+1]);
                    maquinaInicio--;
                } else {
                    offspring->Q[i][j] = offspring->Q[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[jobInicio]][maquinaInicio];
                    maquinaInicio--;
                }
            }
            jobInicio--;  
        }

        itr++;
    }

    // cout << endl;

    offspring->solutionCost.evaluation = offspring->E[parameters->numJobs][parameters->numTools];

    // for (int i = 0; i < parameters->numJobs+1; i++) {
    //     delete[] C[i];
    // }

    // delete[] C;

    if (offspring->solutionCost.evaluation != offspring->E[parameters->numJobs][parameters->numTools] || offspring->solutionCost.evaluation != offspring->Q[parameters->numJobs][parameters->numTools]) {
        cout << "Problema na perturbação!!";
        exit(1);
    }

    offspring->caminhoCritico();
    offspring->jobsInfo.resize(parameters->numJobs);

}

void Genetic::greedy_insertion_adjacent_swap(Individual * parent1) {

    offspring->chromosome = parent1->chromosome;
    offspring->E = parent1->E;
    offspring->Q = parent1->Q;

    int jobsMovidos = 5;
    //const int minJobs = 8;
    //const int maxJobs = parameters->numJobs / 4;
    //int jobsMovidos = minJobs + rand () % (maxJobs - 1);
    
    int menorJobMovido, menorAux;
    int maiorJobMovido, maiorAux;
    vector <int> LJ;
    vector <int> auxiliar(parameters->numJobs);
    
    for (int i = 0; i < offspring->chromosome.size(); i++) {
        auxiliar[offspring->chromosome[i]] = i;
    }

    menorJobMovido = 9999;
    maiorJobMovido = -9999;

    for (int i = 1; i <= jobsMovidos; i++) {
        int seleciona = rand() % offspring->chromosome.size();

        LJ.push_back(offspring->chromosome[seleciona]);

        menorAux = maiorAux = auxiliar[offspring->chromosome[seleciona]];

        if (menorAux < menorJobMovido) {
            menorJobMovido = menorAux;
        }

        if (maiorAux > maiorJobMovido) {
            maiorJobMovido = maiorAux;
        }
       
        offspring->chromosome.erase(offspring->chromosome.begin() + seleciona);
    }

    int p;
    bool improved = false;

    for(int i = menorJobMovido + 1;i < offspring->chromosome.size() + 1; i++) {
        offspring->solutionCost.zeroBlocks = offspring->solutionCost.zeroBlocks - offspring->idleOrBlocked[i-1];
        offspring->idleOrBlocked[i-1] = 0;    

        offspring->E[i][0] = offspring->E[i-1][1];
        for(int j = 1;j < parameters->numTools+1; j++) {
            if (j != parameters->numTools) {
                offspring->E[i][j] = max(offspring->E[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1], offspring->E[i-1][j+1]);
            } else {
                offspring->E[i][j] = offspring->E[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1];
            }
            offspring->idleOrBlocked[i-1] = offspring->idleOrBlocked[i-1] + offspring->E[i][j] - offspring->E[i-1][j] - parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1];
        }
        offspring->solutionCost.zeroBlocks = offspring->solutionCost.zeroBlocks + offspring->idleOrBlocked[i-1];
    }
    
    int jobInicio = offspring->chromosome.size() - 1 - (parameters->numJobs - maiorJobMovido - 1);
    int maquinaInicio, aux;
    
    aux = jobInicio;
    for(int i = offspring->chromosome.size() - aux;i < offspring->chromosome.size() + 1; i++) {
        offspring->Q[i][0] = offspring->Q[i-1][1];
        maquinaInicio = parameters->numTools-1;
        for(int j = 1;j < parameters->numTools+1; j++) {
            if (j != parameters->numTools) {
                offspring->Q[i][j] = max(offspring->Q[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[jobInicio]][maquinaInicio], offspring->Q[i-1][j+1]);
                maquinaInicio--;
            } else {
                offspring->Q[i][j] = offspring->Q[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[jobInicio]][maquinaInicio];
                maquinaInicio--;
            }
        }
        jobInicio--;  
    }

    int ** C = new int * [parameters->numJobs + 1];
    for (int i = 0; i < parameters->numJobs + 1; i++) {
         C[i] = new int [parameters->numTools + 1];
    }

    int auxx;

    int maximo = -9999;
    int vAux, guardaValor;
    while (!LJ.empty()) {
        int mkspPtb = 99999;
        int k = rand() % LJ.size();
        
        for (int i = 0; i < offspring->chromosome.size() + 1; i++) {
                vAux = parameters->numTools;

                for (int j = 1; j <= 1; j++) {
                    guardaValor = max(offspring->Q[offspring->chromosome.size() - i][j] + parameters->jobsToolsMatrix[LJ[k]][parameters->numTools - j], offspring->Q[offspring->chromosome.size()- i][j + 1]);
                    C[i][j] = offspring->E[i][vAux] + guardaValor;
                
                    if (C[i][j] > maximo) {
                        maximo = C[i][j];
                    }

                    vAux--;
                
                }

                for (int j = 2; j < parameters->numTools; j++) {
                    guardaValor = max(guardaValor + parameters->jobsToolsMatrix[LJ[k]][parameters->numTools - j], offspring->Q[offspring->chromosome.size() - i][j + 1]);
                    C[i][j] = offspring->E[i][vAux] + guardaValor;
                    

                    if (C[i][j] > maximo) {
                        maximo = C[i][j];
                    }

                    vAux--;
                
                }

                for (int j = parameters->numTools; j <= parameters->numTools; j++) {
                    C[i][j] = offspring->E[i][vAux] + guardaValor + parameters->jobsToolsMatrix[LJ[k]][parameters->numTools - j];
                
                    if (C[i][j] > maximo) {
                        maximo = C[i][j];
                    }

                    vAux--;

                }
                
                if (maximo < mkspPtb) {
                    improved = true;
                    mkspPtb = maximo;
                    p = i;
                }
                maximo = -9999;
        }
    
        offspring->chromosome.insert(offspring->chromosome.begin() + p, LJ[k]);
        LJ.erase(LJ.begin() + k);

        for (int i = p + 1; i < offspring->chromosome.size() - 1; i++) {
            swap(offspring->chromosome[i], offspring->chromosome[i + 1]);
        }

        if (p + 1 != offspring->chromosome.size() - 1) {
            auxx = offspring->chromosome.size() - 1;
        } else {
            auxx = p;
        }

        for(int i = p + 1; i < offspring->chromosome.size() + 1; i++) {
            offspring->E[i][0] = offspring->E[i-1][1];
            offspring->solutionCost.zeroBlocks = offspring->solutionCost.zeroBlocks - offspring->idleOrBlocked[i-1];
            offspring->idleOrBlocked[i-1] = 0;
            for(int j = 1;j < parameters->numTools+1; j++) {
                if (j != parameters->numTools) {
                    offspring->E[i][j] = max(offspring->E[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1], offspring->E[i-1][j+1]);
                } else {
                    offspring->E[i][j] = offspring->E[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1];
                }

                offspring->idleOrBlocked[i-1] = offspring->idleOrBlocked[i-1] + offspring->E[i][j] - offspring->E[i-1][j] - parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1];
            }
            offspring->solutionCost.zeroBlocks = offspring->solutionCost.zeroBlocks + offspring->idleOrBlocked[i-1];
        }
        
        jobInicio = auxx;
        for(int i = offspring->chromosome.size() - jobInicio;i < offspring->chromosome.size() + 1; i++) {
            offspring->Q[i][0] = offspring->Q[i-1][1];
            maquinaInicio = parameters->numTools-1;
            for(int j = 1;j < parameters->numTools+1; j++) {
                if (j != parameters->numTools) {
                    offspring->Q[i][j] = max(offspring->Q[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[jobInicio]][maquinaInicio], offspring->Q[i-1][j+1]);
                    maquinaInicio--;
                } else {
                    offspring->Q[i][j] = offspring->Q[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[jobInicio]][maquinaInicio];
                    maquinaInicio--;
                }
            }
            jobInicio--;  
        }

    }

    offspring->solutionCost.evaluation = offspring->E[parameters->numJobs][parameters->numTools];

    for (int i = 0; i < parameters->numJobs+1; i++) {
        delete[] C[i];
    }

    delete[] C;

    if (offspring->solutionCost.evaluation != offspring->E[parameters->numJobs][parameters->numTools] || offspring->solutionCost.evaluation != offspring->Q[parameters->numJobs][parameters->numTools]) {
        cout << "Problema na perturbação!!";
        exit(1);
    }

}

void Genetic::crossover2(Individual * parent1, Individual * parent2) {

    // Beginning and end of the crossover zone
    unsigned int begin = rand() % parameters->numJobs;
    unsigned int end = rand() % parameters->numJobs;
    while (end == begin && parameters->numJobs > 1)
        end = rand() % parameters->numJobs;

    if (begin > end) {
        unsigned int temp = begin;
        begin = end;
        end = temp;
    }

    vector <int> LC;

    // Copy part of parent1 to child
    // cout << begin << " " << end << endl << endl;
    offspring->chromosome.clear();
    for (unsigned int i = 0; i < parameters->numJobs; i++) {
        if (i >= begin && i <= end) {
            offspring->chromosome.push_back(parent1->chromosome[i]);
        } else {
            LC.push_back(parent1->chromosome[i]);
        }
    }

    for (int i = 0; i < offspring->chromosome.size(); i++) {
        cout << offspring->chromosome[i] << " ";
    }
    cout << endl;

    for (int i = 0; i < LC.size(); i++) {
        cout << LC[i] << " ";
    }
    cout << endl << endl;

    cout << LC.size() << endl << endl;

    getchar();

    for(int i = 0; i < parameters->numJobs + 1; i++) {
        offspring->E[i][0] = 0;
        offspring->Q[i][0] = 0;
    }

    for(int j = 0; j < parameters->numTools + 1; j++) {
        offspring->E[0][j] = 0;
        offspring->Q[0][j] = 0;
    }
    
    offspring->solutionCost.zeroBlocks = 0; 

    for (int i = 1; i < offspring->chromosome.size() + 1; i++) {
        offspring->E[i][0] = offspring->E[i-1][1];
        offspring->idleOrBlocked.push_back(0);
        for (int j = 1; j < parameters->numTools + 1; j++) {
            if (j != parameters->numTools) {
                offspring->E[i][j] = max(offspring->E[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1], offspring->E[i-1][j+1]);
            } else {
                offspring->E[i][j] = offspring->E[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1];
            }
            offspring->idleOrBlocked[i-1] += offspring->E[i][j] - offspring->E[i-1][j] - parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1];
        }
        offspring->solutionCost.zeroBlocks += offspring->idleOrBlocked[i-1];
    }

    int jobInicio = offspring->chromosome.size() - 1;
    int maquinaInicio;
    for(int i=1;i<offspring->chromosome.size()+1;i++) {
        offspring->Q[i][0] = offspring->Q[i-1][1];
        maquinaInicio = parameters->numTools - 1;
        for(int j=1;j<parameters->numTools+1;j++) {
            
            if (j != parameters->numTools) {
                offspring->Q[i][j] = max(offspring->Q[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[jobInicio]][maquinaInicio], offspring->Q[i-1][j+1]);
                maquinaInicio--;
            } else {
                offspring->Q[i][j] = offspring->Q[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[jobInicio]][maquinaInicio];
                maquinaInicio--;
            }
        }
        jobInicio--;  
    }

    int maximo = -9999;
    int vAux, guardaValor, p;
    while (!LC.empty()) {
        int mkspPtb = 99999;
        
        for (int i = 0; i < offspring->chromosome.size() + 1; i++) {
                vAux = parameters->numTools;

                for (int j = 1; j <= 1; j++) {
                    guardaValor = max(offspring->Q[offspring->chromosome.size() - i][j] + parameters->jobsToolsMatrix[LC[0]][parameters->numTools - j], offspring->Q[offspring->chromosome.size()- i][j + 1]);
                    parameters->C[i][j] = offspring->E[i][vAux] + max(offspring->Q[offspring->chromosome.size() - i][j] + parameters->jobsToolsMatrix[LC[0]][parameters->numTools - j], offspring->Q[offspring->chromosome.size() - i][j + 1]);
                
                    if (parameters->C[i][j] > maximo) {
                        maximo = parameters->C[i][j];
                    }

                    vAux--;
                
                }

                for (int j = 2; j < parameters->numTools; j++) {
                    parameters->C[i][j] = offspring->E[i][vAux] + max(guardaValor + parameters->jobsToolsMatrix[LC[0]][parameters->numTools - j], offspring->Q[offspring->chromosome.size() - i][j + 1]);
                    guardaValor = max(guardaValor + parameters->jobsToolsMatrix[LC[0]][parameters->numTools - j], offspring->Q[offspring->chromosome.size() - i][j + 1]);

                    if (parameters->C[i][j] > maximo) {
                        maximo = parameters->C[i][j];
                    }

                    vAux--;
                
                }

                for (int j = parameters->numTools; j <= parameters->numTools; j++) {
                    parameters->C[i][j] = offspring->E[i][vAux] + guardaValor + parameters->jobsToolsMatrix[LC[0]][parameters->numTools - j];
                
                    if (parameters->C[i][j] > maximo) {
                        maximo = parameters->C[i][j];
                    }

                    vAux--;

                }
                
                if (maximo < mkspPtb) {
                    mkspPtb = maximo;
                    p = i;
                }
                maximo = -9999;
        }

        offspring->chromosome.insert(offspring->chromosome.begin() + p, LC[0]);
        LC.erase(LC.begin());

        for(int i = p + 1; i < offspring->chromosome.size() + 1; i++) {
            offspring->E[i][0] = offspring->E[i-1][1];
            offspring->solutionCost.zeroBlocks= offspring->solutionCost.zeroBlocks - offspring->idleOrBlocked[i-1];
            offspring->idleOrBlocked[i-1] = 0;
            for(int j = 1;j < parameters->numTools+1; j++) {
                if (j != parameters->numTools) {
                    offspring->E[i][j] = max(offspring->E[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1], offspring->E[i-1][j+1]);
                } else {
                    offspring->E[i][j] = offspring->E[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1];
                }

                offspring->idleOrBlocked[i-1] = offspring->idleOrBlocked[i-1] + offspring->E[i][j] - offspring->E[i-1][j] - parameters->jobsToolsMatrix[offspring->chromosome[i-1]][j-1];
            }
            offspring->solutionCost.zeroBlocks= offspring->solutionCost.zeroBlocks + offspring->idleOrBlocked[i-1];
        }

        jobInicio = p;
        for(int i = offspring->chromosome.size() - p;i < offspring->chromosome.size() + 1; i++) {
            offspring->Q[i][0] = offspring->Q[i-1][1];
            maquinaInicio = parameters->numTools-1;
            for(int j = 1;j < parameters->numTools+1; j++) {
                if (j != parameters->numTools) {
                    offspring->Q[i][j] = max(offspring->Q[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[jobInicio]][maquinaInicio], offspring->Q[i-1][j+1]);
                    maquinaInicio--;
                } else {
                    offspring->Q[i][j] = offspring->Q[i][j-1] + parameters->jobsToolsMatrix[offspring->chromosome[jobInicio]][maquinaInicio];
                    maquinaInicio--;
                }
            }
            jobInicio--;  
        }

    }

    offspring->solutionCost.evaluation = offspring->E[parameters->numJobs][parameters->numTools];

    if (offspring->E[parameters->numJobs][parameters->numTools] != offspring->Q[parameters->numJobs][parameters->numTools]) {
        cout << "problem crossover";
        exit(1);
    }

    // cout << offspring->E[offspring->chromosome.size()][parameters->numTools] << " " << offspring->Q[offspring->chromosome.size()][parameters->numTools];
    // exit(1);

    offspring->jobsInfo.resize(parameters->numJobs);
    offspring->caminhoCritico();

}

void Genetic::crossoverOX(Individual *parent1, Individual *parent2) {
    
    // Beginning and end of the crossover zone
    unsigned int begin = rand() % parameters->numJobs;
    unsigned int end = rand() % parameters->numJobs;
    while (end == begin && parameters->numJobs > 1)
        end = rand() % parameters->numJobs;

    if (begin > end) {
        unsigned int temp = begin;
        begin = end;
        end = temp;
    }

    // cout << begin << " " << end << endl;

    // Copy part of parent1 to child
    offspring->chromosome = parent1->chromosome;
    for (unsigned int i = 0; i < parameters->positionsOffspring.size(); i++) {
        parameters->positionsOffspring[i] = false;
    }
    for (unsigned int i = begin; i <= end; i++) {
        parameters->positionsOffspring[parent1->chromosome[i]] = true;
    }
   
    // Copy unused values of parent2 to child sequentially
    unsigned int pos = end + 1, i = end + 1;
    while (pos < parameters->positionsOffspring.size()) {
        if (!parameters->positionsOffspring[parent2->chromosome[i]]) {
            offspring->chromosome[pos] = parent2->chromosome[i];
            pos++;
        }
        i++;
        if (i == parameters->positionsOffspring.size()) {
            i = 0;
        }
    }
    if (i == parameters->positionsOffspring.size()) {
        i = 0;
    }
    pos = 0;
    while (pos < begin) {
        if (!parameters->positionsOffspring[parent2->chromosome[i]]) {
            offspring->chromosome[pos] = parent2->chromosome[i];
            pos++;
        }
        i++;
        if (i == parameters->positionsOffspring.size()) {
            i = 0;
        }
    }
    offspring->solutionCost.evaluation = offspring->calcCost(-1);
    offspring->jobsInfo.resize(parameters->numJobs);
    offspring->caminhoCritico();
}

void Genetic::crossoverSJOX(Individual *parent1, Individual *parent2) {

    int point = rand() % parameters->numJobs;

    // cout << "Pai 1: ";
    // for (int i = 0; i < parameters->numJobs; i++) {
    //     cout << parent1->chromosome[i] << " ";
    // }
    // cout << endl;
    // cout << "Pai 2: ";
    // for (int i = 0; i < parameters->numJobs; i++) {
    //     cout << parent2->chromosome[i] << " ";
    // }
    // cout << endl << endl;

    // cout << "Ponto: " << point << endl;

    // cout << parameters->positionsOffspring.size() << endl;

    for (unsigned int i = 0; i < parameters->positionsOffspring.size(); i++) {
        parameters->positionsOffspring[i] = false;
    }


    offspring->chromosome = parent1->chromosome;

    //     for (int i = 0; i < parameters->numJobs; i++) {
    //     cout << offspring->chromosome[i] << " ";
    // }
    // cout << endl;

    for (int i = 0; i < parameters->numJobs; i++) {
        if (parent1->chromosome[i] != parent2->chromosome[i]) {
            offspring->chromosome[i] = -1;
        } else {
            parameters->positionsOffspring[parent2->chromosome[i]] = true;
        }
    }

    // for (int i = 0; i < parameters->numJobs; i++) {
    //     cout << offspring->chromosome[i] << " ";
    // }
    // cout << endl;

    for (int i = 0; i <= point; i++) {
        if (offspring->chromosome[i] == -1) {
            offspring->chromosome[i] = parent1->chromosome[i];
            parameters->positionsOffspring[parent1->chromosome[i]] = true;
        }
    }
    //     for (int i = 0; i < parameters->numJobs; i++) {
    //     cout << offspring->chromosome[i] << " ";
    // }

    // cout << endl;


    int pos = point + 1;
    int i = 0;
    while (pos < parameters->numJobs) {
        if (offspring->chromosome[pos] == -1) {
            if (!parameters->positionsOffspring[parent2->chromosome[i]]) {
                offspring->chromosome[pos] = parent2->chromosome[i];
                parameters->positionsOffspring[parent2->chromosome[i]] = true;
                pos++;
            } else {
                i++;
            }
        } else {
            pos++;
        }
    }

    // for (int i = point + 1; i < parameters->numJobs; i++) {
    //     if (!parameters->positionsOffspring[parent2->chromosome[i]]) {
    //         offspring->chromosome[i] = parent2->chromosome[i];
    //         parameters->positionsOffspring[parent2->chromosome[i]] = true;
    //     }
    // }

    //     for (int i = 0; i < parameters->numJobs; i++) {
    //     cout << offspring->chromosome[i] << " ";
    // }
    //                     getchar();

    //         cout << "opa" << endl;
    // getchar();

    // for (int i = 0; i < parameters->numJobs; i++) {
    //     cout << offspring->chromosome[i] << " ";
    // }
    // cout << endl;
    // getchar();

    offspring->solutionCost.evaluation = offspring->calcCost(-1);
    offspring->jobsInfo.resize(parameters->numJobs);
    offspring->caminhoCritico();

    // cout << "oba" << endl;

    //         cout << "opa2" << endl;
    // getchar();

}

Genetic::Genetic(Parameters *parameters, Population *population, clock_t ticks, bool traces) :
        parameters(parameters), population(population), ticks(ticks), traces(traces) {
            
}

Genetic::~Genetic() {
  
    /*delete offspring;
    delete trainer->localSearch;
    delete trainer;*/
}
