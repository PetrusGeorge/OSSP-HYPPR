/**
 * @file LocalSearch.cpp
 *
 * Implements LocalSearch class methods
 *
 * @author Jordana Mecler
 *
 * Contact: jmecler@inf.puc-rio.br
 *
 */

#include "LocalSearch.h"
#include "Individual.h"
#include "Util.h"

bool LocalSearch::searchNeighborhood(unsigned int i) {
       
    if (i == 1) {
        return relocate();
    } else if (i == 2) {
        return twoOpt(i);
    } else if (i == 3) {    
        return twoOpt(i);
    } else if (i == 4) {
        return twoOpt(i);
    } else { 
        return swap();
    }
    
}

void LocalSearch::constructionSearch() {
    // For each local search procedure, the order of
    // the search is shuffled and it stops when cant
    // find a better solution using this procedure
    // in an entire loop in the indices

    bool foundBetter;

    while (true) {
        foundBetter = false;
        // shuffleIndices();
        // shuffleMoves();
        for (unsigned int i = 1; i <= 3 && !foundBetter; i++) {
            foundBetter = searchNeighborhood(i);
        }
        if (!foundBetter) break;
    }
}


void LocalSearch::runSearchTotal() {

    bool foundBetter;
    double r;

    while (true) {
        foundBetter = false;
        // shuffleIndices();
        // shuffleMoves();
        for (unsigned int i = 1; i <= parameters->nViz && !foundBetter; i++) {
            foundBetter = searchNeighborhood(i);

            // if (i == 4 && !foundBetter && rand() / (RAND_MAX + 1.) < 0.01) {
            //     foundBetter = swap();
            // }
        }
        if (!foundBetter) break;
    }

    // foundBetter = swap();

}

bool LocalSearch::relocate() {

    // parameters->viz1.second++;

    // double startTime = cpuTime();

    bool improved = false;
    int lowerBound;
    int teste;

    // int ** C = new int * [parameters->numJobs + 1];
    // for (int i = 0; i < parameters->numJobs + 1; i++) {
    //      C[i] = new int [parameters->numTools + 1];
    // }

    // for (const vector<unsigned int> &pos : parameters->orderLSi) {
    for (int k = 0; k < parameters->numJobs; k++) {

        preInsertion(parameters->F, parameters->P, k, 1);

        //Verifica a maximo posição para reinserir a k-ésima tarefa
        int maximo = -9999;
        int vAux, guardaValor;
        for (int l = 0; l < parameters->numJobs; l++) {
            if (!(k - 1 <= l && l < k + 1)) {

                if (individual->jobsInfo[k].b == individual->jobsInfo[l].b) {  
                    continue;
                }
                
                lowerBound = getInsertionLowerBound(k, l);
                if (lowerBound >= 0) {
                    continue;
                }

                vAux = parameters->numTools;

                for (int j = 1; j < parameters->numTools + 1; j++) {
                    if (j != 1) {
                        if (j != parameters->numTools) {
                            guardaValor = max(guardaValor + (int)parameters->jobsToolsMatrix[individual->chromosome[k]][parameters->numTools - j], parameters->P[parameters->numJobs - l - 1][j + 1]); 
                            parameters->C[l][j] = parameters->F[l][vAux] + guardaValor;
                        } else {
                            parameters->C[l][j] = parameters->F[l][vAux] + guardaValor + (int)parameters->jobsToolsMatrix[individual->chromosome[k]][parameters->numTools - j];
                        } 
                    } else {
                        guardaValor = max(parameters->P[parameters->numJobs - l - 1][j] + (int)parameters->jobsToolsMatrix[individual->chromosome[k]][parameters->numTools - j], parameters->P[parameters->numJobs - l - 1][j + 1]);
                        parameters->C[l][j] = parameters->F[l][vAux] + guardaValor;
                    }
                    
                    if (parameters->C[l][j] >= maximo) {
                        maximo = parameters->C[l][j];
                    }

                    if (parameters->C[l][j] > individual->solutionCost.evaluation) {
                        break;
                    }

                    vAux--;
                }

                if (maximo <= individual->solutionCost.evaluation) {
                    unsigned int temporaryIdleOrBlocked;
            
                    temporaryIdleOrBlocked = calcIdleOrBlockedRelocate(k, l, 1);

                    // if (maximo == individual->solutionCost.evaluation) {
                    //     parameters->ties++;
                    // }
               
                    if (maximo < individual->solutionCost.evaluation || temporaryIdleOrBlocked < individual->solutionCost.zeroBlocks) {
                        // if (maximo < individual->solutionCost.evaluation ) {
                        //     parameters->improvesPrimary++;
                        // } else {
                        //     parameters->improvesSecondary++;
                        // }

                        improved = true;
                
                        vector<int> aux;
                        for(int i = 0; i < 1; i++){
                            int job = individual->chromosome[k];
                            aux.push_back(job);
                            individual->chromosome.erase(individual->chromosome.begin() + k);
                        }
                        reverse(aux.begin(),aux.end());
                        for(int i = 0; i < aux.size(); i++){
                            if(l < k){
                                individual->chromosome.insert(individual->chromosome.begin() + l, aux[i]);
                            }else{
                                int index = l - 1 + 1;
                                individual->chromosome.insert(individual->chromosome.begin() + index, aux[i]);
                            }
                        } 

                        int q;
                        if (k < l) {
                            q = individual->updateCost(k, l);
                        } else {
                            q = individual->updateCost(l, k);
                        }
                
                        individual->solutionCost.evaluation = maximo;
                        individual->caminhoCritico();
                        preInsertion(parameters->F, parameters->P, k, 1);

                    } 
                }
            }

            maximo = -9999;
        }
    
    }

    // if (individual->solutionCost.evaluation != individual->E[parameters->numJobs][parameters->numTools] || individual->solutionCost.evaluation != individual->Q[parameters->numJobs][parameters->numTools]) {
    //     cout << "Problema na vizinhança insertion!!";
    //     exit(1);
    // }

    // cout << cpuTime() - startTime;
    // getchar();
    // parameters->viz1.first = parameters->viz1.first + (cpuTime() - startTime);

    return improved;


}

bool LocalSearch::swap() {

    // parameters->viz2.second++;

    // double startTime = cpuTime();

    // int ** C = new int * [parameters->numJobs + 1];
    // for (int i = 0; i < parameters->numJobs + 1; i++) {
    //      C[i] = new int [parameters->numTools + 1];
    // }
   
    bool improved = false;
    int i, j, indice;
    int lowerBound;
    int temp, modifiedCost;

    // shuffleMoves();

    // for (const vector<unsigned int> &pos : parameters->orderLSs) {
    for (int i = 0; i < parameters->numJobs; i++) {
        //  i = parameters->listaSwap[v];

        parameters->C[i + 1][0] = individual->E[i][1];
        for (int j = i + 1; j < parameters->numJobs; j++) {

            if (individual->jobsInfo[i].b == individual->jobsInfo[j].b) {
                continue;
            }

            // lowerBound = getSwapLowerBound(i, j);
            // if (lowerBound >= 0) {      
            //     continue;
            // }

            for (int k = i + 1; k < j + 2; k++) {

                indice = k - 1;
                if (indice == i) {
                    indice = j;
                    //parameters->C[k][0] = individual->E[k-1][1];
                } else if (indice == j){
                    indice = i;
                    //parameters->C[k][0] = individual->E[k-1][1];
                }

                if (k != i + 1) {
                    for (int l = 1; l < parameters->numTools; l++) {
                        parameters->C[k][l] = max(parameters->C[k][l-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[indice]][l-1], parameters->C[k-1][l+1]);
                    } 
                } else {
                    for (int l = 1; l < parameters->numTools; l++) {
                        parameters->C[k][l] = max(parameters->C[k][l-1] + parameters->jobsToolsMatrix[individual->chromosome[indice]][l-1], individual->E[k-1][l+1]);
                    }
                }

                parameters->C[k][parameters->numTools] = parameters->C[k][parameters->numTools-1] + parameters->jobsToolsMatrix[individual->chromosome[indice]][parameters->numTools-1];

                if (k + 1 <= parameters->numJobs) {
                    parameters->C[k + 1][0] = parameters->C[k][1];
                }
            }

            if (parameters->C[j + 1][parameters->numTools] > individual->solutionCost.evaluation) {
                continue;
            }
            
            int maximo = 0;
            for (int l = 1; l < parameters->numTools + 1; l++) {
            
                if (parameters->C[j + 1][l] + individual->Q[parameters->numJobs - j - 1][parameters->numTools - l + 1] > maximo) {
                    maximo = parameters->C[j + 1][l] + individual->Q[parameters->numJobs - j - 1][parameters->numTools - l + 1];
                } 

                if (maximo > individual->solutionCost.evaluation) {
                    break;
                }
            }
        
            if (maximo <= individual->solutionCost.evaluation) {
                unsigned int temporaryIdleOrBlocked = calcIdleOrBlockedSwap(i, j);

                // if (maximo == individual->solutionCost.evaluation) {
                //     parameters->ties++;
                // }                
                
                if (maximo < individual->solutionCost.evaluation || temporaryIdleOrBlocked < individual->solutionCost.zeroBlocks) {
                    // if (maximo < individual->solutionCost.evaluation) {
                    //     parameters->improvesPrimary++;
                    // } else {
                    //     parameters->improvesSecondary++;
                    // } 
                
                    improved = true;
                    individual->solutionCost.evaluation = maximo;
            
                    std::swap(individual->chromosome[i], individual->chromosome[j]);
                    int q = individual->updateCost(i, j);
                    individual->caminhoCritico();
                } 
            } 
        
            // if (improved) {
            //     break;
            // }
        
        }
        
        // if (improved) {
        //     break;
        // }
    
    }

    // parameters->viz2.first = parameters->viz2.first + (cpuTime() - startTime);

    // if (individual->solutionCost.evaluation != individual->E[parameters->numJobs][parameters->numTools] || individual->solutionCost.evaluation != individual->Q[parameters->numJobs][parameters->numTools]) {
    //     cout << "Problema na vizinhança swap!!";
    //     exit(1);
    // }

    return improved;
}

bool LocalSearch::twoOpt(int I) {

    // parameters->viz3.second++;

    // double startTime = cpuTime();

    /*cout << "vamos ver " << individual->solutionCost.evaluation << endl;
    getchar();*/


    // int ** C = new int * [parameters->numJobs + 1];
    // for (int i = 0; i < parameters->numJobs + 1; i++) {
    //      C[i] = new int [parameters->numTools + 1];
    // }

    // const int I = 2;
    int ** bora = new int * [I];
    for (int v = 0; v < I; v++) {
        bora[v] = new int[parameters->numTools + 1];
    }

    for(int i = 0; i < parameters->numTools + 1; i++) {
        bora[0][i] = 0;
    }

    for(int i = 0; i < I; i++) {
        bora[i][0] = 0;
    }

    bool improved = false;
    int k, i;
    int auxilia = I - 1;

    for (int k = 0; k < (parameters->numJobs + 1 - I); k++) {
    // for (const vector<unsigned int> &pos : parameters->orderLSii) {

        preInsertion(parameters->F, parameters->P, k, I);

        //Verifica a maximo posição para reinserir a k-ésima tarefa
        int maximo = 0;
        int vAux, guardaValor;

        for (int i = 0; i < parameters->numJobs; i++) {
            if (!(k - 1 <= i && i < k + I)) { 
                if (i > k) {
                    if (i < parameters->numJobs - 1) {
                        vAux = parameters->numTools;

                        for (int l = 1; l <= 1; l++) {
                            bora[l][0] = parameters->F[i-auxilia][1];
                            for (int m = 1; m < parameters->numTools + 1; m++) {
                                if (m != parameters->numTools) {
                                    bora[l][m] = max(bora[l][m-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[k + l - 1]][m-1], parameters->F[i - auxilia][m+1]);
                                } else {
                                    bora[l][m] = bora[l][m-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[k + l - 1]][m-1];
                                }                                
                            }
                        }

                        for (int l = 2; l <= I - 1; l++) {
                            bora[l][0] = bora[l-1][1];
                            for (int m = 1; m < parameters->numTools + 1; m++) {
                                if (m != parameters->numTools) {
                                    bora[l][m] = max(bora[l][m-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[k + l - 1]][m-1], bora[l-1][m+1]);
                                } else {
                                    bora[l][m] = bora[l][m-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[k + l - 1]][m-1];
                                }
                            }
                        }

                        guardaValor = max(parameters->P[parameters->numJobs - i - I + auxilia][1] + (int)parameters->jobsToolsMatrix[individual->chromosome[k + I - 1]][parameters->numTools - 1], parameters->P[parameters->numJobs - i - I + auxilia][2]);
                        parameters->C[i][1] = bora[I - 1][vAux] + guardaValor;
                        
                        maximo = parameters->C[i][1];
                        vAux--;

                        for (int j = 2; maximo <= individual->solutionCost.evaluation && j < parameters->numTools; j++) {
                            guardaValor = max(guardaValor + (int)parameters->jobsToolsMatrix[individual->chromosome[k + I - 1]][parameters->numTools - j], parameters->P[parameters->numJobs - i - I + auxilia][j + 1]);
                            parameters->C[i][j] = bora[I - 1][vAux] + guardaValor;
                            
                            if (parameters->C[i][j] >= maximo) {
                                maximo = parameters->C[i][j];
                            }

                            vAux--;
                        
                        }
                        
                        if (maximo <= individual->solutionCost.evaluation) {
                            parameters->C[i][parameters->numTools] = bora[I - 1][vAux] + guardaValor + (int)parameters->jobsToolsMatrix[individual->chromosome[k + I - 1]][0];

                            if (parameters->C[i][parameters->numTools] >= maximo) {
                                maximo = parameters->C[i][parameters->numTools];
                            }                        
                        }

                    } else {
                        int cont = 0;
                        for (int l = parameters->numJobs - I + 1; l < parameters->numJobs + 1; l++) {
                            parameters->F[l][0] = parameters->F[l-1][1];
                            for (int m = 1; m < parameters->numTools + 1; m++) {
                                if (m != parameters->numTools) {
                                    parameters->F[l][m] = max(parameters->F[l][m-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[k + cont]][m-1], parameters->F[l-1][m+1]);
                                } else {
                                    parameters->F[l][m] = parameters->F[l][m-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[k + cont]][m-1];
                                }
                            }
                            cont++;
                        }

                        maximo = parameters->F[parameters->numJobs][parameters->numTools];
                    }
                } else if (i == 0){
                    int cont = I - 1;
                    for (int l = parameters->numJobs - I + 1; l < parameters->numJobs + 1; l++) {
                        parameters->P[l][0] = parameters->P[l-1][1];
                        int maquinaInicio = parameters->numTools-1;
                        for(int m=1;m<parameters->numTools+1;m++) {
                            if (m != parameters->numTools) {
                                parameters->P[l][m] = max(parameters->P[l][m-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[k + cont]][maquinaInicio], parameters->P[l-1][m+1]);
                                maquinaInicio--;
                            } else {
                                parameters->P[l][m] = parameters->P[l][m-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[k + cont]][maquinaInicio];
                                maquinaInicio--;
                            }
                        }            
                        cont--;
                    }

                    maximo = parameters->P[parameters->numJobs][parameters->numTools];
                } else if (i < k) {
                    vAux = parameters->numTools;
                    
                    for (int l = 1; l <= 1; l++) {
                        bora[l][0] = parameters->F[i][1];
                        for (int m = 1; m < parameters->numTools + 1; m++) {
                            if (m != parameters->numTools) {
                                bora[l][m] = max(bora[l][m-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[k + l - 1]][m-1], parameters->F[i][m+1]);
                            } else {
                                bora[l][m] = bora[l][m-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[k + l - 1]][m-1];
                            }
                        }
                    }
                    
                    for (int l = 2; l <= I - 1; l++) {
                        bora[l][0] = bora[l-1][1];
                        for (int m = 1; m < parameters->numTools + 1; m++) {
                            if (m != parameters->numTools) {
                                bora[l][m] = max(bora[l][m-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[k + l - 1]][m-1], bora[l-1][m+1]);
                            } else {
                                bora[l][m] = bora[l][m-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[k + l - 1]][m-1];
                            }
                        }
                    }

                    guardaValor = max(parameters->P[parameters->numJobs - i - I][1] + (int)parameters->jobsToolsMatrix[individual->chromosome[k + I - 1]][parameters->numTools - 1], parameters->P[parameters->numJobs - i - I][2]);
                    parameters->C[i][1] = bora[I - 1][vAux] + guardaValor;
                    
                    maximo = parameters->C[i][1];
                    vAux--;
                    
                    for (int j = 2; maximo <= individual->solutionCost.evaluation && j < parameters->numTools; j++) {
                        guardaValor = max(guardaValor + (int)parameters->jobsToolsMatrix[individual->chromosome[k + I - 1]][parameters->numTools - j], parameters->P[parameters->numJobs - i - I][j + 1]);
                        parameters->C[i][j] = bora[I - 1][vAux] + guardaValor;
                        
                        if (parameters->C[i][j] >= maximo) {
                            maximo = parameters->C[i][j];
                        }

                        vAux--;
                    }

                    if (maximo <= individual->solutionCost.evaluation) {
                        parameters->C[i][parameters->numTools] = bora[I - 1][vAux] + guardaValor + (int)parameters->jobsToolsMatrix[individual->chromosome[k + I - 1]][0];

                        if (parameters->C[i][parameters->numTools] >= maximo) {
                            maximo = parameters->C[i][parameters->numTools];
                        }
                    }


                }
        
                if (maximo <= individual->solutionCost.evaluation) {
                    unsigned int temporaryIdleOrBlocked;
                    temporaryIdleOrBlocked = calcIdleOrBlockedRelocate(k, i, I);

                    // if (maximo == individual->solutionCost.evaluation) {
                    //     parameters->ties++;
                    // }
               
                    if (maximo < individual->solutionCost.evaluation || temporaryIdleOrBlocked < individual->solutionCost.zeroBlocks) {
                        // if (maximo < individual->solutionCost.evaluation ) {
                        //     parameters->improvesPrimary++;
                        // } else {
                        //     parameters->improvesSecondary++;
                        // }            

                        improved = true;
                        individual->solutionCost.evaluation = maximo;
                
                        vector<int> aux;
                        for(int l = 0; l < I; l++){
                            int job = individual->chromosome[k];
                            aux.push_back(job);
                            individual->chromosome.erase(individual->chromosome.begin() + k);
                        }
                        reverse(aux.begin(),aux.end());
                        for(int l = 0; l < aux.size(); l++){
                            if(i < k){
                                individual->chromosome.insert(individual->chromosome.begin() + i, aux[l]);
                            }else{
                                int index = i - I + 1;
                                individual->chromosome.insert(individual->chromosome.begin() + index, aux[l]);
                            }
                        }

                        int q;
                        if (k < i) {
                            q = individual->updateCost(k, i);
                        } else {
                            q = individual->updateCost(i, k + I - 1);
                        }

                        preInsertion(parameters->F, parameters->P, k, I);  
                        individual->caminhoCritico();
                    }
                }
            }
      
        }
    
    }

    // if (individual->solutionCost.evaluation != individual->E[parameters->numJobs][parameters->numTools] || individual->solutionCost.evaluation != individual->Q[parameters->numJobs][parameters->numTools]) {
    //     cout << "Problema na vizinhança insertionGenerico!!";
    //     exit(1);
    // }

    // parameters->viz3.first = parameters->viz3.first + (cpuTime() - startTime);

    // for (int i = 0; i < parameters->numJobs+1; i++) {
    //     delete[] C[i];
    // }

    // delete[] C;

    return improved;

}

LocalSearch::LocalSearch() {
}

LocalSearch::LocalSearch(Parameters *parameters, Individual *individual) : parameters(parameters), individual(individual) {
    tempIndiv = new Individual(parameters);
    tempIndiv->recopyIndividual(tempIndiv, individual);
    runSearchTotal();
}

LocalSearch::~LocalSearch() {
    delete tempIndiv;
}

int genRandom(int i) {
    return std::rand() % i;
}

void LocalSearch::shuffleMoves() {
    random_shuffle(parameters->listaSwap.begin(), parameters->listaSwap.end(), genRandom);
}

void LocalSearch::shuffleIndices() {

    random_shuffle(parameters->moves.begin(), parameters->moves.end(), genRandom);
    
    // Shuffling the jobs order vector
    // random_shuffle(parameters->orderLSi.begin(), parameters->orderLSi.end(), genRandom);
    // random_shuffle(parameters->orderLSs.begin(), parameters->orderLSs.end(), genRandom);
    // random_shuffle(parameters->orderLSii.begin(), parameters->orderLSii.end(), genRandom);
}

unsigned int LocalSearch::calcIdleOrBlockedSwap(int p1, int p2) {

    // double startTime = cpuTime();

    int * aux = new int [parameters->numTools + 1];
    int total = individual->solutionCost.zeroBlocks;

    aux[0] = individual->E[p1][1];
    int i;
    for (int t = p1 + 1; t <= p1 + 1; t++) {
        i = p2 + 1;

        total = total - individual->idleOrBlocked[t-1];
        for (int j = 1; j < parameters->numTools + 1; j++) {
            if (j != parameters->numTools) {
                aux[j] = max(aux[j-1] + parameters->jobsToolsMatrix[individual->chromosome[i-1]][j-1], individual->E[t-1][j+1]);
            } else {
                aux[j] = aux[j-1] + parameters->jobsToolsMatrix[individual->chromosome[i-1]][j-1];
            }
            
            total = total + (aux[j] - individual->E[t-1][j] - parameters->jobsToolsMatrix[individual->chromosome[i-1]][j-1]);            
        }
    }

    for(int t = p1 + 2;t < individual->chromosome.size() + 1; t++) {

        if (t != p2 + 1) {
            i = t;
            total = total - individual->idleOrBlocked[t-1];
        } else {
            i = p1 + 1;
            total = total - individual->idleOrBlocked[t-1];
        }
        
        aux[0] = aux[1];
        for(int j = 1;j < parameters->numTools+1; j++) {
            total = total - aux[j];

            if (j != parameters->numTools) {
                aux[j] = max(aux[j-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[i-1]][j-1], aux[j + 1]);
            } else {
                aux[j] = aux[j - 1] + parameters->jobsToolsMatrix[individual->chromosome[i-1]][j-1];
            }

            total = total + aux[j] - parameters->jobsToolsMatrix[individual->chromosome[i-1]][j-1];
             
        }
    }

    delete[] aux;

    // parameters->tempoIdle = parameters->tempoIdle + (cpuTime() - startTime);

    return total;

}

unsigned int LocalSearch::calcIdleOrBlockedRelocate(int p1, int p2, int I) {

    // double startTime = cpuTime();

    int * aux = new int [parameters->numTools + 1];
    int total = individual->solutionCost.zeroBlocks;
    
    int i, p;

    if (p2 > p1) {
        i = p1 + 1 + I;
        p = p1;
    } else {
        i = p1 + 1;
        p = p2;
    }

    aux[0] = individual->E[p][1];

    for (int t = p + 1; t <= p + 1; t++) {
    
        total = total - individual->idleOrBlocked[t-1];
        for (int j = 1; j < parameters->numTools + 1; j++) {
            if (j != parameters->numTools) {
                aux[j] = max(aux[j-1] + parameters->jobsToolsMatrix[individual->chromosome[i-1]][j-1], individual->E[t-1][j+1]);
            } else {
                aux[j] = aux[j-1] + parameters->jobsToolsMatrix[individual->chromosome[i-1]][j-1];
            }

            total = total + (aux[j] - individual->E[t-1][j] - parameters->jobsToolsMatrix[individual->chromosome[i-1]][j-1]);            
        }
    }

    
    if (p2 > p1) {
        for(int t = p1 + 2;t < individual->chromosome.size() + 1; t++) {

            if (i < p2 + 1) {
                i++;
                //cout << i << " a ";
            } else if (i == p2 + 1) {
                //cout << i << " ";
                i = p1 + 1;
                //cout << i << " b ";
            } else {
                i = t;
                //cout << i << " c ";
            }

            total = total - individual->idleOrBlocked[t-1];
            
            aux[0] = aux[1];
            for(int j = 1;j < parameters->numTools+1; j++) {
                total = total - aux[j];

                if (j != parameters->numTools) {
                    aux[j] = max(aux[j-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[i-1]][j-1], aux[j + 1]);
                } else {
                    aux[j] = aux[j - 1] + parameters->jobsToolsMatrix[individual->chromosome[i-1]][j-1];
                }

                total = total + aux[j] - parameters->jobsToolsMatrix[individual->chromosome[i-1]][j-1];
                
            }

            if (i == p1 + I) {
                i = p2 + 1 + I;
            }
        
        }

        //cout << endl;
        //getchar();

        delete[] aux;
        return total;

    } else {
        for(int t = p2 + 2;t < individual->chromosome.size() + 1; t++) {

            if (I > 1 && t <= p2 + I) {
                i++;
            } else if (t < p1 + 1 + I) {
                i = t - I;
            } else {
                i = t;
            }

            total = total - individual->idleOrBlocked[t-1];

            aux[0] = aux[1];
            for(int j = 1;j < parameters->numTools+1; j++) {
                total = total - aux[j];

                if (j != parameters->numTools) {
                    aux[j] = max(aux[j-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[i-1]][j-1], aux[j + 1]);
                } else {
                    aux[j] = aux[j - 1] + parameters->jobsToolsMatrix[individual->chromosome[i-1]][j-1];
                }

                total = total + aux[j] - parameters->jobsToolsMatrix[individual->chromosome[i-1]][j-1];
                
            }

            if (i == p1 + I) {
                i = p2 + 1 + I;
            }
        
        }

        delete[] aux;
        return total;
    }

    // parameters->tempoIdle = parameters->tempoIdle + (cpuTime() - startTime);

}

void LocalSearch::preInsertion(int ** F, int ** P, unsigned int k, unsigned int I) {

    for (int i = 1; i < k + 1; i++) {
        for (int j = 1; j < parameters->numTools + 1; j++) {
            parameters->F[i][j] = individual->E[i][j];
        }
    }
        
    //Calcula os tempos de saída das máquinas na ordem direta retirando a k-ésima tarefa
    for (int i = k + 1; i < parameters->numJobs + 1 - I; i++) {
        parameters->F[i][0] = parameters->F[i-1][1];
        for (int j = 1; j < parameters->numTools + 1; j++) {
            if (j != parameters->numTools) {
                parameters->F[i][j] = max(parameters->F[i][j-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[i - 1 + I]][j-1], parameters->F[i-1][j+1]);
            } else {
                parameters->F[i][j] = parameters->F[i][j-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[i - 1 + I]][j-1];
            }
        }
    }

    int jobInicio = parameters->numJobs;
    for (int i = 1; jobInicio > k + I; i++) {
        for (int j = 1; j < parameters->numTools + 1; j++) {
            parameters->P[i][j] = individual->Q[i][j];
        }
        jobInicio--;
    }
        
    //Calcula os tempos de saída das máquinas na ordem inversa retirando a k-ésima tarefa
    for (int i = parameters->numJobs - jobInicio + 1; i < parameters->numJobs + 1 - I; i++) {
        parameters->P[i][0] = parameters->P[i-1][1];
        int maquinaInicio = parameters->numTools-1;
        for(int j=1;j<parameters->numTools+1;j++) {
            if (j != parameters->numTools) {
                parameters->P[i][j] = max(parameters->P[i][j-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[jobInicio - I - 1]][maquinaInicio], parameters->P[i-1][j+1]);
                maquinaInicio--;
            } else {
                parameters->P[i][j] = parameters->P[i][j-1] + (int)parameters->jobsToolsMatrix[individual->chromosome[jobInicio - I - 1]][maquinaInicio];
                maquinaInicio--;
            }
        }
        jobInicio--;        
    }


}

int LocalSearch::getInsertionLowerBound(const int & k, const int & i) {

    // int lowerBound = -1; // Valor do lower bound
    int delta = 0; // Delta do lower bound 
    int pi_zero; // pi_zero conforme definição 7 de Ding et al. (2016)*
    int bInsercao; // Bloco onde a inserção será realizada

    // Determina o bloco onde a inserção será realizada
    if (individual->jobsInfo[i].nBlocos == 1 || (i == 0 || i == parameters->numJobs - 1)) {
        bInsercao = individual->jobsInfo[i].i1;
    } else if (individual->jobsInfo[i].nBlocos == 2) {
        if (i < k) {
            bInsercao = individual->jobsInfo[i].i1 + 1;
        } else {
            bInsercao = individual->jobsInfo[i].i1;
        }
    } else {
        if (i < k) {
            bInsercao = individual->jobsInfo[i].i1 + 1;
        } else {
            bInsercao = individual->jobsInfo[i].i1 - 1;
        }
    }

     // Definição 7 de Ding et al. (2016)
    /*if ((bInsercao < individual->jobsInfo[k].b && i != individual->blocos[bInsercao].jf) || (bInsercao > individual->jobsInfo[k].b  && i != individual->blocos[bInsercao].jf - 1)) {
        pi_zero = individual->chromosome[individual->blocos[bInsercao].jf - 1];
    } else {
        pi_zero = individual->chromosome[k];
    }*/

    pi_zero = individual->chromosome[individual->blocos[bInsercao].jf];

    // Teorema 4 de Ding et al. (2016)
    if (individual->blocos[individual->jobsInfo[k].b].tipo == 0) {

        for (int i = individual->blocos[individual->jobsInfo[k].b].sm; i <= individual->blocos[individual->jobsInfo[k].b].em; i++) {
            delta = delta - (int)parameters->jobsToolsMatrix[individual->chromosome[k]][i];
        }

        if (k != 0) {
            for (int i = individual->blocos[individual->jobsInfo[k].b].sm + 1; i <= individual->blocos[individual->jobsInfo[k].b].em; i++) {
                delta = delta + (int)parameters->jobsToolsMatrix[individual->chromosome[k-1]][i];
            }
        } else {
            for (int i = individual->blocos[individual->jobsInfo[k].b].em - 1; i >= individual->blocos[individual->jobsInfo[k].b].sm; i--) {
                delta = delta + (int)parameters->jobsToolsMatrix[individual->chromosome[k+1]][i];
            }
        }
        
        if (individual->blocos[bInsercao].tipo == 1 || (i == 0 || i == parameters->numJobs-1)) {
            
            if (i != parameters->numJobs-1) {
                delta = delta + (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[bInsercao].sm];
            } else {
                delta = delta + (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[bInsercao].em];
            }

            // lowerBound = individual->solutionCost.evaluation + delta;
            return delta;
        } else if (individual->blocos[bInsercao].tipo == -1) {
            delta = delta + (int)parameters->jobsToolsMatrix[pi_zero][individual->blocos[bInsercao].em];
           
            // lowerBound = individual->solutionCost.evaluation + delta;
            return delta;
        }

    } 

    // Teorema 2 de Ding et al. (2016)
    if (individual->blocos[individual->jobsInfo[k].b].tipo == 1) {

        if (individual->blocos[bInsercao].tipo == 1 || (i == 0 || i == parameters->numJobs-1)) {
            
            if (i != parameters->numJobs - 1) {
                delta = (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[bInsercao].sm] - (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[individual->jobsInfo[k].b].sm];
            } else {
                delta = (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[bInsercao].em] - (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[individual->jobsInfo[k].b].sm];
            }

            // lowerBound = individual->solutionCost.evaluation + delta;
            return delta;

        } else if (individual->blocos[bInsercao].tipo == -1) {

            delta = (int)parameters->jobsToolsMatrix[pi_zero][individual->blocos[bInsercao].em] - (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[individual->jobsInfo[k].b].sm];

            // lowerBound = individual->solutionCost.evaluation + delta;
            return delta;

           /* 
            if (k != individual->blocos[individual->jobsInfo[i].b].ji) {
                delta = (int)parameters->jobsToolsMatrix[pi_zero][individual->blocos[bInsercao].em] - (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[individual->jobsInfo[k].b].sm];

                lowerBound = individual->solutionCost.evaluation + delta;
                return lowerBound;
            } else {
                delta = (int)parameters->jobsToolsMatrix[individual->chromosome[k+1]][individual->blocos[individual->jobsInfo[k].b].sm] - (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[individual->jobsInfo[k].b].sm];

                lowerBound = individual->solutionCost.evaluation + delta;
                return lowerBound;
            }*/

            /*if (!(individual->blocos[individual->jobsInfo[k].b].joined == true && k == individual->blocos[individual->jobsInfo[k].b].ji)) {
                
                delta = (int)parameters->jobsToolsMatrix[pi_zero][individual->blocos[bInsercao].em] - (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[individual->jobsInfo[k].b].sm];
            } else {
                cout << "a";
                if (individual->jobsInfo[k].b == bInsercao - 1) {
                    return individual->solutionCost.evaluation;
                }
                
                delta = (int)parameters->jobsToolsMatrix[pi_zero][individual->blocos[bInsercao].em] - (int)parameters->jobsToolsMatrix[individual->chromosome[k+1]][individual->blocos[individual->jobsInfo[k].b].sm];
            }*/
        }
    }

    // Teorema 3 de Ding et al. (2016)
    if (individual->blocos[individual->jobsInfo[k].b].tipo == -1) {

        if (individual->blocos[bInsercao].tipo == -1) {
            pi_zero = individual->chromosome[individual->blocos[bInsercao].jf];
            
            if (individual->blocos[individual->jobsInfo[k].b].joined == true) {
                //cout << "a" << endl;
                delta = (int)parameters->jobsToolsMatrix[pi_zero][individual->blocos[bInsercao].em] - (int)parameters->jobsToolsMatrix[individual->chromosome[individual->blocos[individual->jobsInfo[k].b].ji]][individual->blocos[individual->jobsInfo[k].b].sm];

                // lowerBound = individual->solutionCost.evaluation + delta;
                return delta;
            } else {
                //cout << "b" << endl;
                delta = (int)parameters->jobsToolsMatrix[pi_zero][individual->blocos[bInsercao].em] - (int)parameters->jobsToolsMatrix[individual->chromosome[individual->blocos[individual->jobsInfo[k].b+1].ji]][individual->blocos[individual->jobsInfo[k].b+1].em];

                // lowerBound = individual->solutionCost.evaluation + delta;
                return delta;
            }

        } else if (individual->blocos[bInsercao].tipo == 1 || (i == 0 || i == parameters->numJobs-1)) {
            
            if (individual->blocos[individual->jobsInfo[k].b].joined == true) {
                if (i != parameters->numJobs - 1) {
                    delta = (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[bInsercao].sm] - (int)parameters->jobsToolsMatrix[individual->chromosome[individual->blocos[individual->jobsInfo[k].b].ji]][individual->blocos[individual->jobsInfo[k].b].sm];
                } else {
                    delta = (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[bInsercao].em] - (int)parameters->jobsToolsMatrix[individual->chromosome[individual->blocos[individual->jobsInfo[k].b].ji]][individual->blocos[individual->jobsInfo[k].b].sm];
                }

                // lowerBound = individual->solutionCost.evaluation + delta;
                return delta;
            } else {
                if (i != parameters->numJobs - 1) {
                    delta = (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[bInsercao].sm] - (int)parameters->jobsToolsMatrix[individual->chromosome[individual->blocos[individual->jobsInfo[k].b+1].ji]][individual->blocos[individual->jobsInfo[k].b+1].em];
                } else {
                    delta = (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[bInsercao].em] - (int)parameters->jobsToolsMatrix[individual->chromosome[individual->blocos[individual->jobsInfo[k].b+1].ji]][individual->blocos[individual->jobsInfo[k].b+1].em];
                }

                // lowerBound = individual->solutionCost.evaluation + delta;
                return delta;
            }

            /*if ((bInsercao != individual->jobsInfo[k].b - 1) || (bInsercao == individual->jobsInfo[k].b - 1 && i != individual->blocos[bInsercao].ji)) {
                delta = (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[bInsercao].sm] - (int)parameters->jobsToolsMatrix[individual->blocos[individual->jobsInfo[k].b - 1].jf - 1][individual->blocos[individual->jobsInfo[k].b].sm];
            }*/
        }
    }

    // return lowerBound;

}

int LocalSearch::getSwapLowerBound(const int & k, const int & i) {
    
    // int lowerBound = -1; // Valor do lower bound
    int delta = 0; // Delta do lower bound 

    // Calcula o swap lower bound para um job em uma T-sequence
    if (individual->blocos[individual->jobsInfo[k].b].tipo == 0) {

        // Retira os tempos de processamento da T-sequence original e adiciona o da nova T-sequence
        for (int l = individual->blocos[individual->jobsInfo[k].b].sm; l <= individual->blocos[individual->jobsInfo[k].b].em; l++) {
            delta = delta + (int)parameters->jobsToolsMatrix[individual->chromosome[i]][l] - (int)parameters->jobsToolsMatrix[individual->chromosome[k]][l];
        }

        if (individual->blocos[individual->jobsInfo[i].b].tipo == 1) {
            delta = delta + (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[individual->jobsInfo[i].b].sm] - (int)parameters->jobsToolsMatrix[individual->chromosome[i]][individual->blocos[individual->jobsInfo[i].b].sm];  
        } else if (individual->blocos[individual->jobsInfo[i].b].tipo == 0){
            for (int l = individual->blocos[individual->jobsInfo[i].b].sm; l <= individual->blocos[individual->jobsInfo[i].b].em; l++) {
                delta = delta - (int)parameters->jobsToolsMatrix[individual->chromosome[i]][l];
            }
        }

        // lowerBound = individual->solutionCost.evaluation + delta;
        return delta;

    } 

    // Calcula o swap lower bound para o job em um normal block
    if (individual->blocos[individual->jobsInfo[k].b].tipo == 1) {
        
        // Calcula o delta para um normal block
        if (individual->blocos[individual->jobsInfo[i].b].tipo == 1) {
            
            delta = delta + (int)parameters->jobsToolsMatrix[individual->chromosome[i]][individual->blocos[individual->jobsInfo[k].b].sm] + (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[individual->jobsInfo[i].b].sm] - (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[individual->jobsInfo[k].b].sm] - (int)parameters->jobsToolsMatrix[individual->chromosome[i]][individual->blocos[individual->jobsInfo[i].b].sm];
            
        } else if (individual->blocos[individual->jobsInfo[i].b].tipo == -1) {  // Calcula o delta para um anti-block
            
            delta = delta + (int)parameters->jobsToolsMatrix[individual->chromosome[i]][individual->blocos[individual->jobsInfo[k].b].sm] - (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[individual->jobsInfo[k].b].sm];

        } else { 
            
            for (int l = individual->blocos[individual->jobsInfo[i].b].sm; l <= individual->blocos[individual->jobsInfo[i].b].em; l++) {
                delta = delta - (int)parameters->jobsToolsMatrix[individual->chromosome[i]][l] + (int)parameters->jobsToolsMatrix[individual->chromosome[k]][l];
            }

            delta = delta + (int)parameters->jobsToolsMatrix[individual->chromosome[i]][individual->blocos[individual->jobsInfo[k].b].sm] - (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[individual->jobsInfo[k].b].sm];
            
        }

        // lowerBound = individual->solutionCost.evaluation + delta;
        return delta;
    }

    // Calcula o swap lower bound para o job em um anti-block
    if (individual->blocos[individual->jobsInfo[k].b].tipo == -1) {
        
        if (individual->blocos[individual->jobsInfo[i].b].tipo == 1) {
            
            delta = delta + (int)parameters->jobsToolsMatrix[individual->chromosome[k]][individual->blocos[individual->jobsInfo[i].b].sm] - (int)parameters->jobsToolsMatrix[individual->chromosome[i]][individual->blocos[individual->jobsInfo[i].b].sm];
            
        } else if (individual->blocos[individual->jobsInfo[i].b].tipo == 0) {
            
            for (int l = individual->blocos[individual->jobsInfo[i].b].sm; l <= individual->blocos[individual->jobsInfo[i].b].em; l++) {
                delta = delta - (int)parameters->jobsToolsMatrix[individual->chromosome[i]][l] + (int)parameters->jobsToolsMatrix[individual->chromosome[k]][l];
            }

        }
        
        // lowerBound = individual->solutionCost.evaluation + delta;
        return delta;
    }

    // return lowerBound;

}