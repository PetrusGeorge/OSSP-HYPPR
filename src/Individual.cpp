/**
 * @file Individual.cpp
 *
 * Implements Individual class methods
 *
 * @author Jordana Mecler
 *
 * Contact: jmecler@inf.puc-rio.br
 *
 */

#include "Individual.h"

Individual::Individual(Parameters *parameters, int elementoInicial) : parameters(parameters) {

    age = 0;
    isFitnessComputed = false;

    vector<vector<unsigned int>> E(parameters->numJobs + 1, vector<unsigned int>(parameters->numTools + 1));
    vector<vector<unsigned int>> Q(parameters->numJobs + 1, vector<unsigned int>(parameters->numTools + 1));
    
    this->E = E;
    this->Q = Q;

    vector<unsigned int> IB(parameters->numJobs, 0);
    idleOrBlocked = IB;
    solutionCost.zeroBlocks = 0; 
    
    //randomSolution();

    PF_NEH(elementoInicial);

    // int escolha = rand() % 2;
    // if (escolha == 0) {
    //     PF_NEH(elementoInicial);
    // } else {
    //     NEH(elementoInicial);
    // }
} 

void Individual::NEH(int elementoInicial) {

    vector <int> tProcessamento;
    vector <int> auxiliar;

    for (int i = 0; i < parameters->numJobs; i++) {
        tProcessamento.push_back(0);
        for (int j = 0; j < parameters->numTools; j++) {
            tProcessamento[i] += parameters->jobsToolsMatrix[i][j];
        }
    }

    auxiliar = tProcessamento;
    sort(auxiliar.begin(), auxiliar.end());

    int elementosAlocados = 0;
    int elemento = find(tProcessamento.begin(), tProcessamento.end(), auxiliar[elementoInicial]) - tProcessamento.begin();
    chromosome.push_back(elemento);
    edgesIndividuals.push_back(0);

    tProcessamento[elemento] = -1;
    elementosAlocados++;

    for (int i = 1; i < chromosome.size() + 1; i++) {
        this->E[i][0] = this->E[i-1][1];
        idleOrBlocked[i-1] = 0;
        for (int j = 1; j < parameters->numTools + 1; j++) {
            if (j != parameters->numTools) {
                this->E[i][j] = max(this->E[i][j-1] + parameters->jobsToolsMatrix[chromosome[i-1]][j-1], this->E[i-1][j+1]);
            } else {
                this->E[i][j] = this->E[i][j-1] + parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
            }      

            idleOrBlocked[i-1] += this->E[i][j] - this->E[i-1][j] - parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
        }
        solutionCost.zeroBlocks += idleOrBlocked[i - 1];
    }

    int jobInicio = chromosome.size() - 1;
    int maquinaInicio;
    for(int i=1;i<chromosome.size()+1;i++) {
        this->Q[i][0] = this->Q[i-1][1];
        maquinaInicio = parameters->numTools - 1;
        for(int j=1;j<parameters->numTools+1;j++) {
            
            if (j != parameters->numTools) {
                this->Q[i][j] = max(this->Q[i][j-1] + parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio], this->Q[i-1][j+1]);
                maquinaInicio--;
            } else {
                this->Q[i][j] = this->Q[i][j-1] + parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio];
                maquinaInicio--;
            }
        }
        jobInicio--;  
    }

    int ** C = new int * [parameters->numJobs + 1];
    for (int i = 0; i < parameters->numJobs + 1; i++) {
        C[i] = new int[parameters->numTools + 1];
    }

    //limpaMatriz(C, parameters->numJobs + 1, parameters->numTools + 1);   

    int maximo = -9999;
    int p, vAux, guardaValor;
    while (elementosAlocados < parameters->numJobs) {
        int mkspAux = 99999;
        elemento = max_element(tProcessamento.begin(), tProcessamento.end()) - tProcessamento.begin();

        tProcessamento[elemento] = -1;
        elementosAlocados++;

        for (int i = 0; i < chromosome.size() + 1; i++) {
                vAux = parameters->numTools;
                for (int j = 1; j < parameters->numTools + 1; j++) {
                    if (j == 1) {
                        guardaValor = max(this->Q[chromosome.size() - i][j] + parameters->jobsToolsMatrix[elemento][parameters->numTools - j], this->Q[chromosome.size()- i][j + 1]);
                        C[i][j] = this->E[i][vAux] + guardaValor;
                    } else {
                        guardaValor = max(guardaValor + parameters->jobsToolsMatrix[elemento][parameters->numTools - j], this->Q[chromosome.size() - i][j + 1]);
                        C[i][j] = this->E[i][vAux] + guardaValor;
                    }

                    if (C[i][j] > maximo) {
                        maximo = C[i][j];
                    }

                    vAux--;
                }
                
                if (maximo < mkspAux) {
                    mkspAux = maximo;
                    p = i;
                }
                maximo = -9999;
        }
        
        if (elementosAlocados == parameters->numJobs) {
            solutionCost.evaluation = mkspAux;
        }

        chromosome.insert(chromosome.begin() + p, elemento);
        edgesIndividuals.insert(edgesIndividuals.begin() + p, 0);

        for (int i = p+1; i < chromosome.size() + 1; i++) {
            this->E[i][0] = this->E[i-1][1];
            solutionCost.zeroBlocks -= idleOrBlocked[i-1];
            idleOrBlocked[i-1] = 0;
            for (int j = 1; j < parameters->numTools + 1; j++) {
                if (j != parameters->numTools) {
                    this->E[i][j] = max(this->E[i][j-1] + parameters->jobsToolsMatrix[chromosome[i-1]][j-1], this->E[i-1][j+1]);
                } else {
                    this->E[i][j] = this->E[i][j-1] + parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
                }

                idleOrBlocked[i-1] += this->E[i][j] - this->E[i-1][j] - parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
            }
        
            solutionCost.zeroBlocks += idleOrBlocked[i - 1];

        }

        jobInicio = p;
        maquinaInicio;
        for(int i = chromosome.size() - p;i < chromosome.size() + 1; i++) {
            this->Q[i][0] = this->Q[i-1][1];
            maquinaInicio = parameters->numTools-1;
            for(int j = 1;j < parameters->numTools+1; j++) {
                if (j != parameters->numTools) {
                    this->Q[i][j] = max(this->Q[i][j-1] + parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio], this->Q[i-1][j+1]);
                    maquinaInicio--;
                } else {
                    this->Q[i][j] = this->Q[i][j-1] + parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio];
                    maquinaInicio--;
                }
            }
            jobInicio--;  
        }

    }
    
    for (int i = 0; i < parameters->numJobs+1; i++) {
        delete[] C[i];
    }

    delete[] C;

    jobsInfo.resize(parameters->numJobs);
    caminhoCritico();

    if (this->E[parameters->numJobs][parameters->numTools] != solutionCost.evaluation || this->Q[parameters->numJobs][parameters->numTools] != solutionCost.evaluation) {
        cout << this->E[parameters->numJobs][parameters->numTools] << " " << this->Q[parameters->numJobs][parameters->numTools] << " " << solutionCost.evaluation << endl;
        cout << "PROBLEMA NA NEH!!!" << endl << endl;
        exit(1);
    }


}

void Individual::random_swap() {

    int jobsTrocados = 4;
    int contador = 1;
    //const int minJobs = 8;
    //const int maxJobs = parameters->numJobs / 4;
    //int jobsMovidos = minJobs + rand () % (maxJobs - 1);

    int menorJob = 9999;
    int maiorJob = -9999;

    while (contador <= jobsTrocados) {
        int pos1 = rand() % chromosome.size();
        int pos2 = rand() % chromosome.size();

        while (pos2 == pos1) {
            pos2 = rand() % chromosome.size();
        }

        if (pos1 < pos2) {
            if (pos1 < menorJob) {
                menorJob = pos1;
            }

            if (pos2 > maiorJob) {
                maiorJob = pos2;
            }
        } else {
            if (pos2 < menorJob) {
                menorJob = pos2;
            }

            if (pos1 > maiorJob) {
                maiorJob = pos1;
            }            
        }

        swap(chromosome[pos1], chromosome[pos2]);
        contador++;
    }

    solutionCost.evaluation = updateCost(menorJob, maiorJob);
    caminhoCritico();

}

void Individual::profileFitting(int elemento) {

    unsigned int ** C = new unsigned int * [parameters->numJobs + 1];
    for (int i = 0; i < parameters->numJobs + 1; i++) {
         C[i] = new unsigned int [parameters->numTools + 1];
    }

    for(int i = 0; i < parameters->numJobs + 1; i++) {
        this->E[i][0] = 0;
        this->Q[i][0] = 0;
        C[i][0] = 0;
    }

    for(int j = 0; j < parameters->numTools + 1; j++) {
        this->E[0][j] = 0;
        this->Q[0][j] = 0;
        C[0][j] = 0;
    }
    
    chromosome.push_back(elemento);
    edgesIndividuals.push_back(0);
    
    vector <int> LC;
    for (int i = 0; i < parameters->numJobs; i++) {
        if (i != elemento) {
            LC.push_back(i);
        }
    }

    int elementosAlocados = 1;

    for(int i = 1;i <= elementosAlocados; i++) {
        idleOrBlocked[i-1] = 0;
        this->E[i][0] = this->E[i-1][1];
        for(int j = 1;j < parameters->numTools+1; j++) {
            if (j != parameters->numTools) {
                this->E[i][j] = max(this->E[i][j-1] + (int)parameters->jobsToolsMatrix[elemento][j-1], this->E[i-1][j+1]);
            } else {
                this->E[i][j] = this->E[i][j-1] + (int)parameters->jobsToolsMatrix[elemento][j-1];
            }
            idleOrBlocked[i-1] = idleOrBlocked[i-1] + this->E[i][j] - this->E[i-1][j] - (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
        }
        solutionCost.zeroBlocks = solutionCost.zeroBlocks + idleOrBlocked[i-1];
    }
    
    int idle_blocked, decisao;
    while (elementosAlocados < parameters->numJobs) {
        idle_blocked = 99999;
        for (int i = 1; i < LC.size() + 1; i++) {
            decisao = 0;
            for (int j = 1; j < parameters->numTools + 1; j++) {
                if (j != parameters->numTools) {
                    C[i][j] = max(max(this->E[elementosAlocados][j], C[i][j-1]) + (int)parameters->jobsToolsMatrix[LC[i-1]][j-1], this->E[elementosAlocados][j+1]);
                } else {
                    C[i][j] = C[i][j-1] + (int)parameters->jobsToolsMatrix[LC[i-1]][j-1];
                }

                decisao += C[i][j] - this->E[elementosAlocados][j] - (int)parameters->jobsToolsMatrix[LC[i-1]][j-1];
                
            }

            if (decisao < idle_blocked) {
                idle_blocked = decisao;
                elemento = i - 1;
            }
            
        }
        
        chromosome.push_back(LC[elemento]);
        LC.erase(LC.begin() + elemento);

        edgesIndividuals.push_back(0);
        elementosAlocados++;

        for (int j = 1; j < parameters->numTools + 1; j++) {
            this->E[elementosAlocados][j] = C[elemento + 1][j];

            idleOrBlocked[elementosAlocados - 1] = idleOrBlocked[elementosAlocados - 1] + this->E[elementosAlocados][j] - this->E[elementosAlocados-1][j] - (int)parameters->jobsToolsMatrix[chromosome[elementosAlocados - 1]][j-1];
        }

        solutionCost.zeroBlocks = solutionCost.zeroBlocks + idleOrBlocked[elementosAlocados - 1];

        if (elementosAlocados == parameters->numJobs) {
            solutionCost.evaluation = C[elemento + 1][parameters->numTools];
        }

    }

    int jobInicio = parameters->numJobs - 1;
    int maquinaInicio;
    for(int i=1;i<parameters->numJobs+1;i++) {
        this->Q[i][0] = this->Q[i-1][1];
        maquinaInicio = parameters->numTools - 1;
        for(int j=1;j<parameters->numTools+1;j++) {
            
            if (j != parameters->numTools) {
                this->Q[i][j] = max(this->Q[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio], this->Q[i-1][j+1]);
                maquinaInicio--;
            } else {
                this->Q[i][j] = this->Q[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio];
                maquinaInicio--;
            }
        }
        jobInicio--;  
    }

    for (int i = 0; i < parameters->numJobs+1; i++) {
        delete[] C[i];
    }

    delete[] C;

    // jobsInfo.resize(parameters->numJobs);
    // caminhoCritico();

    if (this->E[parameters->numJobs][parameters->numTools] != solutionCost.evaluation || this->Q[parameters->numJobs][parameters->numTools] != solutionCost.evaluation) {
        cout << "PROBLEMA NA PROFILE FITTING!!!" << endl << endl;
        exit(1);
    }

}

void Individual::randomSolution() {

    vector <int> lJobs;
    for (int i = 0; i < parameters->numJobs; i++) {
        lJobs.push_back(i);
        edgesIndividuals.push_back(0);
    }

    int seleciona;

    while(!lJobs.empty()) {
        seleciona = rand() % lJobs.size();

        chromosome.push_back(lJobs[seleciona]);
        lJobs.erase(lJobs.begin() + seleciona);
    }

    for(int i = 0; i < parameters->numJobs + 1; i++) {
        this->E[i][0] = 0;
        this->Q[i][0] = 0;
    }

    for(int j = 0; j < parameters->numTools + 1; j++) {
        this->E[0][j] = 0;
        this->Q[0][j] = 0;
    }

    solutionCost.zeroBlocks = 0;

    for (int i = 1; i < parameters->numJobs + 1; i++) {
        idleOrBlocked[i-1] = 0;
        this->E[i][0] = this->E[i-1][1];
        for (int j = 1; j < parameters->numTools + 1; j++) {
            if (j != parameters->numTools) {
                this->E[i][j] = max(this->E[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1], this->E[i-1][j+1]);
            } else {
                this->E[i][j] = this->E[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
            }
            
            idleOrBlocked[i-1] = idleOrBlocked[i-1] + this->E[i][j] - this->E[i-1][j] - (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
        }
        solutionCost.zeroBlocks = solutionCost.zeroBlocks + idleOrBlocked[i-1];
    }

    int jobInicio = parameters->numJobs - 1;
    int maquinaInicio;
    for(int i=1;i<parameters->numJobs+1;i++) {
        this->Q[i][0] = this->Q[i-1][1];
        maquinaInicio = parameters->numTools - 1;
        for(int j=1;j<parameters->numTools+1;j++) {
            
            if (j != parameters->numTools) {
                this->Q[i][j] = max(this->Q[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio], this->Q[i-1][j+1]);
                maquinaInicio--;
            } else {
                this->Q[i][j] = this->Q[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio];
                maquinaInicio--;
            }
        }
        jobInicio--;  
    }

    solutionCost.evaluation = this->E[parameters->numJobs][parameters->numTools];

    if (this->E[parameters->numJobs][parameters->numTools] != this->Q[parameters->numJobs][parameters->numTools]) {
        cout << "OPA!!!" << endl;
        exit(1);
    }

        jobsInfo.resize(parameters->numJobs);
    caminhoCritico();


}

void Individual::profileFittingT(int elemento) {

    unsigned int ** C = new unsigned int * [parameters->numJobs + 1];
    for (int i = 0; i < parameters->numJobs + 1; i++) {
         C[i] = new unsigned int [parameters->numTools + 1];
    }

    for(int i = 0; i < parameters->numJobs + 1; i++) {
        this->E[i][0] = 0;
        this->Q[i][0] = 0;
        C[i][0] = 0;
    }

    for(int j = 0; j < parameters->numTools + 1; j++) {
        this->E[0][j] = 0;
        this->Q[0][j] = 0;
        C[0][j] = 0;
    }
    
    chromosome.push_back(elemento);
    edgesIndividuals.push_back(0);
    
    vector <int> LC;
    for (int i = 0; i < parameters->numJobs; i++) {
        if (i != elemento) {
            LC.push_back(i);
        }
    }

    int elementosAlocados = 1;

    for(int i = 1;i <= elementosAlocados; i++) {
        idleOrBlocked[i-1] = 0;
        this->E[i][0] = this->E[i-1][1];
        for(int j = 1;j < parameters->numTools+1; j++) {
            if (j != parameters->numTools) {
                this->E[i][j] = max(this->E[i][j-1] + (int)parameters->jobsToolsMatrix[elemento][j-1], this->E[i-1][j+1]);
            } else {
                this->E[i][j] = this->E[i][j-1] + (int)parameters->jobsToolsMatrix[elemento][j-1];
            }
            idleOrBlocked[i-1] = idleOrBlocked[i-1] + this->E[i][j] - this->E[i-1][j] - (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
        }
        solutionCost.zeroBlocks = solutionCost.zeroBlocks + idleOrBlocked[i-1];
    }
    
    int idle_blocked, decisao;
    while (elementosAlocados < parameters->numJobs) {
        idle_blocked = 99999;
        for (int i = 1; i < LC.size() + 1; i++) {
            decisao = 0;
            for (int j = 1; j < parameters->numTools + 1; j++) {
                if (j != parameters->numTools) {
                    C[i][j] = max(max(this->E[elementosAlocados][j], C[i][j-1]) + (int)parameters->jobsToolsMatrix[LC[i-1]][j-1], this->E[elementosAlocados][j+1]);
                } else {
                    C[i][j] = C[i][j-1] + (int)parameters->jobsToolsMatrix[LC[i-1]][j-1];
                }

                decisao += C[i][j] - this->E[elementosAlocados][j] - (int)parameters->jobsToolsMatrix[LC[i-1]][j-1];
                
            }

            decisao = ((parameters->numJobs - elementosAlocados + 1 - 2) * decisao) + C[i][parameters->numTools]; 

            if (decisao < idle_blocked) {
                idle_blocked = decisao;
                elemento = i - 1;
            }
            
        }
        
        chromosome.push_back(LC[elemento]);
        LC.erase(LC.begin() + elemento);

        edgesIndividuals.push_back(0);
        elementosAlocados++;

        for (int j = 1; j < parameters->numTools + 1; j++) {
            this->E[elementosAlocados][j] = C[elemento + 1][j];

            idleOrBlocked[elementosAlocados - 1] = idleOrBlocked[elementosAlocados - 1] + this->E[elementosAlocados][j] - this->E[elementosAlocados-1][j] - (int)parameters->jobsToolsMatrix[chromosome[elementosAlocados - 1]][j-1];
        }

        solutionCost.zeroBlocks = solutionCost.zeroBlocks + idleOrBlocked[elementosAlocados - 1];

        if (elementosAlocados == parameters->numJobs) {
            solutionCost.evaluation = C[elemento + 1][parameters->numTools];
        }

    }

    int jobInicio = parameters->numJobs - 1;
    int maquinaInicio;
    for(int i=1;i<parameters->numJobs+1;i++) {
        this->Q[i][0] = this->Q[i-1][1];
        maquinaInicio = parameters->numTools - 1;
        for(int j=1;j<parameters->numTools+1;j++) {
            
            if (j != parameters->numTools) {
                this->Q[i][j] = max(this->Q[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio], this->Q[i-1][j+1]);
                maquinaInicio--;
            } else {
                this->Q[i][j] = this->Q[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio];
                maquinaInicio--;
            }
        }
        jobInicio--;  
    }

    for (int i = 0; i < parameters->numJobs+1; i++) {
        delete[] C[i];
    }

    delete[] C;

    // jobsInfo.resize(parameters->numJobs);
    // caminhoCritico();

    if (this->E[parameters->numJobs][parameters->numTools] != solutionCost.evaluation || this->Q[parameters->numJobs][parameters->numTools] != solutionCost.evaluation) {
        cout << "PROBLEMA NA PROFILE FITTING!!!" << endl << endl;
        exit(1);
    }

}

void Individual::PFT_NEH(int e) {

    profileFittingT(e);
    
    int gamma;
    if (parameters->numJobs > 20) {
        gamma = 25;
    } else {
        gamma = 21;
    }

    int menorJobMovido = parameters->numJobs - gamma + 1;
    int maiorJobMovido = parameters->numJobs - 1;

    vector <int> LJ;

    for (int i = parameters->numJobs - 1; i >= parameters->numJobs - gamma; i--) {
        LJ.push_back(chromosome[i]);

        solutionCost.zeroBlocks = solutionCost.zeroBlocks - idleOrBlocked[i];
        chromosome.erase(chromosome.begin() + i);
        idleOrBlocked[i] = 0;
    }

    int p;
    bool improved = false;

    for(int i = menorJobMovido + 1;i < chromosome.size() + 1; i++) {
        solutionCost.zeroBlocks = solutionCost.zeroBlocks - idleOrBlocked[i-1];
        idleOrBlocked[i-1] = 0;
        
        this->E[i][0] = this->E[i-1][1];
        for(int j = 1;j < parameters->numTools+1; j++) {
            if (j != parameters->numTools) {
                this->E[i][j] = max(this->E[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1], this->E[i-1][j+1]);
            } else {
                this->E[i][j] = this->E[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
            }

            idleOrBlocked[i-1] = idleOrBlocked[i-1] + this->E[i][j] - this->E[i-1][j] - (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
        }
        
        solutionCost.zeroBlocks = solutionCost.zeroBlocks + idleOrBlocked[i-1];
        // sol.jobsInfo[i-1].posicao = i - 1;
    }

    int jobInicio = chromosome.size() - 1 - (parameters->numJobs - maiorJobMovido - 1);
    int maquinaInicio, aux;
    
    aux = jobInicio;
    for(int i = chromosome.size() - aux;i < chromosome.size() + 1; i++) {
        this->Q[i][0] = this->Q[i-1][1];
        maquinaInicio = parameters->numTools-1;
        for(int j = 1;j < parameters->numTools+1; j++) {
            if (j != parameters->numTools) {
                this->Q[i][j] = max(this->Q[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio], this->Q[i-1][j+1]);
                maquinaInicio--;
            } else {
                this->Q[i][j] = this->Q[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio];
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
        
        for (int i = 0; i < chromosome.size() + 1; i++) {
                vAux = parameters->numTools;

                for (int j = 1; j <= 1; j++) {
                    guardaValor = max(this->Q[chromosome.size() - i][j] + (int)parameters->jobsToolsMatrix[LJ[k]][parameters->numTools - j], this->Q[chromosome.size() - i][j + 1]);
                    C[i][j] = this->E[i][vAux] + guardaValor;
                
                    if (C[i][j] > maximo) {
                        maximo = C[i][j];
                    }

                    vAux--;
                
                }

                for (int j = 2; j < parameters->numTools; j++) {
                    guardaValor = max(guardaValor + parameters->jobsToolsMatrix[LJ[k]][parameters->numTools - j], this->Q[chromosome.size() - i][j + 1]);
                    C[i][j] = this->E[i][vAux] + guardaValor;
                    
                    if (C[i][j] > maximo) {
                        maximo = C[i][j];
                    }

                    vAux--;
                
                }

                for (int j = parameters->numTools; j <= parameters->numTools; j++) {
                    C[i][j] = this->E[i][vAux] + guardaValor + (int)parameters->jobsToolsMatrix[LJ[k]][parameters->numTools - j];
                
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
    
        chromosome.insert(chromosome.begin() + p, LJ[k]);
        LJ.erase(LJ.begin() + k);

        for(int i = p + 1; i < chromosome.size() + 1; i++) {
            this->E[i][0] = this->E[i-1][1];
            solutionCost.zeroBlocks = solutionCost.zeroBlocks - idleOrBlocked[i-1];
            idleOrBlocked[i-1] = 0;
            for(int j = 1;j < parameters->numTools+1; j++) {
                if (j != parameters->numTools) {
                    this->E[i][j] = max(this->E[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1], this->E[i-1][j+1]);
                } else {
                    this->E[i][j] = this->E[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
                }

                idleOrBlocked[i-1] = idleOrBlocked[i-1] + this->E[i][j] - this->E[i-1][j] - (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
            }
            // sol.jobsInfo[i-1].posicao = i-1;
            solutionCost.zeroBlocks = solutionCost.zeroBlocks + idleOrBlocked[i-1];
            // sol.jobsInfo[i-1].tempo_sistema = this->E[i][parameters->numTools] - this->E[i][1];
        }
        
        jobInicio = p;
        for(int i = chromosome.size() - p;i < chromosome.size() + 1; i++) {
            this->Q[i][0] = this->Q[i-1][1];
            maquinaInicio = parameters->numTools-1;
            for(int j = 1;j < parameters->numTools+1; j++) {
                if (j != parameters->numTools) {
                    this->Q[i][j] = max(this->Q[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio], this->Q[i-1][j+1]);
                    maquinaInicio--;
                } else {
                    this->Q[i][j] = this->Q[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio];
                    maquinaInicio--;
                }
            }
            jobInicio--;  
        }

    }

     solutionCost.evaluation = this->E[parameters->numJobs][parameters->numTools];

    for (int i = 0; i < parameters->numJobs+1; i++) {
        delete[] C[i];
    }

    delete[] C;

    jobsInfo.resize(parameters->numJobs);
    caminhoCritico();

    if (this->E[parameters->numJobs][parameters->numTools] !=  solutionCost.evaluation || this->Q[parameters->numJobs][parameters->numTools] !=  solutionCost.evaluation) {
        cout << "PROBLEMA NA PF_NEH!!!" << endl;
        exit(1);
    }



}

void Individual::PF_NEH(int e) {

    profileFitting(e);
    
    int gamma;
    if (parameters->numJobs > 20) {
        gamma = 25;
    } else {
        gamma = 21;
    }

    int menorJobMovido = parameters->numJobs - gamma + 1;
    int maiorJobMovido = parameters->numJobs - 1;

    vector <int> LJ;

    for (int i = parameters->numJobs - 1; i >= parameters->numJobs - gamma; i--) {
        LJ.push_back(chromosome[i]);

        solutionCost.zeroBlocks = solutionCost.zeroBlocks - idleOrBlocked[i];
        chromosome.erase(chromosome.begin() + i);
        idleOrBlocked[i] = 0;
    }

    int p;
    bool improved = false;

    for(int i = menorJobMovido + 1;i < chromosome.size() + 1; i++) {
        solutionCost.zeroBlocks = solutionCost.zeroBlocks - idleOrBlocked[i-1];
        idleOrBlocked[i-1] = 0;
        
        this->E[i][0] = this->E[i-1][1];
        for(int j = 1;j < parameters->numTools+1; j++) {
            if (j != parameters->numTools) {
                this->E[i][j] = max(this->E[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1], this->E[i-1][j+1]);
            } else {
                this->E[i][j] = this->E[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
            }

            idleOrBlocked[i-1] = idleOrBlocked[i-1] + this->E[i][j] - this->E[i-1][j] - (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
        }
        
        solutionCost.zeroBlocks = solutionCost.zeroBlocks + idleOrBlocked[i-1];
        // sol.jobsInfo[i-1].posicao = i - 1;
    }

    int jobInicio = chromosome.size() - 1 - (parameters->numJobs - maiorJobMovido - 1);
    int maquinaInicio, aux;
    
    aux = jobInicio;
    for(int i = chromosome.size() - aux;i < chromosome.size() + 1; i++) {
        this->Q[i][0] = this->Q[i-1][1];
        maquinaInicio = parameters->numTools-1;
        for(int j = 1;j < parameters->numTools+1; j++) {
            if (j != parameters->numTools) {
                this->Q[i][j] = max(this->Q[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio], this->Q[i-1][j+1]);
                maquinaInicio--;
            } else {
                this->Q[i][j] = this->Q[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio];
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
        
        for (int i = 0; i < chromosome.size() + 1; i++) {
                vAux = parameters->numTools;

                for (int j = 1; j <= 1; j++) {
                    guardaValor = max(this->Q[chromosome.size() - i][j] + (int)parameters->jobsToolsMatrix[LJ[k]][parameters->numTools - j], this->Q[chromosome.size() - i][j + 1]);
                    C[i][j] = this->E[i][vAux] + guardaValor;
                
                    if (C[i][j] > maximo) {
                        maximo = C[i][j];
                    }

                    vAux--;
                
                }

                for (int j = 2; j < parameters->numTools; j++) {
                    guardaValor = max(guardaValor + parameters->jobsToolsMatrix[LJ[k]][parameters->numTools - j], this->Q[chromosome.size() - i][j + 1]);
                    C[i][j] = this->E[i][vAux] + guardaValor;
                    
                    if (C[i][j] > maximo) {
                        maximo = C[i][j];
                    }

                    vAux--;
                
                }

                for (int j = parameters->numTools; j <= parameters->numTools; j++) {
                    C[i][j] = this->E[i][vAux] + guardaValor + (int)parameters->jobsToolsMatrix[LJ[k]][parameters->numTools - j];
                
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
    
        chromosome.insert(chromosome.begin() + p, LJ[k]);
        LJ.erase(LJ.begin() + k);

        for(int i = p + 1; i < chromosome.size() + 1; i++) {
            this->E[i][0] = this->E[i-1][1];
            solutionCost.zeroBlocks = solutionCost.zeroBlocks - idleOrBlocked[i-1];
            idleOrBlocked[i-1] = 0;
            for(int j = 1;j < parameters->numTools+1; j++) {
                if (j != parameters->numTools) {
                    this->E[i][j] = max(this->E[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1], this->E[i-1][j+1]);
                } else {
                    this->E[i][j] = this->E[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
                }

                idleOrBlocked[i-1] = idleOrBlocked[i-1] + this->E[i][j] - this->E[i-1][j] - (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
            }
            // sol.jobsInfo[i-1].posicao = i-1;
            solutionCost.zeroBlocks = solutionCost.zeroBlocks + idleOrBlocked[i-1];
            // sol.jobsInfo[i-1].tempo_sistema = this->E[i][parameters->numTools] - this->E[i][1];
        }
        
        jobInicio = p;
        for(int i = chromosome.size() - p;i < chromosome.size() + 1; i++) {
            this->Q[i][0] = this->Q[i-1][1];
            maquinaInicio = parameters->numTools-1;
            for(int j = 1;j < parameters->numTools+1; j++) {
                if (j != parameters->numTools) {
                    this->Q[i][j] = max(this->Q[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio], this->Q[i-1][j+1]);
                    maquinaInicio--;
                } else {
                    this->Q[i][j] = this->Q[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio];
                    maquinaInicio--;
                }
            }
            jobInicio--;  
        }

    }

     solutionCost.evaluation = this->E[parameters->numJobs][parameters->numTools];

    for (int i = 0; i < parameters->numJobs+1; i++) {
        delete[] C[i];
    }

    delete[] C;

    jobsInfo.resize(parameters->numJobs);
    caminhoCritico();

    if (this->E[parameters->numJobs][parameters->numTools] !=  solutionCost.evaluation || this->Q[parameters->numJobs][parameters->numTools] !=  solutionCost.evaluation) {
        cout << "PROBLEMA NA PF_NEH!!!" << endl;
        exit(1);
    }



}

void Individual::profileFittingAux(int elemento) {

    unsigned int ** C = new unsigned int * [parameters->numJobs + 1];
    for (int i = 0; i < parameters->numJobs + 1; i++) {
         C[i] = new unsigned int [parameters->numTools + 1];
    }

    for(int i = 0; i < parameters->numJobs + 1; i++) {
        this->E[i][0] = 0;
        this->Q[i][0] = 0;
        C[i][0] = 0;
    }

    for(int j = 0; j < parameters->numTools + 1; j++) {
        this->E[0][j] = 0;
        this->Q[0][j] = 0;
        C[0][j] = 0;
    }

    chromosome.push_back(elemento);
    edgesIndividuals.push_back(0);

    vector <jobMap> LC(parameters->numJobs);
    for (int i = 0; i < parameters->numJobs; i++) {
        if (i != elemento) {
            LC[i].indice = i;
        }
    }


    int elementosAlocados = 1;

    for(int i = 1;i <= elementosAlocados; i++) {
        idleOrBlocked[i-1] = 0;
        this->E[i][0] = this->E[i-1][1];
        for(int j = 1;j < parameters->numTools+1; j++) {
            if (j != parameters->numTools) {
                this->E[i][j] = max(this->E[i][j-1] + (int)parameters->jobsToolsMatrix[elemento][j-1], this->E[i-1][j+1]);
            } else {
                this->E[i][j] = this->E[i][j-1] + (int)parameters->jobsToolsMatrix[elemento][j-1];
            }
            idleOrBlocked[i-1] = idleOrBlocked[i-1] + this->E[i][j] - this->E[i-1][j] - (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
        }
        solutionCost.zeroBlocks = solutionCost.zeroBlocks + idleOrBlocked[i-1];
    }


    int idle_blocked, decisao; 
    while (elementosAlocados < parameters->numJobs) {
        idle_blocked = 99999;
        for (int i = 1; i < LC.size() + 1; i++) {
            LC[i-1].decisao = 0;
            for (int j = 1; j < parameters->numTools + 1; j++) {
                if (j != parameters->numTools) {
                    C[i][j] = max(max(this->E[elementosAlocados][j], C[i][j-1]) + (int)parameters->jobsToolsMatrix[LC[i-1].indice][j-1], this->E[elementosAlocados][j+1]);
                } else {
                    C[i][j] = C[i][j-1] + (int)parameters->jobsToolsMatrix[LC[i-1].indice][j-1];
                }

                LC[i-1].decisao += C[i][j] - this->E[elementosAlocados][j] - (int)parameters->jobsToolsMatrix[LC[i-1].indice][j-1];
                
            }

            // if (decisao < idle_blocked) {
            //     idle_blocked = decisao;
            //     elemento = i - 1;
            // }
            
        }
        
        sort(LC.begin(), LC.end(), comparaDecisao);
        int aux = static_cast<int>(ceil((LC.size() / 2) + 0.5));
        elemento = rand() % aux;

        chromosome.push_back(LC[elemento].indice);
        LC.erase(LC.begin() + elemento);

        edgesIndividuals.push_back(0);
        elementosAlocados++;

        for (int j = 1; j < parameters->numTools + 1; j++) {
            this->E[elementosAlocados][j] = C[elemento + 1][j];

            idleOrBlocked[elementosAlocados - 1] = idleOrBlocked[elementosAlocados - 1] + this->E[elementosAlocados][j] - this->E[elementosAlocados-1][j] - (int)parameters->jobsToolsMatrix[chromosome[elementosAlocados - 1]][j-1];
        }

        solutionCost.zeroBlocks = solutionCost.zeroBlocks + idleOrBlocked[elementosAlocados - 1];

        if (elementosAlocados == parameters->numJobs) {
            solutionCost.evaluation = C[elemento + 1][parameters->numTools];
        }

    }

    int jobInicio = parameters->numJobs - 1;
    int maquinaInicio;
    for(int i=1;i<parameters->numJobs+1;i++) {
        this->Q[i][0] = this->Q[i-1][1];
        maquinaInicio = parameters->numTools - 1;
        for(int j=1;j<parameters->numTools+1;j++) {
            
            if (j != parameters->numTools) {
                this->Q[i][j] = max(this->Q[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio], this->Q[i-1][j+1]);
                maquinaInicio--;
            } else {
                this->Q[i][j] = this->Q[i][j-1] + (int)parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio];
                maquinaInicio--;
            }
        }
        jobInicio--;  
    }

    for (int i = 0; i < parameters->numJobs+1; i++) {
        delete[] C[i];
    }

    delete[] C;

    jobsInfo.resize(parameters->numJobs);
    caminhoCritico();

    if (this->E[parameters->numJobs][parameters->numTools] != solutionCost.evaluation || this->Q[parameters->numJobs][parameters->numTools] != solutionCost.evaluation) {
        for (int i = 0; i < parameters->numJobs; i++) {
            cout << chromosome[i] << " ";
        }
        cout << endl;
        cout << this->E[parameters->numJobs][parameters->numTools] << " " << solutionCost.evaluation << " " << this->Q[parameters->numJobs][parameters->numTools] << endl;
        cout << "PROBLEMA NA PROFILE FITTING!!!" << endl << endl;
        exit(1);
    }

}

Individual::Individual(Parameters *parameters) : parameters(parameters) {

    age = 0;
    isFitnessComputed = false;

    vector<vector<unsigned int>> E(parameters->numJobs + 1, vector<unsigned int>(parameters->numTools + 1));
    vector<vector<unsigned int>> Q(parameters->numJobs + 1, vector<unsigned int>(parameters->numTools + 1));
    
    this->E = E;
    this->Q = Q;

    vector<unsigned int> IB(parameters->numJobs, 0);
    idleOrBlocked = IB;
    solutionCost.zeroBlocks = 0; 

    // Start individual as permutation of jobs
    for (unsigned int i = 0; i < parameters->numJobs; i++) {
        chromosome.push_back(i);
        edgesIndividuals.push_back(0);
    }

    unsigned int temp, jj;
    // shuffling
    for (unsigned int i = 0; i <= (unsigned int)chromosome.size() - 1; i++) {
        jj = i + rand() % ((unsigned int)chromosome.size() - i);
        temp = chromosome[i];
        chromosome[i] = chromosome[jj];
        chromosome[jj] = temp;
    }

    jobsInfo.resize(parameters->numJobs);

    solutionCost.evaluation = calcCost(-1);
    caminhoCritico();

    // for (int i = 0; i < parameters->numJobs; i++) {
    //     cout << chromosome[i] << " ";
    // }

    // cout << endl << endl;
    
    // for (int i = 0; i < parameters->numJobs; i++) {
    //     cout << jobsInfo[i].b << " " << jobsInfo[i].i1 << " " << jobsInfo[i].nBlocos << endl;
    // }
    // cout << endl;

    // for (int i = 0; i < parameters->numJobs; i++) {
    //     cout << blocos[jobsInfo[i].b].em << " " << blocos[jobsInfo[i].b].sm << " " << blocos[jobsInfo[i].b].indice << " " << blocos[jobsInfo[i].b].jf << " " << blocos[jobsInfo[i].b].ji << " " << blocos[jobsInfo[i].b].joined << " " << blocos[jobsInfo[i].b].tipo << endl;
    // }
    // getchar();
    //solutionCost.zeroBlocks = calcZeroBlocks();*/
}

Individual::~Individual() {}

unsigned int Individual::calcCost(int jump) {

    for(int i = 0; i < parameters->numJobs + 1; i++) {
        E[i][0] = 0;
        Q[i][0] = 0;
    }

    for(int j = 0; j < parameters->numTools + 1; j++) {
        E[0][j] = 0;
        Q[0][j] = 0;
    }
    
    solutionCost.zeroBlocks = 0; 

    for (int i = 1; i < parameters->numJobs + 1; i++) {
        E[i][0] = E[i-1][1];
        idleOrBlocked.push_back(0);
        for (int j = 1; j < parameters->numTools + 1; j++) {
            if (j != parameters->numTools) {
                E[i][j] = max(E[i][j-1] + parameters->jobsToolsMatrix[chromosome[i-1]][j-1], E[i-1][j+1]);
            } else {
                E[i][j] = E[i][j-1] + parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
            }
            idleOrBlocked[i-1] += E[i][j] - E[i-1][j] - parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
        }
        solutionCost.zeroBlocks += idleOrBlocked[i-1];
    }

    int jobInicio = parameters->numJobs - 1;
    int maquinaInicio;
    for(int i=1;i<parameters->numJobs+1;i++) {
        Q[i][0] = Q[i-1][1];
        maquinaInicio = parameters->numTools - 1;
        for(int j=1;j<parameters->numTools+1;j++) {
            
            if (j != parameters->numTools) {
                Q[i][j] = max(Q[i][j-1] + parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio], Q[i-1][j+1]);
                maquinaInicio--;
            } else {
                Q[i][j] = Q[i][j-1] + parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio];
                maquinaInicio--;
            }
        }
        jobInicio--;  
    }

/*    
    // This method implements the KTNS policy
    // (keep tools needed soonest)
    // It is a dynamic programming approach to calculate
    // the number of tool switches given a fixed order of
    // jobs in polynomial time

    // Fills L matrix (auxiliary structure used by KTNS)
    for (int i = (parameters->numJobs - 1); i >= 0; i--) {
        for (unsigned int j = 0; j < parameters->numTools; j++) {
            if (i != jump && parameters->jobsToolsMatrix[chromosome[i]][j] == 1) {
                parameters->L[j][i] = (unsigned int)i;
            }
            else if (i < (int)parameters->numJobs - 1) {
                parameters->L[j][i] = parameters->L[j][i + 1];
            }
            else {
                parameters->L[j][i] = parameters->numJobs;
            }
            parameters->used[j] = false;
        }
    }

    unsigned int switches = 0;
    unsigned int capacity = 0;
    unsigned int tool = 0;
    double minVal;

    // Fills auxiliary vector W_n that corresponds to the
    // matrix row at instant n
    for (unsigned int i = 0; i < parameters->numTools; i++) {
        if (parameters->L[i][0] == 0) {
            parameters->W_n[i] = 1;
            parameters->used[i] = true;
            capacity++;
        }
        else {
            parameters->W_n[i] = 0;
        }
    }

    while (capacity < parameters->maxCapacity) {
        minVal = numeric_limits<double>::infinity();
        for (unsigned int i = 0; i < parameters->numTools; i++) {
            if (!parameters->used[i] && (parameters->L[i][0] < minVal)) {
                tool = i;
                minVal = parameters->L[i][0];
            }
        }
        parameters->used[tool] = true;
        parameters->W_n[tool] = 1;
        capacity++;
    }

    parameters->loadedMatrix[0] = parameters->W_n;

    unsigned int maxVal;
    for (unsigned int n = 1; n < parameters->numJobs; n++) {
        for (unsigned int i = 0; i < parameters->numTools; i++) {
            if (parameters->W_n[i] != 1 && parameters->L[i][n] == n) {
                parameters->W_n[i] = 1;
                capacity++;
            }
        }
        while (capacity > parameters->maxCapacity) {
            maxVal = n;
            for (unsigned int i = 0; i < parameters->numTools; i++) {
                if (parameters->W_n[i] == 1 && parameters->L[i][n] > maxVal) {
                    tool = i;
                    maxVal = parameters->L[i][n];
                }
            }
            parameters->W_n[tool] = 0;
            capacity--;
            switches++;
        }
        parameters->loadedMatrix[n] = parameters->W_n;
    }*/

    return E[parameters->numJobs][parameters->numTools];
}

unsigned int Individual::updateCost(int p1, int p2) {

    for(int i = p1 + 1;i < chromosome.size() + 1; i++) {
        solutionCost.zeroBlocks -= idleOrBlocked[i-1];
        idleOrBlocked[i-1] = 0;
        for(int j = 1;j < parameters->numTools+1; j++) {
            E[i][0] = E[i-1][1];
            if (j != parameters->numTools) {
                E[i][j] = max(E[i][j-1] + parameters->jobsToolsMatrix[chromosome[i-1]][j-1], E[i-1][j+1]);
            } else {
                E[i][j] = E[i][j-1] + parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
            }

            idleOrBlocked[i-1] += E[i][j] - E[i-1][j] - parameters->jobsToolsMatrix[chromosome[i-1]][j-1];
        }
        solutionCost.zeroBlocks += idleOrBlocked[i-1];
    }
    
    
    int jobInicio = chromosome.size() - 1 - (parameters->numJobs - p2 - 1);
    int maquinaInicio;
    for(int i = chromosome.size() - jobInicio;i < chromosome.size() + 1; i++) {
        Q[i][0] = Q[i-1][1];
        maquinaInicio = parameters->numTools-1;
        for(int j = 1;j < parameters->numTools+1; j++) {
            if (j != parameters->numTools) {
                Q[i][j] = max(Q[i][j-1] + parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio], Q[i-1][j+1]);
                maquinaInicio--;
            } else {
                Q[i][j] = Q[i][j-1] + parameters->jobsToolsMatrix[chromosome[jobInicio]][maquinaInicio];
                maquinaInicio--;
            }
        }
        jobInicio--;  
    }

    return E[parameters->numJobs][parameters->numTools];

}

//unsigned int Individual::calcOneBlocks() {
//
//    unsigned int oneBlocks = 0;
//
//    for (unsigned int i = 0; i < parameters->numTools; i++) {
//        if (parameters->jobsToolsMatrix[chromosome[0]][i] == 1) {
//            oneBlocks++;
//        }
//        for (unsigned int j = 1; j < parameters->numJobs; j++) {
//            if (parameters->jobsToolsMatrix[chromosome[j - 1]][i] == 0 && parameters->jobsToolsMatrix[chromosome[j]][i] == 1) {
//                oneBlocks++;
//            }
//        }
//    }
//    return oneBlocks;
//}

double Individual::calcZeroBlocks() {
    // This method calculates the second (auxiliary) objective
    // used in the local search procedures:
    // The sum of the square roots of zero block sizes of the solution

    double zeroBlocks = 0;

    for (unsigned int i = 0; i < parameters->numTools; i++) {
        double sizeOfBlock = 0;
        for (unsigned int j = 1; j < parameters->numJobs; j++) {
            if ((parameters->loadedMatrix[j - 1][i] == 1 || sizeOfBlock > 0) && parameters->loadedMatrix[j][i] == 0) {
                sizeOfBlock++;
            }
            if (sizeOfBlock > 0 && parameters->loadedMatrix[j][i] == 1) {
                zeroBlocks += sqrt(sizeOfBlock);
                sizeOfBlock = 0;
            }
        }
        zeroBlocks += sqrt(sizeOfBlock);
    }

    return zeroBlocks;
}

void Individual::recopyIndividual(Individual *destination, Individual *source) {

    destination->chromosome = source->chromosome;
    destination->solutionCost.evaluation = source->solutionCost.evaluation;
    destination->solutionCost.zeroBlocks = source->solutionCost.zeroBlocks;
    destination->closest = source->closest;
    destination->isFitnessComputed = source->isFitnessComputed;
    destination->age = 0;
    destination->E = source->E;
    destination->Q = source->Q;
    destination->jobsInfo = source->jobsInfo;
    destination->blocos = source->blocos;
}

unsigned int Individual::distance(Individual *indiv) {
    // This method computed the edges distance between two solutions
    // If an element of solution 1 has a different neighbor (behind
    // or ahead) from the same element of solution 2, then the distance
    // is incremented by one

    unsigned int dist = 0;
    for (unsigned int i = 0; i < chromosome.size() - 2; i++) {
        edgesIndividuals[chromosome[i]] = (int) chromosome[i + 1];
    }
    for (unsigned int i = 0; i < indiv->chromosome.size() - 2; i++) {
        if (edgesIndividuals[indiv->chromosome[i]] != indiv->chromosome[i + 1] && edgesIndividuals[indiv->chromosome[i + 1]] != indiv->chromosome[i]) {
            dist++;
        }
    }
    return dist;
}

void Individual::addClose(Individual *indiv) {
    // Add an individual in the structure of proximity
    list<ProxData>::iterator it;
    ProxData data;
    data.individual = indiv;
    data.dist = distance(indiv);

    if (closest.empty())
        closest.push_back(data);
    else {
        it = closest.begin();
        while (it != closest.end() && it->dist < data.dist)
            ++it;
        closest.insert(it, data);
    }
}

void Individual::removeClose(Individual *indiv) {
    // Remove an individual in the structure of proximity
    list<ProxData>::iterator last = closest.end();
    for (list<ProxData>::iterator first = closest.begin(); first != last;)
        if (first->individual == indiv)
            first = closest.erase(first);
        else
            ++first;
}

double Individual::distToClosests(int n) {
    // Compute the average distance with the n close elements
    double result = 0;
    double compte = 0;
    list<ProxData>::iterator it = closest.begin();

    for (int i = 0; i < n && it != closest.end(); i++) {
        result += it->dist;
        compte += 1.0;
        ++it;
    }
    return result / compte;
}

void Individual::caminhoCritico() {

    Job job;
    Bloco bloco;

    //Inicia a construção do caminho crítico pelo último nó
    int j = parameters->numTools; 
    int i = parameters->numJobs;

    // cout << endl << j << " " << i;
    // getchar();

    blocos.clear();

    int cont = 0; //Variável para contar o número de blocos

    int subSeqAtual = 0; //Determina o tipo de sub-sequência atual. Pode ser: 0, 1, -1.

    // Variáveis para armazenar o último nó visitado
    int ultima_maquina = j;
    int ultimo_job = i; 

    // cout << endl << endl;

    // Inicializa o último bloco do caminho crítico, que termina no nó (n, m)
    bloco.indice = cont;
    bloco.ji = parameters->numJobs - 1;
    bloco.jf = parameters->numJobs - 1;
    bloco.em = parameters->numTools - 1;
    bloco.tipo = subSeqAtual;

    //Construção do caminho crítico
    while (i != 1 || j != 1) {

        // cout << i << " " << j << endl;

        // if (i < 1) {
        //     cout << i << endl;
        //     for (int l = 0; l < parameters->numJobs; l++) {
        //         cout << this->chromosome[l] << " ";
        //     }

        //     cout << endl << endl;

        //     for (int l = 0; l <= parameters->numJobs + 1; l++) {
        //         for (int m = 0; m <= parameters->numTools + 1; m++) {
        //             cout << this->E[l][m] << " ";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        //     exit(1);
        // }
       
        if (i != 1 && j != parameters->numTools && E[i][j] == E[i-1][j+1]) {
            i = i - 1;
            j = j + 1;

            // if (i == parameters->numJobs - 1) {
            //     cout << "a";
            //     getchar();
            // }

            if (subSeqAtual != -1) {

                bloco.sm = ultima_maquina - 1;
                bloco.ji = ultimo_job - 1;
                
                if (subSeqAtual == 1) {
                    //cout << "aaaaa";
                    //getchar();
                    
                    bloco.joined = true;
                    blocos.push_back(bloco);

                    jobsInfo[i].nBlocos = 2;
                    jobsInfo[i].i1 = cont;
                    jobsInfo[i].b = cont + 1;

                } else {
                    bloco.joined = false;
                    blocos.push_back(bloco);

                    jobsInfo[i].nBlocos = 3;
                    jobsInfo[i].i1 = cont; 
                    jobsInfo[i].b = cont;
                    //cout << "bbxcxcxcxc";
                    //getchar();

                }

                ++cont;

                ++bloco.indice;
                bloco.em = ultima_maquina - 1;
                bloco.jf = ultimo_job - 1;
                bloco.tipo = -1;

                subSeqAtual = -1;
                
            } else {
                jobsInfo[i].nBlocos = 1;
                jobsInfo[i].i1 = cont;
                jobsInfo[i].b = cont;
            }

        } else if (i != 1 && (E[i][j] == E[i-1][j] + (int)parameters->jobsToolsMatrix[chromosome[i-1]][j-1])) {
            i = i - 1;
            j = j;

            // if (i == parameters->numJobs - 1) {
            //     cout << "b";
            //     getchar();
            // }


            if (subSeqAtual != 1) {

                bloco.sm = ultima_maquina - 1;
                bloco.ji = ultimo_job - 1;
                
                if (subSeqAtual == -1) {

                    bloco.joined = true;
                    blocos.push_back(bloco);

                    jobsInfo[i].nBlocos = 2;
                    jobsInfo[i].i1 = cont;
                    jobsInfo[i].b = cont + 1;

                } else {
                    bloco.joined = false;
                    blocos.push_back(bloco);

                    jobsInfo[i].nBlocos = 3;
                    jobsInfo[i].i1 = cont;  
                    jobsInfo[i].b = cont;

                }

                ++cont;

                ++bloco.indice;
                bloco.em = ultima_maquina - 1;
                bloco.jf = ultimo_job - 1;
                bloco.tipo = 1;

                subSeqAtual = 1;

            } else {
                jobsInfo[i].nBlocos = 1;
                jobsInfo[i].i1 = cont;
                jobsInfo[i].b = cont;
            }

        } else {
            i = i;
            j = j - 1;

            // if (i == parameters->numJobs - 1) {
            //     cout << "c";
            //     getchar();
            // }


            if (subSeqAtual != 0) {

                bloco.sm = ultima_maquina - 1;
                bloco.ji = ultimo_job - 1;
                bloco.joined = false;

                ++cont;

                blocos.push_back(bloco);

                ++bloco.indice;
                bloco.em = ultima_maquina - 1;
                bloco.jf = ultimo_job - 1;
                bloco.tipo = 0;

                subSeqAtual = 0;

            }
        } 

        ultima_maquina = j;
        ultimo_job = i;
            
    }

    bloco.sm = ultima_maquina - 1;
    bloco.ji = ultimo_job - 1;

    if (bloco.tipo == 0) {
        jobsInfo[i-1].nBlocos = 2;
        jobsInfo[i-1].b = bloco.indice;
        jobsInfo[i-1].i1 = bloco.indice;

    } else {
        jobsInfo[i-1].nBlocos = 1;
        jobsInfo[i-1].b = bloco.indice;
        jobsInfo[i-1].i1 = bloco.indice;
    }

    blocos.push_back(bloco);

    // cout << "oi " << auxx << endl;
    // getchar();

}

bool comparaDecisao(jobMap j1, jobMap j2) {

    return j1.decisao < j2.decisao;

}