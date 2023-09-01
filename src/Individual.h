/**
 * @class Individual
 *
 * @brief Handle individual information
 *
 * Store individual (solution) information such as its chromosome,
 * cost, zero blocks and fitness.
 *
 * @author Jordana Mecler
 *
 * Contact: jmecler@inf.puc-rio.br
 *
 */
#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <vector>
#include <list>
#include <iostream>
#include "Parameters.h"
#include "LocalSearch.h"

using namespace std;

struct jobMap {

    int indice;
    int decisao;

};

bool comparaDecisao(jobMap j1, jobMap j2);

class LocalSearch;

/// Structure used to keep the cost of a solution
struct SolutionCost {

    // Cost of the solution
    unsigned int evaluation;

    // Sum of square root of size of zero blocks
    double zeroBlocks;

    // Constructor
    SolutionCost() {
        evaluation = 0;
        zeroBlocks = 0;
    }
};

/// Preliminary declaration, because proxData depends on Individual
class Individual;

/// Structure about the proximity of an Individual with regards to the others in the population
struct ProxData {
    // Individual
    Individual *individual;

    // Its distance to the others
    double dist;
};

struct Bloco {

    int indice; // O índice do bloco
    int tipo; // Tipo do bloco (GH = 1, ST = -1, T = 0)
    //int tamanho; // Tamanho do bloco
    int sm; // Máquina em que o bloco começa
    int em; // Máquina em que o bloco termina
    int ji; // Job que inicia o bloco
    int jf; // Job que termina o bloco
    bool joined;

};

struct Job {

    int nBlocos; // Quantidade de blocos que o job se encontra (<= 2)
    int b;  // Índice do bloco a ser considerado
    int i1; // Índice do 1º bloco que o job se encontra
    int i2; // Índice do último bloco que o job se encontra (-1 se estiver em apenas um bloco)
    // int idleOrBlocked = 0;
    // int tempo_sistema;
    // int posicao;
    //int indice; // Índice do job na sequência

};

class Individual
{
public:

    Parameters *parameters; ///< problem parameters

    int age; ///< age of individual

    double fitnessExt; ///< fitness of individual considering its cost and contribution to diversity

    float divRank; ///< rank of individual in terms of diversity

    float fitRank; ///< rank of individual in terms of cost

    SolutionCost solutionCost; ///< KTNS cost of solution

    vector<int> chromosome; ///< chromosome of individual stores a permutation of jobs
    vector< vector<unsigned int> > E; ///< departure times in direct order
    vector< vector<unsigned int> > Q; ///< departure times in reverse order
    vector<unsigned int> idleOrBlocked; ///< sum of idle and blocking times related to a job

    vector <Job> jobsInfo;
    vector <Bloco> blocos;

    void caminhoCritico(); 
    void random_swap();
    
    vector<int> edgesIndividuals; ///< temporary structure to compute distance

    bool isFitnessComputed; ///< if fitness has been computed

    /**
     * Calculate distance (edges distance) to an individual.
     *
     * @param indiv individual to be compared with current
     *
     * @return distance
     */
    unsigned int distance(Individual *indiv);

    list<ProxData> closest; ///< Population individuals ranked by proximity

    /**
     * Add individual to proximity structure.
     *
     * @param indiv individual to be added
     */
    void addClose(Individual *indiv);

    /**
     * Remove individual from proximity structure.
     *
     * @param indiv individual to be removed.
     */
    void removeClose(Individual *indiv);

    /**
     * Calculate distance to n closest individuals
     *
     * @param n number of individuals to be considered
     *
     * @return distance
     */
    double distToClosests(int n);

    LocalSearch *localSearch; ///< structure to perform local search

    /**
     * Calculates the cost of the solution using KTNS
     * (keep tools needed soonest).
     *
     * @param jump position of job to skip
     *
     * @return value of solution cost
     */
    unsigned int calcCost(int jump);
    unsigned int updateCost(int, int);
    
    /**
     * Calculate the auxiliary objective used in local search:
     * sum of the square roots of the sizes of zero blocks in the solution.
     *
     * @return value of auxiliary objective
     */
    double calcZeroBlocks();

    /**
     * Copy one individual attributes to another.
     *
     * @param destination destination individual
     * @param source source individual
     * jobsTools */
    void recopyIndividual(Individual *destination, Individual *source);

    /**
     * Constructor initializes a random individual.
     *
     * @param parameters problem parameters
     */
    explicit Individual(Parameters *parameters);
    explicit Individual(Parameters *parameters, int);

    void profileFitting(int);
    void profileFittingT(int);
    void profileFittingAux(int);
    void NEH(int);
    void PF_NEH(int);
    void PFT_NEH(int);
    void randomSolution();

    /// Destructor
    ~Individual();

};

#endif //INDIVIDUAL_H