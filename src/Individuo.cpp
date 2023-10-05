#include "Individuo.h"

Individuo::Individuo(Parameters *parameters){
    this->parameters = parameters;
    this->makespan = BICH_MIH(parameters, endTimeOperations, chromosome);
}

unsigned int Individuo::calcMakespan(){
    vector<int> M(parameters->numTools, 0);//Machines acumulated time
    vector<int> J(parameters->numJobs, 0);//Jobs acumulated time
    vector<int> U = chromosome;
    //U is the sequence of operations
    while(!U.empty()){
     
        int op = U[0];
        U.erase(U.begin());
        
        pair<unsigned, unsigned> indexJM = decryptJobMachineIndex(op, parameters->numJobs);

        int acumulatedJ = (M[indexJM.second] - J[indexJM.first]) > EPSILON ? (M[indexJM.second] - J[indexJM.first]) : 0;
        int acumulatedM = (J[indexJM.first] - M[indexJM.second]) > EPSILON ? J[indexJM.first] - M[indexJM.second] : 0;

        J[indexJM.first] += acumulatedJ + parameters->jobsToolsMatrix[indexJM.first][indexJM.second];
        M[indexJM.second] += acumulatedM + parameters->jobsToolsMatrix[indexJM.first][indexJM.second];
    }

    return *max_element(M.begin(), M.end()) ;
}

// unsigned int Individuo::distance(Individuo *indiv) {
//     // This method computed the edges distance between two solutions
//     // If an element of solution 1 has a different neighbor (behind
//     // or ahead) from the same element of solution 2, then the distance
//     // is incremented by one

//     unsigned int dist = 0;
//     for (unsigned int i = 0; i < chromosome.size() - 2; i++) {
//         edgesIndividuals[chromosome[i]] = (int) chromosome[i + 1];
//     }
//     for (unsigned int i = 0; i < indiv->chromosome.size() - 2; i++) {
//         if (edgesIndividuals[indiv->chromosome[i]] != indiv->chromosome[i + 1] && edgesIndividuals[indiv->chromosome[i + 1]] != indiv->chromosome[i]) {
//             dist++;
//         }
//     }
//     return dist;
// }

// void Individual::addClose(Individual *indiv) {
//     // Add an individual in the structure of proximity
//     list<ProxData>::iterator it;
//     ProxData data;
//     data.individual = indiv;
//     data.dist = distance(indiv);

//     if (closest.empty())
//         closest.push_back(data);
//     else {
//         it = closest.begin();
//         while (it != closest.end() && it->dist < data.dist)
//             ++it;
//         closest.insert(it, data);
//     }
// }

// void Individual::removeClose(Individual *indiv) {
//     // Remove an individual in the structure of proximity
//     list<ProxData>::iterator last = closest.end();
//     for (list<ProxData>::iterator first = closest.begin(); first != last;)
//         if (first->individual == indiv)
//             first = closest.erase(first);
//         else
//             ++first;
// }

// double Individuo::distToClosests(int n) {
//     // Compute the average distance with the n close elements
//     double result = 0;
//     double compte = 0;
//     list<ProxData>::iterator it = closest.begin();

//     for (int i = 0; i < n && it != closest.end(); i++) {
//         result += it->dist;
//         compte += 1.0;
//         ++it;
//     }
//     return result / compte;
// }