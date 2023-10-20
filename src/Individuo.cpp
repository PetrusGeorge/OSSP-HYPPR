#include "Individuo.h"

Individuo::Individuo(Parameters *parameters){
    this->parameters = parameters;
    this->makespan = BICH_MIH(parameters, endTimeOperations, chromosome);
    this->age = 0;
    this->fitRank = 0;
    this->fitnessExt = 0;
    this->divRank = 0;
}

Individuo::Individuo(Parameters *parameters, int c){
    this->parameters = parameters;
    this->makespan = 999999;
    this->chromosome = vector<int>(parameters->numJobs * parameters->numTools, 0);
    this->endTimeOperations = vector<int>(parameters->numJobs * parameters->numTools, 0);
    this->age = 0;
    this->fitRank = 0;
    this->fitnessExt = 0;
    this->divRank = 0;
}

Individuo::Individuo(){
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

        endTimeOperations[op-1] = M[indexJM.second] > J[indexJM.first] ? M[indexJM.second] : J[indexJM.first];
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

void Individuo::recopyIndividual(Individuo *indiv){
    this->makespan = indiv->makespan;
    this->chromosome = indiv->chromosome;
    this->fitRank = indiv->fitRank;
    this->closest = indiv->closest;
    this->fitnessExt = indiv->fitnessExt;
    this->divRank = indiv->divRank;
    this->endTimeOperations = indiv->endTimeOperations;
    this->age = 0;
}


void Individuo::addClose(Individuo *indiv) {
    // Add an individual in the structure of proximity
    list<Proximity>::iterator it;
    Proximity data;
    data.individual = indiv;
    data.distance = calculateDistance(indiv);

    if (closest.empty())
        closest.push_back(data);
    else {
        it = closest.begin();
        while (it != closest.end() && it->distance < data.distance)
            ++it;
        closest.insert(it, data);
    }
}

void Individuo::removeClose(Individuo *indiv) {
    // Remove an individual in the structure of proximity
    list<Proximity>::iterator last = closest.end();
    for (list<Proximity>::iterator first = closest.begin(); first != last;)
        if (first->individual == indiv)
            first = closest.erase(first);
        else
            ++first;
}

double Individuo::distToClosests(int n) {
    // Compute the average distance with the n close elements
    double result = 0;
    double compte = 0;
    list<Proximity>::iterator it = closest.begin();

    for (int i = 0; i < n && it != closest.end(); i++) {
        result += it->distance;
        compte += 1.0;
        ++it;
    }
    return result / compte;
}

int Individuo::calculateDistance(Individuo *individual){

    int distance = 0;
    
    vector<unsigned> nextOperaration(chromosome.size(), 0);

    for(int i = 0; i < chromosome.size() - 2; i++){
        nextOperaration[chromosome[i] - 1] = chromosome[i + 1];
    }

    for (unsigned int i = 1; i < individual->chromosome.size() - 2; i++) {
        if (nextOperaration[individual->chromosome[i] - 1] != individual->chromosome[i + 1] && nextOperaration[individual->chromosome[i - 1] - 1] != individual->chromosome[i]) {
            distance++;
        }
    }

    return distance;
}

bool Individuo::verifySequence(){
    vector<int> machines;

    for(int i = 0; i < chromosome.size(); i++){
        machines = vector<int>(parameters->numTools, 0);

        for(int j = 0; j < parameters->numTools; j++){
            int m = decryptJobMachineIndex(chromosome[i+j], parameters->numJobs).second;
            //cout << m << endl;
            if(machines[m] == 1){
                return false;
            }
            machines[m] = 1;
        }

        i = i + parameters->numTools - 1;
    }

    return true;

}