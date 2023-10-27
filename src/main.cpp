/**
 * @file main.cpp
 *
 * @author Jordana Mecler
 *
 * Contact: jmecler@inf.puc-rio.br
 *
 */
#include <chrono>

#include "Individual.h"
#include "Genetic.h"
#include "Population.h"

using namespace std;

int main(int argc, char* argv[]) {
	
	srand(time(NULL));

	Parameters *parameters;

	try{
		parameters = new Parameters(argv[1]);
	} catch(string error){
		cout << error << endl;
		return 0;
	}

	const clock_t startTime = clock();

	LocalSearch *BL = new LocalSearch(parameters);

	Population *p = new Population(parameters, BL);

	cout << p->getBestIndividual()->makespan << endl;

	for(int i =0; i < p->getBestIndividual()->chromosome.size(); i++){
			cout << p->getBestIndividual()->chromosome[i] << " ";
	} cout << endl;

	Genetic *GE = new Genetic(parameters, p, BL);

	// Individual *i1 = GE->RR(p->getBestIndividual());

	// cout << i1->verifySequence() << endl;

	GE->evolve(parameters->numJobs * 20);

	cout << p->getBestIndividual()->makespan << endl;

	for(int i =0; i < p->getBestIndividual()->chromosome.size(); i++){
			cout << p->getBestIndividual()->chromosome[i] << " ";
	} cout << endl;

	cout << p->getBestIndividual()->verifySequence() << endl;

	double totalTime = (float(clock() - startTime) / CLOCKS_PER_SEC);

	cout << totalTime << endl;
	delete parameters;

	return 0;
}