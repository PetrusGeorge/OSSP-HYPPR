/**
 * @file main.cpp
 *
 * @author Jordana Mecler
 *
 * Contact: jmecler@inf.puc-rio.br
 *
 */
#include <chrono>

#include "Individuo.h"
#include "Genetico.h"
#include "Populacao.h"

using namespace std;

int main(int argc, char* argv[]) {
	
	srand(time(NULL));

	Parameters *parameters;

	const clock_t startTime = clock();

	BuscaLocal *BL = new BuscaLocal(parameters);

	Populacao *p = new Populacao(parameters, BL);

	cout << p->getBestIndividual()->makespan << endl;

	for(int i =0; i < p->getBestIndividual()->chromosome.size(); i++){
			cout << p->getBestIndividual()->chromosome[i] << " ";
	} cout << endl;

	Genetico *GE = new Genetico(parameters, p, BL);

	// Individuo *i1 = GE->RR(p->getBestIndividual());

	// cout << i1->verifySequence() << endl;

	GE->evolve(parameters->numJobs * 20);

	cout << p->getBestIndividual()->makespan << endl;

	for(int i =0; i < p->getBestIndividual()->chromosome.size(); i++){
			cout << p->getBestIndividual()->chromosome[i] << " ";
	} cout << endl;

	cout << p->getBestIndividual()->verifySequence() << endl;

	totalTime = (float(clock() - startTime) / CLOCKS_PER_SEC);

	cout << totalTime << endl;
	delete parameters;

	return 0;
}