#include "Parameters.h"

void Parameters::setMethodParams() {
    populationSize = 20;
    maxPopulationSize = 40;
    numberElite = 10;
    numberCloseIndividuals = 3;
    maxDiversify = 1000;
    terminate = false;

    this->positionsOffspring = vector<bool>(numJobs*numTools, false);

    if (numJobs != 500) {
        nViz = 5;
    } else {
        nViz = 4;
    }
}