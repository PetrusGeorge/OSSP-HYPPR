#include "Parameters.h"

Parameters::Parameters(char *argv){
    this->instancePath = string(argv);

    this->readInstance();
    this->setMethodParams();
}

void Parameters::readInstance(){
    ifstream instanceFile(this->instancePath.c_str());

    if(!instanceFile.is_open()){
        stringstream errorStream;
        errorStream << "Unable to open file " << this->instancePath;
        
        throw(errorStream.str());
    }

    string line;
    bool firstLine = true;
    while(getline(instanceFile, line)){

        stringstream lineStream(line);
        std::vector<unsigned int> lineVector;

        if(firstLine){
            firstLine = false;
            lineStream >> this->numJobs;
            lineStream >> this->numTools;
            continue;
        }

        unsigned int value;
        
        while(lineStream >> value){

            lineVector.push_back(value);
        }

        this->jobsToolsMatrix.push_back(lineVector);
    }
}

void Parameters::setAlpha() {
    // Calculate alpha parameter
    switch (this->numTools){
        case 3:
        case 4:
        case 15:
            this->alpha = 0.89;
            break;
        case 5:
            this->alpha = 0.59;
            break;
        case 6:
            this->alpha = 0.11;
            break;
        case 7:
        case 8:
        case 20:
            this->alpha = 0.22;
            break;
        case 9:
        case 10:
            this->alpha = 0.67;
            break;
        default:
            this->alpha = 0.1;
    }
}

void Parameters::setMethodParams() {
    this->populationSize = 20;
    this->maxPopulationSize = 40;
    this->numberElite = 10;
    this->numberCloseIndividuals = 3;
    this->maxDiversify = 1000;
    this->terminate = false;

    this->positionsOffspring = vector<bool>(this->numJobs*this->numTools, false);

    if (this->numJobs != 500) {
        this->nIterNeighborhood = 5;
    } else {
        this->nIterNeighborhood = 4;
    }

    setAlpha();
}