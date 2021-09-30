//Author: Elliot Eklund
//Title: dataWorkup
//Description: File for working up data
//Last Updated: March 20, 2018

//Include statements for file handling, numerical mathematics, and timing
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include "mpi.h"
#include <vector>

#include "includes/nr3.h"
#include "includes/ran.h"
#include "includes/gamma.h"
#include "includes/deviates.h"

using namespace std;

int main(int argc, char **argv) {
    
    int num_procs = 10;
    int nBead = 32;
    double averageEnergy = 0;
    double standardDev = 0;
    
    vector <string> files; //vector of file names
    vector <double> energies;
    vector <double> totals;

    ifstream filez;
    
    for (int j=0; j<num_procs;j++){
        string file = "configOUTfile";
        string Result;
        ostringstream convert;
        convert << j;
        Result = convert.str();
        string fullName = file + Result;
        files.push_back(fullName);
    }

    for (int i =0; i<num_procs; i++) {
        filez.open(files[i].c_str());
        for (int j=0; j<nBead; j++) {
            double dummy;
            filez >> dummy;
        }
        double ENERGY;
        double MCTOTAL;
        filez >> ENERGY;
        filez >> MCTOTAL;
        energies.push_back(ENERGY);
        totals.push_back(MCTOTAL);
        
        filez.close();
    }
    
    
    for (int i=0; i<num_procs; i++) {
        energies[i] = energies[i]/totals[i];
        averageEnergy = averageEnergy + energies[i];
    }
    averageEnergy = averageEnergy/num_procs;
    
    for (int i=0; i<num_procs; i++) {
        standardDev = standardDev + (energies[i] - averageEnergy)*(energies[i] - averageEnergy);
    }
    standardDev = sqrt(standardDev/(num_procs-1));
    cout << "Average Energy: " << averageEnergy << endl;
    cout << "Standard Deviation: " << standardDev << endl;
    cout << "std/E: " << 100*standardDev/averageEnergy << endl;
    
    return 0;
}
