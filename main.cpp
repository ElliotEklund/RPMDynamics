//Author: Elliot Eklund
//Title: RPMDynamics
//Description: Main file for running RPMDynamics program
//Last Update: March 20, 2018

//Include statements for file handling, numerical mathematics, and timing
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include "mpi.h"
#include <string>

//Include headers specific to this program
#include "includes/parameters.h"
#include "includes/functions.h"

//Include files from Numerical Recipies text. They  must be included in the order shown. These
//files are necessary for random number generation
#include "includes/nr3.h"
#include "includes/ran.h"
#include "includes/gamma.h"
#include "includes/deviates.h"

//Global parameters
double const k = 1, hbar =1; //boltzmann & reduced planck's const [a.u.]
double m, T, BETA, BETAN, nucStep, dt, preFactor, mww; 
//mass [a.u.], temp [K], beta (1/kT), beta/nbead, MC step size, dynamics time interval, spring term prefactor
//prefactor for VV calculations
int nBead, trajs; //# of beads, # of dynamics trajectories

using namespace std;

int main(int argc, char **argv) {
    
    //Read in Files ----------------------------------------------------------------------------------------//
    double dynTime, dynSteps; //dynamics run time, steps need in VV algorithm 
    unsigned long long mcSteps, mcOutRate, dynDecorrRate;  
    //# of MC moves, when to output MC calculations, decorrelation "length" for sampling trajectories

    ifstream vars;
    vars.open("parameters.txt"); //open parameter.txt file
    vars >> m;
    vars >> T;
    vars >> nBead;

    vars >> mcSteps;
    vars >> mcOutRate;
    vars >> nucStep;

    vars >> trajs;
    vars >> dynTime;
    vars >> dynDecorrRate;
    vars >> dt;
    vars.close();  //close parameter.txt file
    
    int dynOption, EEOption; //run dynamics option, output energy estimator option
    int configIN, configOUT; //read in config file, output config file

    ifstream controls;
    controls.open("control.txt"); //open control file
    controls >> dynOption;
    controls >> EEOption;
    controls >> configIN;
    controls >> configOUT;
    controls.close(); //close control file
    
    //Main variables ---------------------------------------------------------------------------------------//
    int  mcMove, dynTrajsSampCount; //which bead to move, loop counter for  sampling dyn trajs
    unsigned long long mcTotal, accptRatio; //counts total MC moves, tracks accepted MC moves
    double stdev, avgEnerg, energi, energf, estimatori, estimatorf;
    //standard deviation for momentum sampling, average energy, initial and final energy, initial and final
    //energy estimators
    double * energiP = &energi; //pointer to initial energy
    double * energfP = &energf; //pointer to final energy
    double * estiiP = &estimatori; //pointer to initial energy estimator
    double * estifP = &estimatorf; //pointer to final energy estimator
    double qi [nBead]; //initial nuclear coordinates
    double qf [nBead]; //final nuclear coordinates
    double* qDynam = new double[nBead*trajs]; //coordinate vector for dynamics
    double* pDynam = new double[nBead*trajs]; //momentum vector for dynamics
    double  cent, duration; //centroid coordinate, code execution time
    
    //Main Variables ----------------------------------------------------------------------------------------//
    int ierr, num_procs, my_id, root_process;
    double sum, partial_sum;

    //Initialize variables
    mcTotal=0, accptRatio =0.0, dynTrajsSampCount=0;
    BETA=1.0/(k*T), BETAN=BETA/nBead;
    stdev=sqrt(m/BETAN);
    preFactor=0.5*m/(BETAN*BETAN*hbar*hbar);
    mww=m/(BETAN*BETAN*hbar*hbar);

    //Set up MPI ----------------------------------------------------------------------------------------------//
    root_process=0;
    MPI_Status status;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    ierr = MPI_Bcast(&mcSteps, 1, MPI_INT, root_process,MPI_COMM_WORLD);

    //Random Number Generator set up -----------------------------------------------------------------------//
    srand((unsigned)time(NULL)+my_id*num_procs); //Num Rec rng must be seeded with C++ standard rng
    Ran myRand(rand()); //position distribution from uniform distribution 
    Normaldev_BM mom(0,stdev,rand()); //momentum distribution from Gaussian(mu,sigma,seed)
    
    //Set up configurations file handler -------------------------------------------------------------------//
    
    vector <string> files; //vector of file names
    ifstream filez;
    if(configIN==1||configOUT==1){ //only fill files vector if we are writing or reading
        for (int j=0; j<num_procs;j++){
            string file = "configOUTfile";
            string Result;
            ostringstream convert;
            convert << my_id;
            Result = convert.str();
            string fullName = file + Result;
            files.push_back(fullName);
        }
    }
    else{}
    
    //Initialize variables for MC --------------------------------------------------------------------------//
    
    //read configOUTfile to nuclear bead positions
    if (configIN==1){
        for (int j=0; j<num_procs;j++){
            if(my_id==j){
                filez.open(files[j].c_str());
                for(int k=0; k<nBead; k++){
                    filez >> qi[k];
                    qf[k] = qi[k];
                }
                filez.close();
            }
            else {}
        }
    }
    
    else {
        for (int i=0; i<nBead; i++) {
            qi[i] = gen_rand(-nucStep,nucStep,myRand.doub()); //create distribution of initial coordinates
            qf[i] = qi[i];
        }
    }
    
    //Initialize Energy
    getEnergy(qf,energiP,estiiP);
    if(configIN==1){ //read in configOUTfile for energy and MCtotal variables
        for(int i=0;i<num_procs;i++){
            if(my_id==i){
                filez.open(files[i].c_str());
                for(int j=0;j<nBead;j++){ //skip down to energy
                    double dummy;
                    filez >> dummy;
                }
                double MCtotalPlace; //place holder to read in MCtotal
                filez >> avgEnerg; //read in avgEnergy
                filez >> MCtotalPlace; //read in MCtotal
                mcTotal = MCtotalPlace;
                filez.close();
            }
            else {}
        }
    }
    else{
        avgEnerg = estimatori;
    }
    
    //File Handling for MC ---------------------------------------------------------------------------------//
    ofstream mcCalculations;
    if (EEOption==1) {
        mcCalculations.open("mcCalculations.csv"); //create file for MC output
    }
    else {}

    //Performance clock stetup -----------------------------------------------------------------------------//
    clock_t start; //measure code performance 
    start = clock(); //begin timing
    
    
    //Run Monte Carlo --------------------------------------------------------------------------------------//
 
    for (unsigned long long j=0; j<mcSteps; j++) {
        for (int i=0; i<nBead; i++) {
            mcMove = myRand.int64() % nBead; //pick bead to move
            qf[mcMove] = qi[mcMove] + gen_rand(-nucStep,nucStep,myRand.doub());
        }
        getEnergy(qf, energfP, estifP);
        if (energf < energi) {
            for (int i=0; i<nBead; i++) {
                qi[i] = qf[i];
            }
            energi = energf;
            estimatori = estimatorf;
            accptRatio = accptRatio + 1;
        }
        
        else if (((myRand.int64()% 10000000)/10000000.0) <= exp(-BETAN*(energf-energi))) {
            for (int i=0; i<nBead; i++) {
                qi[i] = qf[i];
            }
            energi = energf;
            estimatori = estimatorf;
            accptRatio = accptRatio + 1;
        }
        
        else {
            for (int i=0; i<nBead; i++) {
                qf[i] = qi[i];
            }
            estimatorf = estimatori;
        }
        mcTotal = mcTotal + 1;
        avgEnerg = avgEnerg + estimatorf;
        
        if (mcTotal % mcOutRate == 0 && EEOption==1) {
            cent = 0;
            mcCalculations << mcTotal << ",";
            for (int j=0; j<nBead; j++) {
                cent = cent + qi[j];
            }
            mcCalculations << avgEnerg/mcTotal << "," << cent/nBead << endl;
        }
    }
    mcCalculations.close();
    cout << "Nuclear Acceptance Ratio: " << accptRatio/double(mcSteps) << endl << endl;
    cout << "Energy:" << avgEnerg/double(mcTotal) << endl;
    
    //Output final configurations from MC ------------------------------------------------------------------//
    if (configOUT==1){
        for (int j=0;j<num_procs;j++) {
            if (my_id==j) {
                ofstream filez;
                filez.open(files[j].c_str());
                for (int j=0; j<nBead; j++) {
                    filez << qf[j] << endl;
                }
                filez << avgEnerg << endl;
                filez << double(mcTotal) << endl;
                filez.close();
            }
            else {}
        }
    }
    else {}
    
    //Run Dynamics -----------------------------------------------------------------------------------------//   
    if (dynOption==1) {
    dynSteps = dynTime/dt; //initialize number of VV steps

    //Write samples to array for dynamics ------------------------------------------------------------------//	
        do {
            for (int j=0; j< ((dynDecorrRate + myRand.int64()%(dynDecorrRate/100)) - (dynDecorrRate/200)); j++) {
                for (int i=0; i<nBead; i++) {
                    mcMove = myRand.int64() % nBead; //pick bead to move
                    qf[mcMove] = qi[mcMove] + gen_rand(-nucStep,nucStep,myRand.doub());
                }
             
                getEnergy(qf ,energfP, estifP);
                if (energf < energi) {
                    for (int i=0; i<nBead; i++) {
                        qi[i] = qf[i];
                    }
                    energi = energf;
                }
                
                else if (((myRand.int64()% 10000000)/10000000.0) <= exp(-BETAN*(energf-energi))) {
                    for (int i=0; i<nBead; i++) {
                        qi[i] = qf[i];
                    }
                    energi = energf;
                }
                
                else {
                    for (int i=0; i<nBead; i++) {
                        qf[i] = qi[i];
                    }
                }
            }

            for (int l=0; l<nBead; l++) {
       
		qDynam [dynTrajsSampCount*nBead+l] = qi[l];
                pDynam [dynTrajsSampCount*nBead+l] = mom.dev();
            }
            dynTrajsSampCount = dynTrajsSampCount+1;
        } while(dynTrajsSampCount < trajs);
        
        //File handling for dynamics -----------------------------------------------------------------------//
        ofstream dynCalculations; //create object for dynamics calculations file
        dynCalculations.open("dynCalculations.csv"); //open dynamics calculations file
        
        vector <string> sideExpectFiles;
        for (int j=0; j<num_procs;j++){
            string file = "SideExpect";
            string Result;
            ostringstream convert;
            convert << my_id;
            Result = convert.str();
            string fullName = file + Result;
            sideExpectFiles.push_back(fullName);
        }
        ofstream outputFile[num_procs]; //array of ofstream objects
        for (int i = 0;i<num_procs; i++){
            outputFile[i].open(sideExpectFiles[i].c_str());
        }
        
        //Centroid set up ----------------------------------------------------------------------------------//
        double qCent[trajs], corrVal; //centroid, sum of centroids, current correlation fxn value
        double * corr = &corrVal; //pointer to current correlation fxn value
        
        //Initialize centroid
        for (int i=0; i<trajs; i++) {
            cent = 0;
            for (int j=0; j<nBead; j++) {
                cent = cent + qDynam[nBead*i + j];
            }
            qCent[i] = cent/nBead;
        }
        
        //Run Velocity Verlet and correlation fxn algorithms -----------------------------------------------//
        for (int i=0; i<dynSteps; i++) {
            VV(qDynam,pDynam); //call VV algorithm
            corrPos(qCent, qDynam, corr); //call correlation fxn algorithm
            //dynCalculations << i*dt << "," << corrVal << endl; //write outputs to file
            
            for(int j=0;j<num_procs;j++){
                if(my_id==j){
                    outputFile[j] << i*dt <<"," << corrVal << endl;
                }
                else{}
            }
        }
        
        dynCalculations.close(); //close dynamics calculations file

	delete[] qDynam; //delete pointer
	delete[] pDynam; //delete pointer
    }

    else {}

    //Close up MPI------------------------------------------------------------------------------------------//
    
    ierr = MPI_Reduce(&partial_sum, &sum, 1, MPI_FLOAT, MPI_SUM, root_process, MPI_COMM_WORLD); //bring all procs together
    ierr = MPI_Finalize();

    duration = ( clock() - start ) / (double) CLOCKS_PER_SEC; //end performance timing
    cout << "Time: "<< duration << endl;
    return 0;
}
