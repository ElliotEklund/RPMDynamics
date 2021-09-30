//Author: Elliot Eklund
//Title: RPMDynamics
//Description: File containing functions needed for MC and dynamics algorithms

#include "includes/parameters.h"
#include "includes/functions.h"
#include <iostream>
#include <vector>

using namespace std;

//Computes instantaneous energy and energy estimator -------------------------------------------------------//
void getEnergy(double q[], double *energy, double *esti) {

    double tempEnergy = 0.0; //instantaneous energy
    double pot = 0.0; //potential energy
    double virial = 0.0; //virial energy estimator

    for (int i = 0; i<nBead; i++) {
        tempEnergy = tempEnergy + preFactor*(q[i]-q[(i+1) % nBead])*(q[i]-q[(i+1) % nBead]); //spring term
        pot = pot + potential(q[i]); //sum potential terms
	virial = virial + q[i]*dVdq(q[i]); 
    }
   
    * energy = tempEnergy + pot;
    * esti = 0.5*virial/nBead + pot/nBead;
}

//Potential energy function -------------------------------------------------------------------------------//
double potential(double q) {
//    return 0.5*q*q + 0.1*q*q*q + 0.01*q*q*q*q;
//	return 0.25*q*q*q*q;
    return 0.5*q*q;
}

//Derivative of potential energy function -----------------------------------------------------------------//
double dVdq (double q) {
//    return q + 0.3*q*q + 0.04*q*q*q;
//	return q*q*q;
    return q;
    }

//Velocity Verlety Algorithm ------------------------------------------------------------------------------//
void VV(double q[], double p[]) {

    int ii; //trajectory index
    double pHalfStep [nBead*trajs]; //momentum half step vector
    double halftime = 0.5*dt; //half time step
    
    for (int i=0; i<trajs; i++) {
        ii = nBead*i;
        
	//compute half step in momentum
        for (int j = 0; j<nBead; j++) {
            pHalfStep[ii+j] = p[ii+j] + halftime*(-mww*(2*q[ii+j] - q[ii+((j-1+nBead)%nBead)] -
                              q[ii+((j+1)%nBead)]) -dVdq(q[ii+j]));
        }
	//compute full step in position
        for (int j=0; j<nBead; j++) {
            q[ii+j] = q[ii+j] + dt*pHalfStep[ii+j]/m;
        }
	//compute final half step in momentum
        for (int j=0; j<nBead; j++) {
            p[ii+j] = pHalfStep[ii+j] + halftime*(-mww*(2*q[ii+j] - q[ii+((j-1+nBead)%nBead)] -
                    q[ii+((j+1)%nBead)]) -dVdq(q[ii+j]));
        }
    }
}

//Velocity Verlety Algorithm Vector Version-------------------------------------------------------------------//
void VV(vector <double> &q, vector <double> &p) {

   int ii; //trajectory index
   vector <double> pHalfStep;
   double halftime = 0.5*dt; //half time step

   for (int i=0; i<trajs; i++) {
        ii = nBead*i;
        for (int j = 0; j<nBead; j++) {
            pHalfStep.push_back(p[ii+j] + halftime*(-mww*(2*q[ii+j] - q[ii+((j-1+nBead)%nBead)] -
                              q[ii+((j+1)%nBead)]) -dVdq(q[ii+j])));
        }

        for (int j=0; j<nBead; j++) {
            q[ii+j] = q[ii+j] + dt*pHalfStep[ii+j]/m;
        }
        for (int j=0; j<nBead; j++) {
            p[ii+j] = pHalfStep[ii+j] + halftime*(-mww*(2*q[ii+j] - q[ii+((j-1+nBead)%nBead)] -
                    q[ii+((j+1)%nBead)]) -dVdq(q[ii+j]));
        }
    }
}



//Position Auto-Correlation function ----------------------------------------------------------------------//
void corrPos(double q0[], double qt[], double * corr) {

  double tCorr; //correlation function at time t
  double tCorrFinal = 0; //final correlation function at time t

  for (int i = 0; i<trajs; i++) {
      tCorr = 0;
      for (int j=0; j<nBead; j++) {
          tCorr = tCorr + qt[nBead*i + j];
      }
      tCorr = tCorr*q0[i]/nBead; //average over beads
      tCorrFinal = tCorrFinal + tCorr; //add current trajectory to total correlation function
  }
  *corr = tCorrFinal/trajs; //return average correlation function at time t
}

//Uniform distribution for nuclear step size --------------------------------------------------------------// 
double gen_rand(double min, double max, double random) {
    return (random * (max - min)) + min;
}
