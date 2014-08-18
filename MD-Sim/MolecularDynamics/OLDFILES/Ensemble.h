//Ensemble.h

#ifndef ENSEMBLE_H_
#define ENSEMBLE_H_

#include <cmath>
#include <cstring>
#include "integrator.h"

class Ion; // dealing with circular dependens same in Ion header.


class Ensemble
{
private:
	int NumberOfIons;
	Ion *ions; // pointer to array of ion objects
	double ***histogram; // pointer to 3d histogram array of size ?
	double VzSecRMS; // RMS secular z-velocity of ions - used to rescale velocity distribution
	double SteadyStateTemperature;
public:
	// constructors
	Ensemble(int N, int TimeSteps);
	// member functions
	int GetNumberOfIons();
	void RescaleVelocityDistribution(int TimeStep); // rescaling velocity distribution
	void SetSteadyStateTemperature(double Val);
	double Mass(int N); // returning mass of ion N
	double Position(int dim, int N, int TimeStep); // returning position
	double Velocity(int dim, int N, int TimeStep); // returning velocity
	void CrystalGenerator(); // set ions initial positions in grid and set initial velocities - under construction
	void CrystalGenerator(double Vrf, double Vend); // set ions initial positions in grid and set initial velocities using the plasma model in bcc structure - under construction
	double Ekin(int TimeStep); // return total kinetic energy of crystal for given time step
	double Ttot(int TimeStep); // return total temperature of crystal for given time step
	double *Ttot(int Tstart, int Tend); // Returns array with total temperature from time step Tstart to Tend
	double *Tsec(int Tstart, int Tend); // Returns array with secular temperature from time step Tstart to Tend - under construction
	double Tsec(int TimeStep); // Returns Tsec calculated from last RF periode - so for the first periode it just returns 0;
	void SavePositionToFile(); // Writing out position to .xyz file which can be read by vmd
	void SavePositionToFile(int Tend); //write out all time steps
	// friend functions
	friend void LeFrogintegrator(Ensemble & ensemble, int TimeSteps, double Vrf, double Vend);
	friend double distance(Ion & ion1, Ion & ion2);
	friend double Ffriction(Ensemble & ensemble, int N, int TimeStep, int dim);
	friend double Ftrap(Ensemble & ensemble, int N, int TimeStep, int dim, double Vrf, double Vend);
	friend double Fcoulumb(Ensemble & ensemble, int N, int TimeStep, int dim);
};

double asinh(double x); // not in cmath lib
double AspectRatioFunction(double alpha); // from peter herskind thesis
double CalculateAspectRatio(double wratiosquare); // returning aspect ratio by solving aspectratioequation




#endif

