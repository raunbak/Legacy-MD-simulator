//Ion.h

#ifndef ION_H_
#define ION_H_


#include "integrator.h"
#include <cstring>

class Ensemble;

class Ion
{
private:
	double m;	  	// Ion mass
	double **Pos; 	// pointer to array (dim,TimeStep) of size: TimeSteps x 3 but memory is not allocated before member function initialize
	double **Vel; 	// pointer to array (dim,TimeStep) of size: TimeSteps x 3 but memory is not allocated before member function initialize
	double VzSec;
public:
	// constructor

	// member functions
	void initialize(int TimeSteps); // allocate memory for Pos and Vel.
	void SetPosition(int dim, int TimeStep, double Val); // set position
	void SetVelocity(int dim, int TimeStep, double Val); // set velocity
	double GetMass();
	double Position(int dim, int TimeStep); // returning position
	double Velocity(int dim, int TimeStep); // returning velocity
	double Velocity(int TimeStep); // return norm of velocity
	double Ekin(int TimeStep); // return kinetic energy

	// friend functions
	friend void LeFrogintegrator(Ensemble & ensemble, int TimeSteps, double Vrf, double Vend);
	friend double distance(Ion & ion1, Ion & ion2, int TimeStep);
	friend double Ffriction(Ensemble & ensemble, int N, int TimeStep, int dim);
	friend double Ftrap(Ensemble & ensemble, int N, int TimeStep, int dim, double Vrf, double Vend);
	friend double Fcoulumb(Ensemble & ensemble, int N, int TimeStep, int dim);
	friend double FcoulumbNear(Ensemble & ensemble, int N, int TimeStep, int dim, int ****CellIndex, double h, int Ngrid);
};

#endif
