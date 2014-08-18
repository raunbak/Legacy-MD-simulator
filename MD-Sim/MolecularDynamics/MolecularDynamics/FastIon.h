//FastIon.h

#ifndef FASTION_H_
#define FASTION_H_


#include "integrator.h"
#include <cstring>

class FastEnsemble;

class FastIon
{
public:
	double m;	  	// Ion mass
	double *Pos; 	// pointer to array of size: 3 but memory is not allocated before member function initialize
	double *Vel; 	// pointer to array of size: 3 but memory is not allocated before member function initialize
	double VxSec;
	double VySec;
	double VzSec;
	double VrSec;
	double Vavg[3];
	double Vsec;
//public:
	// member functions
	void CleanUpIon();
	double ReturnVzSec();
	double ReturnVxSec();
	double ReturnVySec();
	void initialize(int mass); //mass should always be in u. allocate memory for Pos and Vel.
	void SetPosition(int dim, double Val); // set position
	void SetVelocity(int dim, double Val); // set velocity
	double GetMass();
	double GetVsec();
	double Position(int dim); // returning position
	double Velocity(int dim); // returning velocity
	double Velocity(); // return norm of velocity
	double Ekin(); // return kinetic energy
	// friend functions
	friend void LeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend);
	friend void CudaLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend);
	friend void TestCudaLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend);
	friend void DynamicTemperatureLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend); // Testing temperature control as in articles
	friend double DistanceSquar(FastIon & ion1, FastIon & ion2);
	friend double Distance(FastIon & ion1, FastIon & ion2);
	friend double Ffriction(FastEnsemble & ensemble, int N, int dim);
	friend double Ftrap(FastEnsemble & ensemble, int N, int TimeStep, int dim, double Vrf, double Vend);
	friend double Fcoulumb(FastEnsemble & ensemble, int N, int dim);

};

#endif
