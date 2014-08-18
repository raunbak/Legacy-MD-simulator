//FastEnsemble.h

#ifndef FASTENSEMBLE_H_
#define FASTENSEMBLE_H_

#include <cmath>
#include <cstring>
#include "integrator.h"

// To save files
#include <iostream>
#include <fstream>

class FastIon; // dealing with circular dependens same in Ion header.


class FastEnsemble
{
public:
	int NumberOfIons;
	FastIon *ions; // pointer to array of ion objects
	double *Tsec; // array with all temperatures - needs to be initialized
	double *Tsecrad;
	double *Tsecz;
	double *Trms; // array with all temperatures - needs to be initialized
	double *Trmsz;
	long int ***histogram; // pointer to 3d histogram array of size ?
	long double ***VelHistogram;
	long int ***CountHistogram;
	double VzSecRMS; // RMS secular z-velocity of ions - used to rescale velocity distribution
	double VrSecRMS; // RMS secular radial-velocity of ions - used to rescale velocity distribution
	double SteadyStateTemperature;
	double ActualAvgTemperature;
	double ActualTemperatureSTD;
	double Radius;
	double Length;

	// For rescaling
	int NumberOfRings;
	double* VzSecRMS_in_ring;
	double* VrSecRMS_in_ring;
	int* NumberOfIonsInRadius;
	int* RadiiOfIons;
	
	double* Vz_IonArray;
	double* Vx_IonArray;
	double* Vy_IonArray;

	bool ReachedTempArea;

//public:
	// constructors
	FastEnsemble(int m1, int n1, int m2, int n2);
	// member function
	double GetTrms(int TimeStep);
	double GetTrmsz(int TimeStep);
	double GetTsecular(int TimeStep);
	double GetTsecularz(int TimeStep);
	double GetTsecularrad(int TimeStep);
	double getRho0(double Vrf);
	void InitialiseTemperatureArrays(int TimeSteps);
	void FreeTemperatureArrays();
	void CleanUpEnsemble();
	void PrintSecVel();
	void SetActualTemperature(double val);
	void SetActualTemperatureSTD(double val);
	double GetActualTemperature();
	double GetActualTemperatureSTD();
	double GetCurrentTemperature();
	void RescaleVelocityDistribution(); // rescaling velocity distribution
	void RescaleVelocityDistribution(int TimeStep); // rescaling velocity distribution
	void VelocityKick(int TimeStep); // Velocity Kick
	void InitialiseHistogram();
	void InitialiseVelocityHistogram();
	void InitialiseCountHistogram();
	void UpdateHistogram();
	void UpdateVelocityHistogram();
	void UpdateCountHistogram();
	double ReturnHist(int i, int j, int k); //returning value of bin (i,j,k) in histogram
	double ReturnVelHist(int i, int j, int k); //returning value of bin (i,j,k) in histogram
	double ReturnCountHist(int i, int j, int k);
	void SetSteadyStateTemperature(double Val);
	int GetNumberOfIons();
	double Mass(int N); // returning mass of ion N
	double Position(int dim, int N); // returning position
	double Velocity(int dim, int N); // returning velocity
	void CrystalGenerator(); // set ions initial positions in grid and set initial velocities - under construction
	void CrystalGenerator(double Vrf, double Vend); // set ions initial positions in grid and set initial velocities using the plasma model in bcc structure - under construction
	void CrystalGeneratorForPseudoPotential(double Vrf, double Vend);
	double Ekin(); // return total kinetic energy of crystal for given time step
	double Ttot(); // return total temperature of crystal for given time step
	void SavePositionToFile();
	

	void MyUpdateVelocityHistogram();
	void SaveIonDataToFile();
	void RescaleVelocityDistributionRadial();
	void RescaleVelocityDistributionIndivialIons();
	void RescaleVelocityDumb(double Vrmsz, double Vzrms);
	void RescaleVelocityXYZ(double Total_V_x_rms, double Total_V_y_rms,double Total_V_z_rms);
	// friend function
	friend double DistanceSquar(FastIon & ion1, FastIon & ion2);
	friend double Distance(FastIon & ion1, FastIon & ion2);
	friend double Ftrap(FastEnsemble & ensemble, int N, int TimeStep, int dim, double Vrf, double Vend);
	friend double Fcoulumb(FastEnsemble & ensemble, int N, int dim);
	friend double Ffriction(FastEnsemble & ensemble, int N, int dim);
	friend void LeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend);
	friend void CudaLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend);
	friend void TestCudaLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend);
	friend void DynamicTemperatureLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend); // Testing temperature control as in articles


};

double asinh(double x); // not in cmath lib
double AspectRatioFunction(double alpha); // from peter herskind thesis
double CalculateAspectRatio(double wratiosquare); // returning aspect ratio by solving aspectratioequation

#endif

