//Ion.cpp

#include "stdafx.h" // Pre-compiled header.
//Ion.cpp

#include <cmath>
#include <cstdlib>
#include "Ion.h"
#include "constants.h"

using namespace std;

// Constructer


// Member functions

void Ion::initialize(int TimeSteps) // setting mass and allocating memory
{
	m = 40*u2kg;
	VzSec = 0; // setting secular temperature to zero

	Pos = new double * [3];
		for(int i =0; i < 3 ; i++)
			Pos[i]= new double [TimeSteps];


	Vel = new double * [3];
		for(int i =0; i < 3 ; i++)
			Vel[i]= new double [TimeSteps];

}

void Ion::SetPosition(int dim, int TimeStep, double Val) // set position
{
	Pos[dim][TimeStep] = Val;
}

void Ion::SetVelocity(int dim, int TimeStep, double Val) // set position
{
	Vel[dim][TimeStep] = Val;
}

double Ion::GetMass()
{
	return m;
}

double Ion::Position(int dim, int TimeStep) // change name to get at some point
{
	return Pos[dim][TimeStep];
}


double Ion::Velocity(int dim, int TimeStep)
{
	return Vel[dim][TimeStep];
}

double Ion::Velocity(int TimeStep)
{
	return sqrt(pow(Vel[0][TimeStep],2) + pow(Vel[1][TimeStep],2) + pow(Vel[2][TimeStep],2));
}

double Ion::Ekin(int TimeStep)
{
	return 0.5*m*pow(pow(Vel[0][TimeStep],2) + pow(Vel[1][TimeStep],2) + pow(Vel[2][TimeStep],2),2);
}


// Friend functions
double distance(Ion & ion1, Ion & ion2, int TimeStep)
{
	return sqrt(pow(ion1.Pos[0][TimeStep]-ion2.Pos[0][TimeStep],2) + pow(ion1.Pos[1][TimeStep]-ion2.Pos[1][TimeStep],2) + pow(ion1.Pos[2][TimeStep]-ion2.Pos[2][TimeStep],2));
}
