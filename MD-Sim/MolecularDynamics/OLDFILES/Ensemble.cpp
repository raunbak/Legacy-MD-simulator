//Ensemble.cpp
#include "stdafx.h" // Pre-compiled header.

#include <cmath>
#include <cstdlib>
#include "Ensemble.h"
#include "Ion.h"
#include "constants.h"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>


using namespace std;

// Constructor
Ensemble::Ensemble(int N, int TimeSteps)
{
	NumberOfIons=N;
	ions = new Ion [NumberOfIons];

	for(int n = 0; n < NumberOfIons; n++)
		ions[n].initialize(TimeSteps); // allocating memory for the ions pos and vel arrays

	// remember to delete after testing
	if (N==1)
			ions[0].SetPosition(0,0,0);

}
// Member functions
int Ensemble::GetNumberOfIons()
{
	return NumberOfIons;
}

void Ensemble::RescaleVelocityDistribution(int TimeStep)
{
	double CurrentTemp;
	// calculate current temperature
	CurrentTemp = pow(VzSecRMS,2)*Mass(0)/Kb;
	cout << "Current Secular z-Temperature is " << CurrentTemp*1000 << " mK \n";

	// rescale velocity distribution
	double a = sqrt(SteadyStateTemperature/CurrentTemp);
	cout << "Rescaling factor is " << a << '\n';

	if (a > 1.1)
	{
		a = 1.1;
		cout << "To big temperature jump rescaling factor is set to " << a << '\n';
	}
	for(int n = 0; n < NumberOfIons; n++)
	{
		ions[n].SetVelocity(0, TimeStep, ions[n].Velocity(0, TimeStep)*a);
		ions[n].SetVelocity(1, TimeStep, ions[n].Velocity(1, TimeStep)*a);
		ions[n].SetVelocity(2, TimeStep, ions[n].Velocity(2, TimeStep)*a);
	}
}

void Ensemble::SetSteadyStateTemperature(double Val)
{
	SteadyStateTemperature = Val;
}

double Ensemble::Mass(int N)
{
	return ions[N].GetMass();
}

void Ensemble::CrystalGenerator()
{
	for(int N = 0; N < NumberOfIons ; N++)
	{
		// setting initial positions
		ions[N].SetPosition(0, 0, pow((double) -1,(double) N)*1e-7);
		ions[N].SetPosition(1, 0, pow((double) -1,(double) N+1)*1e-7);
		ions[N].SetPosition(2, 0, GridSpacing*(N+1) - 0.5*GridSpacing*(NumberOfIons-1)); // making a string of ions

		// setting initial velocities
		double Vmean = sqrt(Kb*Tinitial/ions[N].GetMass());
		ions[N].SetVelocity(0, 0, 0); //(rand()/RAND_MAX-0.5)*Vmean
		ions[N].SetVelocity(1, 0, Vmean);
		ions[N].SetVelocity(2, 0, Vmean);
	}
}

void Ensemble::CrystalGenerator(double Vrf, double Vend)
{
	// Calculating pseudo trap frequencies
	double wz2=2*eta*e*Vend/pow(z0,2)/Mass(0); // ERROR in formula from Magnus' thesis! he means potential and not electric potential.
	double wr2=pow(e*Vrf/Mass(0)/OmegaRF,2)/2/pow(r0,4)-eta*e*Vend/Mass(0)/pow(z0,2);

	// Solving plasma model aspect-ratio equation
	double alpha = CalculateAspectRatio(wz2/wr2);
	//double alpha = CalculateAspectRatio(2);
	cout << "aspect ratio is " << alpha << '\n';
	/*cout << "out putting ratio(alpha)\n";
	for (double ratio = 0.01; ratio < 2; ratio+=0.01)
		if (ratio != 1)
			cout << CalculateAspectRatio(ratio*ratio) << ';' << '\n';
*/
	// Calculating number density
	double rho0 = eps0*Vrf*Vrf/(Mass(0)*pow(r0,4)*OmegaRF*OmegaRF);

	// Calculating volume of crystal
	double V = ((double) NumberOfIons)/rho0;
	cout << "volume of crystal is " << V << '\n';

	// Calculating length and radius of crystal
	double L = pow(6*V/(PI*alpha*alpha), ((double) 1) /  ((double) 3));
	bool CrystalNotDone = true;
	while(CrystalNotDone)
	{
		double R = alpha*L/2;

		cout << "radius is " << R << '\n';
		cout << "Length is " << L << '\n';
		cout << "2R/L=" << 2*R/L << '\n';

		// unit cell cube length for bcc structure
		double a = pow( 2 / rho0, ((double) 1) / ((double) 3));

		cout << "unit cell length is " << a << '\n';

		cout << "Number of ions " << NumberOfIons << '\n';

		int IonNumber = 0;
		for(int i = -1*((int) ceil(R/a)); i <= ((int) ceil(R/a));i++)
			for(int j = -1*((int) ceil(R/a)); j <= ((int) ceil(R/a));j++)
				for(int k = -1*((int) ceil(L/a/2)); k <= ((int) ceil(L/a/2));k++)
				{
					// placing ion in unit cell "origo" corner and checking if all ions have a position
					if((IonNumber < NumberOfIons) && (pow(a*((double) i)/R,2) + pow(a*((double) j)/R,2) + pow(a*((double) k)*2/L,2) <= 1)) //if ion is inside ellipsoid
					{
						ions[IonNumber].SetPosition(0, 0, ((double) i)*a);
						ions[IonNumber].SetPosition(1, 0, ((double) j)*a);
						ions[IonNumber].SetPosition(2, 0, ((double) k)*a);
						IonNumber++;
					}


					// placing ion in unit cell center and checking if all ions have a position
					if((IonNumber < NumberOfIons) && (pow(a*(((double) i) + 0.5)/R,2) + pow(a*(((double) j) + 0.5)/R,2) + pow(a*(((double) k) + 0.5)*2/L,2) <= 1)) //if ion is inside ellipsoid
					{
						ions[IonNumber].SetPosition(0, 0, (((double) i) + 0.5)*a);
						ions[IonNumber].SetPosition(1, 0, (((double) j) + 0.5)*a);
						ions[IonNumber].SetPosition(2, 0, (((double) k) + 0.5)*a);
						IonNumber++;
					}

				}

		if(IonNumber < NumberOfIons)
		{
			cout << NumberOfIons - IonNumber << " ion(s) have not been positioned\n";
			L=L*1.005; // increasing crystal length with 0.5 percent
		}
		else
		{
			cout << "All ions have been positioned\n";
			CrystalNotDone=false;
		}

	}

	// setting all ions to zero-velocity
	for (int i = 0; i < NumberOfIons; i++)
			for (int dim = 0; dim < 3; dim++)
				ions[i].SetVelocity(dim, 0, 0);

}

double Ensemble::Ekin(int TimeStep)
{
	double Ekin = 0;
	for(int n = 0; n < NumberOfIons; n++)
		Ekin += ions[n].Ekin(TimeStep);

	return Ekin;
}

double Ensemble::Ttot(int TimeStep)
{
	return Ekin(TimeStep)/(1.5*NumberOfIons*Kb);
}

double * Ensemble::Ttot(int Tstart, int Tend)
{
	if(Tstart < Tend)
	{
		double *Temp = new double[Tend-Tstart];

		for(int T = Tstart; T < Tend; T++)
			Temp[T] = Ttot(T);

		return Temp;
	}
	else
		return NULL; // invalid time interval
}

double * Ensemble::Tsec(int Tstart, int Tend) // under construction
{
	if(Tstart < Tend)
	{
		double *Temp = new double[Tend-Tstart];
		for(int T = Tstart; T < Tend; T++)
		{
			Temp[T] = 0;
			for(int N = 0; N < NumberOfIons; N++)
				Temp[T] +=42;
		}

		return Temp;
	}
	else
		return NULL;
}

double Ensemble::Tsec(int TimeStep)
{
	if (TimeStep < StepsPrPeriode) // in the beginning there was nothing
	{
		return 0;
	}
	else
	{
		double Ekin = 0;
		for(int N = 0; N < NumberOfIons; N++) // summing up kinetic energy for all ions
		{
			for(int dim = 0; dim < 3; dim++)  // sumining up kinetic energy for x y z
			{
				double vmean = 0; // calculating mean velocity
				for (int t = TimeStep - StepsPrPeriode; t < TimeStep; t++)
					vmean+=Velocity( dim, N, t);
				Ekin+=0.5*Mass(N)*pow(vmean/StepsPrPeriode,2);
			}
		}

		return 2*Ekin/Kb/NumberOfIons/3;
	}
}

double Ensemble::Position(int dim, int N, int TimeStep)
{
	return ions[N].Position(dim, TimeStep);
}


double Ensemble::Velocity(int dim, int N, int TimeStep)
{
	return ions[N].Velocity(dim, TimeStep);
}


void Ensemble::SavePositionToFile()
{
	FILE *f;

	char FilNavn[14];
	_snprintf(FilNavn, 14, "%s%d%s", "MDPos",NumberOfIons, ".xyz"); //Change to _snprintf because windows. 

	//f=fopen("MDpos.xyz","w");
	f=fopen(FilNavn, "w");

	fprintf(f, "%d\n%s\n", NumberOfIons, "Ca");
	for(int N=1; N <= NumberOfIons; N++)
	{
		fprintf(f, "%s%d\t%f\t%f\t%f\n ", "Ca", N, Position(0, N-1, 0)*1e6, Position(1, N-1, 0)*1e6, Position(2, N-1, 0)*1e6);
	}
	fclose(f);
}
// Operators

void Ensemble::SavePositionToFile(int Tend)
{
	FILE *f;

	f=fopen("MDpos.xyz","w");

	fprintf(f, "%d\n%s\n", NumberOfIons, "Ca");
	for(int Tstep = 0; Tstep < Tend; Tstep++)
		for(int N=1; N <= NumberOfIons; N++)
		{
			fprintf(f, "%s%d\t%f\t%f\t%f\n ", "Ca", N, Position(0, N-1, Tstep)*1e6, Position(1, N-1, Tstep)*1e6, Position(2, N-1, Tstep)*1e6);
		}
	fclose(f);

}


double asinh(double x) // making invers sine hyperbolic function (not in cmath lib)
{
	return log(x + sqrt(x*x+1));
}

struct my_f_params // gsl solver need struct with function parameter
{
	double wsquar;
};

double AspectRatioEquation(double alpha, void *p) // from peter herskind thesis
{
	struct my_f_params * params = (struct my_f_params *)p; // pointer magic :-)
	double wsquar = params->wsquar; // aspect ratio

	if (alpha < 1)
		return  -2*(asinh(sqrt(fabs(pow(alpha,-2) - 1))) - alpha*sqrt(fabs(pow(alpha,-2) - 1)))/(asinh(sqrt(fabs(pow(alpha,-2) - 1))) - sqrt(fabs(pow(alpha,-2) - 1))/alpha) - wsquar;
	else
		return  -2*(asin(sqrt(fabs(pow(alpha,-2) - 1))) - alpha*sqrt(fabs(pow(alpha,-2) - 1)))/(asin(sqrt(fabs(pow(alpha,-2) - 1))) - sqrt(fabs(pow(alpha,-2) - 1))/alpha) - wsquar;

}

double CalculateAspectRatio(double wratiosquare)
{
	struct my_f_params params = {wratiosquare};
	gsl_function F;
	F.function = &AspectRatioEquation;
	F.params = &params;

	int status;
	int iter = 0, max_iter = 1000; // maximum number of iterations
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 0;
	double x_lo = 0.0001, x_hi = 100;	// upper and lower bounds
	T = gsl_root_fsolver_bisection; // root finding algorithm
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &F, x_lo, x_hi);

	//printf ("using %s method\n",
		//	gsl_root_fsolver_name (s));


	do
	{
		iter++;
		status = gsl_root_fsolver_iterate (s);
		r = gsl_root_fsolver_root (s);
		x_lo = gsl_root_fsolver_x_lower (s);
		x_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);

		//if (status == GSL_SUCCESS)
		//{
		//	printf ("Converged:\n");
		//}

	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (s); // freeing memory

	//cout << "root is " << r << '\n';

	return r;
}

