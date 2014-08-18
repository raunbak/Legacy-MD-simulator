#include "stdafx.h" // Pre-compiled header.
#include "FastEnsemble.h"
#include "FastIon.h"
#include "constants.h"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <stdarg.h>


using namespace std;

// Constructor
FastEnsemble::FastEnsemble(int m1, int n1, int m2, int n2)
{
	SteadyStateTemperature = SteadyStateTzSec;
	NumberOfIons=n1+n2;
	ions = new FastIon [NumberOfIons];

	for(int n = 0; n < NumberOfIons; n++) {
		if (n<n1 ){
		
  		ions[n].initialize(m1); // setting mass, secular velocity and allocating memory for the ions pos and vel arrays
		}
		else {
		
  		ions[n].initialize(m2); // setting mass, secular velocity and allocating memory for the ions pos and vel arrays
		}
	}


	Vz_IonArray = new double[NumberOfIons];
	Vx_IonArray = new double[NumberOfIons];
	Vy_IonArray = new double[NumberOfIons];

	NumberOfRings = 30;
	VzSecRMS_in_ring = new double[NumberOfRings];
	VrSecRMS_in_ring = new double[NumberOfRings];
	NumberOfIonsInRadius = new int[NumberOfRings];
	RadiiOfIons = new int [NumberOfIons];

	// Setting them all to zero at first.
	for (int i = 0; i < NumberOfRings ; i ++ ) {
		VzSecRMS_in_ring[i]			= 0.0;
		VrSecRMS_in_ring[i]			= 0.0;
		NumberOfIonsInRadius[i]		= 0;
	}
	for(int n = 0; n < NumberOfIons; n++)
	{
		RadiiOfIons[n] = 1;
		Vz_IonArray[n] = 0.0;
		Vx_IonArray[n] = 0.0;
		Vy_IonArray[n] = 0.0;
	}

	ReachedTempArea = false;
}
	

// Member functions
double FastEnsemble::GetTrms(int TimeStep)
{
	return Trms[TimeStep];
}

double FastEnsemble::GetTrmsz(int TimeStep)
{
	return Trmsz[TimeStep];
}

double FastEnsemble::GetTsecular(int TimeStep)
{
	return Tsec[TimeStep];
}

double FastEnsemble::GetTsecularrad(int TimeStep)
{
	return Tsecrad[TimeStep];
}

double FastEnsemble::GetTsecularz(int TimeStep)
{
	return Tsecz[TimeStep];
}

void FastEnsemble::InitialiseTemperatureArrays(int TimeSteps)
{
	Tsec = new double [TimeSteps];
	Tsecrad = new double [TimeSteps];
	Tsecz = new double [TimeSteps];
	Trms = new double [TimeSteps];
	Trmsz = new double [TimeSteps];

	for(int t = 0; t < TimeSteps;t++)
	{
		Tsec[t] = 0;
		Tsecrad[t] = 0;
		Tsecz[t] = 0;
		Trms[t] = 0;
		Trmsz[t] = 0;
	}
}

void FastEnsemble::FreeTemperatureArrays()
{
	delete Tsec;
	delete Trms;
	delete Trmsz;
}

void FastEnsemble::CleanUpEnsemble()
{
	for(int n = 0; n  < NumberOfIons; n++)
		ions[n].CleanUpIon();


	for(int i = 0; i < HistNx;i++)
	{
		for(int j = 0; j < HistNy;j++)
		{
			delete [] histogram[i][j];
			delete [] VelHistogram[i][j];
		}

		delete [] histogram[i];
		delete [] VelHistogram[i];
	}
	delete histogram;
	delete VelHistogram;

}

void FastEnsemble::PrintSecVel()
{
	cout << "secular velocities\n";

	for(int n = 0; n < NumberOfIons; n++)
		cout << ions[n].ReturnVzSec()/StepsPrPeriode << ";\n";
}

void FastEnsemble::SetActualTemperature(double val)
{
	ActualAvgTemperature = val;
}

void FastEnsemble::SetActualTemperatureSTD(double val)
{
	ActualTemperatureSTD = val;
}

double FastEnsemble::GetActualTemperature()
{
	return ActualAvgTemperature;
}

double FastEnsemble::GetActualTemperatureSTD()
{
	return ActualTemperatureSTD;
}

double FastEnsemble::GetCurrentTemperature()
{
	return (pow(VzSecRMS,2)+pow(VrSecRMS,2))*Mass(0)/Kb/3; // REMEBER TO CHANGE THIS FOR MASSES
}

void FastEnsemble::RescaleVelocityDistributionIndivialIons()
{
	double* a_x = new double[NumberOfIons];
	double* a_y = new double[NumberOfIons];
	double* a_z = new double[NumberOfIons];

	for (int i = 0; i < NumberOfIons ; i ++ ) 
	{


		double CurrentTempx	= pow(Vx_IonArray[i],2)*ions[i].GetMass()/(Kb);
		double CurrentTempy	= pow(Vy_IonArray[i],2)*ions[i].GetMass()/(Kb);
		double CurrentTempz	= pow(Vz_IonArray[i],2)*ions[i].GetMass()/(Kb);

		a_x[i]						= sqrt((SteadyStateTemperature/3)/CurrentTempx);
		a_y[i]						= sqrt((SteadyStateTemperature/3)/CurrentTempy);
		a_z[i]						= sqrt((SteadyStateTemperature/3)/CurrentTempz);

		// This is so goddamn ugly it makes my eyes hurt.
		if (a_x[i] > 1.05)
		{
			a_x[i] = 1.05;
		} 
		if (a_x[i] < 0.95)
		{
			a_x[i] = 0.95;
		}

		if (a_y[i] > 1.05)
		{
			a_y[i] = 1.05;
		} 
		if (a_y[i] < 0.95)
		{
			a_y[i] = 0.95;
		}

		if (a_z[i] > 1.05)
		{
			a_z[i] = 1.05;
		} 
		if (a_z[i] < 0.95)
		{
			a_z[i] = 0.95;
		}


		double VelX = ions[i].Velocity(0)*a_x[i];
		double VelY = ions[i].Velocity(1)*a_y[i];
		double VelZ = ions[i].Velocity(2)*a_z[i];

		
		ions[i].SetVelocity(0, VelX);
		ions[i].SetVelocity(1, VelY);
		ions[i].SetVelocity(2, VelZ);

	}
	
}

void FastEnsemble::RescaleVelocityDistributionRadial()
{

	double* CurrentTemp			= new double[NumberOfRings];
	double* CurrentRadialTemp	= new double[NumberOfRings];
	double* a					= new double[NumberOfRings];	
	double* ar					= new double[NumberOfRings];


	for (int i = 0; i < NumberOfRings ; i ++ ) 
	{
		
		if ( NumberOfIonsInRadius[i] != 0) 
		{

		//double Tz				= pow(VzSecRMS_in_ring[i],2)*Mass(0)/Kb;
		//double Tr				= pow(VrSecRMS_in_ring[i],2)*Mass(0)/Kb/2;

		double Tz				= VzSecRMS_in_ring[i]*Mass(0)/Kb;
		double Tr				= VrSecRMS_in_ring[i]*Mass(0)/Kb/2;
		
		CurrentTemp[i]			= Tz;
		CurrentRadialTemp[i]	= Tr;

		a[i]					= sqrt(SteadyStateTemperature/CurrentTemp[i]);
		ar[i]					= sqrt(SteadyStateTemperature/CurrentRadialTemp[i]);

		if (a[i] > 1.05)
		{
			a[i] = 1.05;
		} 
		if (a[i] < 0.95)
		{
			a[i] = 0.95;
		}
		if (ar[i] > 1.05)
		{
			ar[i] = 1.05;
		}
		if (ar[i] < 0.95)
		{
			ar[i] = 0.95;
		}
		}
		else
		{
			a[i] = 1.00;
			ar[i] = 1.00;
		}

	}

	// And now we rescale the velocity.
	double VelX = 0.0;
	double VelY = 0.0;
	double VelZ = 0.0;

	for(int n = 0; n < NumberOfIons; n++)
	{

		double arNow = (double) ar[RadiiOfIons[n]-1];
		double aNow = (double) a[RadiiOfIons[n]-1];

		VelX = ions[n].Velocity(0)*arNow;
		VelY = ions[n].Velocity(1)*arNow;
		VelZ = ions[n].Velocity(2)*aNow;

		ions[n].SetVelocity(0, VelX);
		ions[n].SetVelocity(1, VelY);
		ions[n].SetVelocity(2, VelZ);
	}

}

void FastEnsemble::RescaleVelocityXYZ(double Total_V_x_rms,double Total_V_y_rms,double Total_V_z_rms)
{
		// calculate current temperature
	double CurrentTempZ = (pow(Total_V_z_rms,2)*Mass(0))/(Kb);
	double CurrentTempX = (pow(Total_V_x_rms,2)*Mass(0))/(Kb);
	double CurrentTempY = (pow(Total_V_y_rms,2)*Mass(0))/(Kb);

	// rescale velocity distribution (Maybe change which fraction of SST is given)
	 
	cout << sqrt(pow(CurrentTempZ,2)+pow(CurrentTempX,2) + pow(CurrentTempY,2) )<< endl;
	double T_fraction = SteadyStateTemperature;

	double az = sqrt((T_fraction)/CurrentTempZ);
	double ax = sqrt((T_fraction)/CurrentTempX);
	double ay = sqrt((T_fraction)/CurrentTempY);

	cout << ax <<" " << ay << " " << az << endl; 

	double lowerLimit = 0.98;//sqrt(100/102.5); // Hvornår der ikke skal køles mere.
	double upperLimit = 1.02;//sqrt(100/97.5); // Hvornår der ikke skal varmes op mere.

	//cout << lowerLimit << " " << upperLimit << endl;
	// Good value for raise temp is 1.005

		if (ax > upperLimit)
		{
			ax =1.0000002;
		} 
		else if (ax < lowerLimit)
		{
			ax = 0.9990;
		}
		else {
			ax = 1.0;
		}

		if (ay > upperLimit)
		{
			ay =1.0000002;
		}
		else if (ay < lowerLimit)
		{
			ay = 0.9990;
		}
		else {

			ay = 1.0;
		}

		if (az > upperLimit)
		{
			az = 1.0000002;
		}
		else if (az < lowerLimit)
		{
			az = 0.9990;
		}
		else {
			az = 1.0;
		}
		
		/*
		if (ax < lowerLimit)
		{
			ax = 0.9990;
		}
		else {

			if ( ReachedTempArea == true)
			{
			ax = 0.9999;
			ReachedTempArea = false;
			}
			else
			{
			ax = 1.0001;
			ReachedTempArea = true;
			}
		}

		if (ay < lowerLimit)
		{
			ay = 0.9990;
		}
		else {

			if ( ReachedTempArea == true)
			{
			ay = 0.9999;
			ReachedTempArea = false;
			}
			else
			{
			ay = 1.0001;
			ReachedTempArea = true;
			}
		}

		if (az < lowerLimit)
		{
			az = 0.9990;
		}
		else {
			if ( ReachedTempArea == true)
			{
			az = 0.99999;
			ReachedTempArea = false;
			}
			else
			{
			az = 1.00001;
			ReachedTempArea = true;
			}
		}
		
		if (ax < lowerLimit)
		{
			ax = 0.9990;
		}
		else {

			ax = 1.0000;
		}

		if (ay < lowerLimit)
		{
			ay = 0.9990;
		}
		else {

	       ay = 1.0000;
		}

		if (az < lowerLimit)
		{
			az = 0.9990;
		}
		else {
			az = 1.0000;
		}
		*/
		cout <<"ax = " << ax << "  ay = " << ay << " az= " << az  << endl;

	for(int n = 0; n < NumberOfIons; n++)
	{
		ions[n].SetVelocity(0, ions[n].Velocity(0)*ax);
		ions[n].SetVelocity(1, ions[n].Velocity(1)*ay);
		ions[n].SetVelocity(2, ions[n].Velocity(2)*az);
	}
}

void FastEnsemble::RescaleVelocityDumb(double VrmsZ, double VrmsR)
{
		// calculate current temperature
	double CurrentTemp = (pow(VrmsZ,2)*Mass(0))/Kb;
	double CurrentRadialTemp = (pow(VrmsR,2)*Mass(0))/2*Kb;

	//cout << CurrentTemp << " V   Z" << endl;
	//cout << CurrentRadialTemp << " V   R" << endl;

	// rescale velocity distribution (Maybe change which fraction of SST is given)
	double a = sqrt(SteadyStateTemperature/CurrentTemp);
	double ar = sqrt(SteadyStateTemperature/CurrentRadialTemp);

	if (a > 1.05)
	{
		a = 1.05;
	}

	if (a < 0.95)
		a = 0.95;

	if (ar > 1.05)
	{
		ar = 1.05;
	}

	if (ar < 0.95)
		ar = 0.95;

	for(int n = 0; n < NumberOfIons; n++)
	{
		ions[n].SetVelocity(0, ions[n].Velocity(0)*ar);
		ions[n].SetVelocity(1, ions[n].Velocity(1)*ar);
		ions[n].SetVelocity(2, ions[n].Velocity(2)*a);
	}
}

void FastEnsemble::RescaleVelocityDistribution()
{
	// calculate current temperature
	//double CurrentTemp = pow(VzSecRMS,2)*Mass(0)/Kb;
	//double CurrentRadialTemp = pow(VrSecRMS,2)*Mass(0)/Kb/2;
	double CurrentTemp = pow(VzSecRMS,2)*Mass(0)/Kb;
	double CurrentRadialTemp = pow(VrSecRMS,2)*Mass(0)/Kb/2;
	//cout << "Current Secular z-Temperature is " << CurrentTemp*1000 << " mK \n" << endl;
	//cout << CurrentTemp*1000 << '\n' << endl;

	// rescale velocity distribution
	double a = sqrt((SteadyStateTemperature/3)/CurrentTemp);
	double ar = sqrt((SteadyStateTemperature * (2/3))/CurrentRadialTemp);

	
	if (a > 1.05)
	{
		a = 1.05;
	}
	if (a < 0.95)
	{
		a = 0.95;
	}
	if (ar > 1.05)
	{
		ar = 1.05;
	}
	if (ar < 0.95)
	{
		ar = 0.95;
	}

	for(int n = 0; n < NumberOfIons; n++)
	{
		ions[n].SetVelocity(0, ions[n].Velocity(0)*ar);
		ions[n].SetVelocity(1, ions[n].Velocity(1)*ar);
		ions[n].SetVelocity(2, ions[n].Velocity(2)*a);
	}

}

void FastEnsemble::RescaleVelocityDistribution(int TimeStep)
{
	// calculate current temperature
	double CurrentTemp = Tsecz[TimeStep];
	double CurrentRadialTemp = Tsecrad[TimeStep];
	//double CurrentRadialTemp = Trms[TimeStep]; // temperature based on avg velocity speed
	//cout << "Current Secular z-Temperature is " << CurrentTemp*1000 << " mK \n";
	//cout << CurrentTemp*1000 << '\n';

	// rescale velocity distribution
	double a = sqrt(SteadyStateTemperature/CurrentTemp);
	double ar = sqrt(SteadyStateTemperature/CurrentRadialTemp);

	if (a > 1.05)
	{
		a = 1.05;
		//cout << "To big temperature jump rescaling factor is set to " << a << '\n';
	}

	if (a < 0.95)
		a = 0.95;

	if (ar > 1.05)
	{
		ar = 1.05;
		//cout << "To big temperature jump rescaling factor is set to " << a << '\n';
	}

	if (ar < 0.95)
		ar = 0.95;

	for(int n = 0; n < NumberOfIons; n++)
	{
		ions[n].SetVelocity(0, ions[n].Velocity(0)*ar);
		ions[n].SetVelocity(1, ions[n].Velocity(1)*ar);
		ions[n].SetVelocity(2, ions[n].Velocity(2)*a);
	}
}

void FastEnsemble::VelocityKick(int TimeStep)
{
	// calculate random direction for each ion
	double theta = ((double) rand())/((double) RAND_MAX) * 2 * PI;
	double phi = ((double) rand())/((double) RAND_MAX) * PI;

	for(int n = 0; n < NumberOfIons; n++)
	{
		ions[n].SetVelocity(0, ions[n].Velocity(0) + cos(theta)*sin(phi)*vkick);
		ions[n].SetVelocity(1, ions[n].Velocity(1) + sin(theta)*sin(phi)*vkick);
		ions[n].SetVelocity(2, ions[n].Velocity(2) + cos(phi)*vkick);
	}

}

void FastEnsemble::InitialiseVelocityHistogram()
{
	//allocating memory to histogram
	VelHistogram = new long double ** [HistNx];
	for(int i = 0; i < HistNx;i++)
	{
		VelHistogram[i] = new long double *[HistNy];
		for(int j = 0; j < HistNy;j++)
			VelHistogram[i][j] = new long double [HistNz];
	}
	for(int i=0; i < HistNx; i++)
		for(int j=0; j < HistNy; j++)
			for(int k=0; k < HistNz; k++)
				VelHistogram[i][j][k] = 0;
}

void FastEnsemble::InitialiseCountHistogram()
{
	//allocating memory to histogram
	CountHistogram = new long int ** [HistNx];
	for(int i = 0; i < HistNx;i++)
	{
		CountHistogram[i] = new long int *[HistNy];
		for(int j = 0; j < HistNy;j++)
			CountHistogram[i][j] = new long int [HistNz];
	}
	for(int i=0; i < HistNx; i++)
		for(int j=0; j < HistNy; j++)
			for(int k=0; k < HistNz; k++)
				CountHistogram[i][j][k] = 0;
}

void FastEnsemble::InitialiseHistogram()
{
	//allocating memory to histogram
	histogram = new long int ** [HistNx];
	for(int i = 0; i < HistNx;i++)
	{
		histogram[i] = new long int *[HistNy];
		for(int j = 0; j < HistNy;j++)
			histogram[i][j] = new long int [HistNz];
	}
	for(int i=0; i < HistNx; i++)
		for(int j=0; j < HistNy; j++)
			for(int k=0; k < HistNz; k++)
				histogram[i][j][k] = 0;

}

void FastEnsemble::UpdateHistogram()
{
	for (int i = 0; i < NumberOfIons; i++)
	{
		int Nx = ((int) ((ions[i].Position(0))/PixelToDistance+((double) HistNx)/2));
		int Ny = ((int) ((ions[i].Position(1))/PixelToDistance+((double) HistNy)/2));
		int Nz = ((int) ((ions[i].Position(2))/PixelToDistance+((double) HistNz)/2));

		if (Nx < HistNx && Ny < HistNy && Nz < HistNz && Nx > 0 && Ny > 0 && Nz > 0)
		{
			histogram[Nx][Ny][Nz]++;
		}
	}
}

void FastEnsemble::UpdateVelocityHistogram()
{
	for (int i = 0; i < NumberOfIons; i++)
	{
		int Nx = ((int) ((ions[i].Position(0))/PixelToDistance+((double) HistNx)/2));
		int Ny = ((int) ((ions[i].Position(1))/PixelToDistance+((double) HistNy)/2));
		int Nz = ((int) ((ions[i].Position(2))/PixelToDistance+((double) HistNz)/2));

		if (Nx < HistNx && Ny < HistNy && Nz < HistNz && Nx > 0 && Ny > 0 && Nz > 0)
		{
			VelHistogram[Nx][Ny][Nz]+= ((long double) pow(ions[i].GetVsec(),2));
		}
	}
}

void FastEnsemble::MyUpdateVelocityHistogram()
{
	for (int i = 0; i < NumberOfIons; i++)
	{
		int Nx = ((int) ((ions[i].Position(0))/PixelToDistance+((double) HistNx)/2));
		int Ny = ((int) ((ions[i].Position(1))/PixelToDistance+((double) HistNy)/2));
		int Nz = ((int) ((ions[i].Position(2))/PixelToDistance+((double) HistNz)/2));

		if (Nx < HistNx && Ny < HistNy && Nz < HistNz && Nx > 0 && Ny > 0 && Nz > 0)
		{
			VelHistogram[Nx][Ny][Nz]+= ((long double) pow(ions[i].Velocity(),2));
		}
	}
}

void FastEnsemble::UpdateCountHistogram()
{
	for (int i = 0; i < NumberOfIons; i++)
	{
		int Nx = ((int) ((ions[i].Position(0))/PixelToDistance+((double) HistNx)/2));
		int Ny = ((int) ((ions[i].Position(1))/PixelToDistance+((double) HistNy)/2));
		int Nz = ((int) ((ions[i].Position(2))/PixelToDistance+((double) HistNz)/2));

		if (Nx < HistNx && Ny < HistNy && Nz < HistNz && Nx > 0 && Ny > 0 && Nz > 0)
		{
			CountHistogram[Nx][Ny][Nz]++;
		}
	}
}

double FastEnsemble::ReturnHist(int i, int j, int k)
{
	return histogram[i][j][k];
}

double FastEnsemble::ReturnVelHist(int i, int j, int k)
{
	return ((double) VelHistogram[i][j][k]);
}

double FastEnsemble::ReturnCountHist(int i, int j, int k)
{
	return ((double) CountHistogram[i][j][k]);
}

void FastEnsemble::SetSteadyStateTemperature(double Val)
{
	SteadyStateTemperature = Val;
}

int FastEnsemble::GetNumberOfIons()
{
	return NumberOfIons;
}

double FastEnsemble::Mass(int N)
{
	return ions[N].GetMass();
}

void FastEnsemble::CrystalGenerator()
{
	// Maybe change the initial velocities.
	for(int N = 0; N < NumberOfIons ; N++)
	{
		// setting initial positions
		ions[N].SetPosition(0, pow((double) -1,(double) N)*1e-7);
		ions[N].SetPosition(1, pow((double) -1,(double) N+1)*1e-7);
		ions[N].SetPosition(2, GridSpacing*(N+1) - 0.5*GridSpacing*(NumberOfIons-1)); // making a string of ions

		// setting initial velocities
		double Vmean = 0.0;
		if (Tinitial != 0.0) 
		{
			Vmean = sqrt(Kb*Tinitial/ions[N].GetMass());
		}
		ions[N].SetVelocity(0, 0); //(rand()/RAND_MAX-0.5)*Vmean
		ions[N].SetVelocity(1, Vmean);
		ions[N].SetVelocity(2, Vmean);
	}
}

double FastEnsemble::getRho0(double Vrf){

	return eps0*Vrf*Vrf/(Mass(0)*pow(r0,4)*OmegaRF*OmegaRF);
}

void FastEnsemble::CrystalGenerator(double Vrf, double Vend)
{
	// Calculating pseudo trap frequencies
	double wz2=2*eta*e*Vend/pow(z0,2)/Mass(0); // ERROR in formula from Magnus' thesis! he means potential and not electric potential.
	double wr2=pow(e*Vrf/Mass(0)/OmegaRF,2)/2/pow(r0,4)-eta*e*Vend/Mass(0)/pow(z0,2);

	// Solving plasma model aspect-ratio equation
	double alpha = CalculateAspectRatio(wz2/wr2);
	//double alpha = CalculateAspectRatio(2);
	//cout << "aspect ratio is " << alpha << '\n';
	/*cout << "out putting ratio(alpha)\n";
	for (double ratio = 0.01; ratio < 2; ratio+=0.01)
		if (ratio != 1)
			cout << CalculateAspectRatio(ratio*ratio) << ';' << '\n';
	 */
	// Calculating number density
	double rho0 = (eps0*Vrf*Vrf)/(Mass(0)*pow(r0,4)*OmegaRF*OmegaRF);

	// Calculating volume of crystal
	double V = ((double) NumberOfIons)/rho0;
	//cout << "volume of crystal is " << V << '\n';

	// Calculating length and radius of crystal
	double L = pow(6*V/(PI*alpha*alpha), ((double) 1) /  ((double) 3));
	bool CrystalNotDone = true;
	while(CrystalNotDone)
	{
		double R = alpha*L/2;

		//cout << "radius is " << R << '\n';
		//cout << "Length is " << L << '\n';
		//cout << "2R/L=" << 2*R/L << '\n';

		// unit cell cube length for bcc structure
		double a = pow( 2 / rho0, ((double) 1) / ((double) 3));

		//cout << "unit cell length is " << a << '\n';

		//cout << "Number of ions " << NumberOfIons << '\n';

		int IonNumber = 0;
		for(int i = -1*((int) ceil(R/a)); i <= ((int) ceil(R/a));i++)
			for(int j = -1*((int) ceil(R/a)); j <= ((int) ceil(R/a));j++)
				for(int k = -1*((int) ceil(L/a/2)); k <= ((int) ceil(L/a/2));k++)
				{
					// placing ion in unit cell "origo" corner and checking if all ions have a position
					if((IonNumber < NumberOfIons) && (pow(a*((double) i)/R,2) + pow(a*((double) j)/R,2) + pow(a*((double) k)*2/L,2) <= 1)) //if ion is inside ellipsoid
					{
						ions[IonNumber].SetPosition(0, ((double) i)*a);
						ions[IonNumber].SetPosition(1, ((double) j)*a);
						ions[IonNumber].SetPosition(2, ((double) k)*a);
						IonNumber++;
					}


					// placing ion in unit cell center and checking if all ions have a position
					if((IonNumber < NumberOfIons) && (pow(a*(((double) i) + 0.5)/R,2) + pow(a*(((double) j) + 0.5)/R,2) + pow(a*(((double) k) + 0.5)*2/L,2) <= 1)) //if ion is inside ellipsoid
					{
						ions[IonNumber].SetPosition(0, (((double) i) + 0.5)*a);
						ions[IonNumber].SetPosition(1, (((double) j) + 0.5)*a);
						ions[IonNumber].SetPosition(2, (((double) k) + 0.5)*a);
						IonNumber++;
					}

				}

		if(IonNumber < NumberOfIons)
		{
			//cout << NumberOfIons - IonNumber << " ion(s) have not been positioned\n";
			L=L*1.005; // increasing crystal length with 0.5 percent
		}
		else
		{
			Radius = R;
			Length =L;
			//cout << "All ions have been positioned\n";
			CrystalNotDone=false;
		}


		//cout << a << " a_WS" <<endl;
	}

	// NORMALT ER DEN NUL! LIGE NU SVARER DEN TIL 1K !
	// setting all ions to zero-velocity
	for (int i = 0; i < NumberOfIons; i++)
	{
		for (int dim = 0; dim < 3; dim++)
		{
			double vel = rand() % 10;
			if(vel < 5 ) {
				vel = -1.0;
			}
			else {
				vel = 1.0;
			}
			ions[i].SetVelocity(dim, vel*5.5838);
		}

	}
	cout << Radius << " = Radius of crystal " << Length << " = length of crystal"<< endl; //REMOVE THIS!!
}

void FastEnsemble::CrystalGeneratorForPseudoPotential(double Vrf, double Vend)
{
	// Calculating pseudo trap frequencies
	double wz2=2*eta*e*Vend/pow(z0,2)/Mass(0); // ERROR in formula from Magnus' thesis! he means potential and not electric potential.
	double wr2=pow(e*Vrf/Mass(0)/OmegaRF,2)/2/pow(r0,4)-eta*e*Vend/Mass(0)/pow(z0,2);

	// Solving plasma model aspect-ratio equation
	double alpha = CalculateAspectRatio(wz2/wr2);
	//double alpha = CalculateAspectRatio(2);
	//cout << "aspect ratio is " << alpha << '\n';
	/*cout << "out putting ratio(alpha)\n";
	for (double ratio = 0.01; ratio < 2; ratio+=0.01)
		if (ratio != 1)
			cout << CalculateAspectRatio(ratio*ratio) << ';' << '\n';
	 */
	// Calculating number density
	double rho0 = eps0*Vrf*Vrf/(Mass(0)*pow(r0,4)*OmegaRF*OmegaRF);

	cout << "this is rho0 :" << rho0 << endl;

	// Calculating volume of crystal
	double V = ((double) NumberOfIons)/rho0;
	//cout << "volume of crystal is " << V << '\n';

	// Calculating length and radius of crystal
	double L = pow(6*V/(PI*alpha*alpha), ((double) 1) /  ((double) 3));
	bool CrystalNotDone = true;
	while(CrystalNotDone)
	{
		double R = alpha*L/2;

		//cout << "radius is " << R << '\n';
		//cout << "Length is " << L << '\n';
		//cout << "2R/L=" << 2*R/L << '\n';

		// unit cell cube length for bcc structure
		double a = pow( 2 / rho0, ((double) 1) / ((double) 3));

		//cout << "unit cell length is " << a << '\n';

		//cout << "Number of ions " << NumberOfIons << '\n';

		int IonNumber = 0;
		for(int i = -1*((int) ceil(R/a)); i <= ((int) ceil(R/a));i++)
			for(int j = -1*((int) ceil(R/a)); j <= ((int) ceil(R/a));j++)
				for(int k = -1*((int) ceil(L/a/2)); k <= ((int) ceil(L/a/2));k++)
				{
					// placing ion in unit cell "origo" corner and checking if all ions have a position
					if((IonNumber < NumberOfIons) && (pow(a*((double) i)/R,2) + pow(a*((double) j)/R,2) + pow(a*((double) k)*2/L,2) <= 1)) //if ion is inside ellipsoid
					{
						ions[IonNumber].SetPosition(0, ((double) i)*a);
						ions[IonNumber].SetPosition(1, ((double) j)*a);
						ions[IonNumber].SetPosition(2, ((double) k)*a);
						IonNumber++;
					}


					// placing ion in unit cell center and checking if all ions have a position
					if((IonNumber < NumberOfIons) && (pow(a*(((double) i) + 0.5)/R,2) + pow(a*(((double) j) + 0.5)/R,2) + pow(a*(((double) k) + 0.5)*2/L,2) <= 1)) //if ion is inside ellipsoid
					{
						ions[IonNumber].SetPosition(0, (((double) i) + 0.5)*a);
						ions[IonNumber].SetPosition(1, (((double) j) + 0.5)*a);
						ions[IonNumber].SetPosition(2, (((double) k) + 0.5)*a);
						IonNumber++;
					}

				}

		if(IonNumber < NumberOfIons)
		{
			//cout << NumberOfIons - IonNumber << " ion(s) have not been positioned\n";
			L=L*1.005; // increasing crystal length with 0.5 percent
		}
		else
		{

			//cout << "All ions have been positioned\n";
			CrystalNotDone=false;
		}

	}

	

	// setting all ions to steadystate temperature
	for (int i = 0; i < NumberOfIons; i++)
	{

		// setting initial velocities
		double Vmean = sqrt(Kb*SteadyStateTemperature/ions[i].GetMass());
		double RandomAngle=rand()/RAND_MAX*PI;
		ions[i].SetVelocity(0, 0); //(rand()/RAND_MAX-0.5)*Vmean
		ions[i].SetVelocity(1, Vmean*sin(RandomAngle));
		ions[i].SetVelocity(2, Vmean*cos(RandomAngle));
	}

}

double FastEnsemble::Ekin()
{
	double Ekin = 0;
	for(int n = 0; n < NumberOfIons; n++)
		Ekin += ions[n].Ekin();

	return Ekin;
}

double FastEnsemble::Ttot()
{
	return Ekin()/(1.5*NumberOfIons*Kb);
}

double FastEnsemble::Position(int dim, int N)
{
	return ions[N].Position(dim);
}

double FastEnsemble::Velocity(int dim, int N)
{
	return ions[N].Velocity(dim);
}

void FastEnsemble::SavePositionToFile()
{
	FILE *f;

	f=fopen("MDpos.xyz","w");

	fprintf(f, "%d\n%s\n", NumberOfIons, "Ca");
	for(int N=1; N <= NumberOfIons; N++)
		fprintf(f, "%s%d\t%f\t%f\t%f\n ", "Ca", N, Position(0, N-1)*1e6, Position(1, N-1)*1e6, Position(2, N-1)*1e6);

	fclose(f);
}

void FastEnsemble::SaveIonDataToFile()
{
	// Saves ion data to file
	ofstream Ionfile ("IonData.txt");

	Ionfile << "N \t mass \t x \t y \t z \t Vx \t Vy \t Vz \t Vsec \t norm(V)" << endl; 
	Ionfile << "R: "<< Radius <<" L: " << Length << endl; 
	for(int N = 0; N < GetNumberOfIons(); N++)
	{
		if (Ionfile.is_open())
		{
			Ionfile <<			N	+ 1			<< ",  " <<
								ions[N].GetMass()		<< ",  " <<
								ions[N].Position(0)		<< ",  " <<
								ions[N].Position(1)		<< ",  " <<
								ions[N].Position(2)		<< ",  " <<
								ions[N].Velocity(0)		<< ",  " <<
								ions[N].Velocity(1)		<< ",  " <<
								ions[N].Velocity(2)		<< ",  " <<
								ions[N].ReturnVzSec()	<< ",  " <<
								ions[N].Velocity()		<< ",  " <<
								endl;
		}   

	}
	Ionfile.close();
    

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
