// integrator.cpp
#include "stdafx.h" // Pre-compiled header.
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "integrator.h"
#include "forces.h"
#include "constants.h"
#include "cudaforces.cuh"


#include <time.h>

using namespace std;

// simpel leap frog algorithm as in Matthey thesis - written out for readability, NOT OPTIMIZED!

void LeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend)
{
	// allocate memory to the histogram
	ensemble.InitialiseHistogram();
	ensemble.InitialiseVelocityHistogram();

	// making tempeary position and velocities
	double **TempPos;
	double **TempVel;


	TempPos = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempPos[i] = new double [3];

	TempVel = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempVel[i] = new double [3];


	int counter = 1;

	double AvgTemp=0;
	double TempSquare=0;
	int RandomNumber = 1;
	double TotalVzSquared = 0;
	double TotalTemp = 0;
	double TotalTempStd = 0;

	for(int t=1; t < TimeSteps ; t++)
	{
		TotalVzSquared = 0;
		// updating position, looping through ions and x-y-z
		for(int N=0; N < ensemble.NumberOfIons; N++)
		{
			for(int dim=0; dim < 3; dim++)
			{
				TempPos[N][dim]  = ensemble.ions[N].Pos[dim] + dt*dt*(Ftot(ensemble, N, t-1, dim, Vrf, Vend))/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[dim];
				TempVel[N][dim]  = ensemble.ions[N].Vel[dim] + dt/ensemble.ions[N].m*Ftot(ensemble, N, t-1, dim, Vrf, Vend);

				// for the velocity histogram
				ensemble.ions[N].Vavg[dim] += TempVel[N][dim];
			}
			ensemble.ions[N].VzSec = ensemble.ions[N].VzSec + TempVel[N][2];
			TotalVzSquared += pow(TempVel[N][2],2);

		}

		TotalTemp += ensemble.ions[0].m/Kb*TotalVzSquared/ensemble.NumberOfIons;
		TotalTempStd += pow(ensemble.ions[0].m/Kb*TotalVzSquared/ensemble.NumberOfIons,2);

		if (((t-1) % ((int) StepsPrPeriode)) == 0)
			for(int N=0; N < ensemble.NumberOfIons; N++)
			{
				// calculating secular velocity
				for(int dim=0; dim < 3; dim++)
					ensemble.ions[N].Vavg[dim] = ensemble.ions[N].Vavg[dim]/StepsPrPeriode;
				// setting ions secular velocity
				ensemble.ions[N].Vsec = sqrt(pow(ensemble.ions[N].Vavg[0],2) + pow(ensemble.ions[N].Vavg[1],2) + pow(ensemble.ions[N].Vavg[2],2));

				// setting to zero
				for(int dim=0; dim < 3; dim++)
					ensemble.ions[N].Vavg[dim] = 0;
			}


		// saving time step
		for(int N=0; N < ensemble.NumberOfIons; N++)
			for(int dim=0; dim < 3; dim++)
			{
				ensemble.ions[N].Pos[dim] = TempPos[N][dim];
				ensemble.ions[N].Vel[dim] = TempVel[N][dim];
			}

		if (t > StartRecordingHistogram)
		{
			ensemble.UpdateHistogram();
			ensemble.UpdateVelocityHistogram();
		}

		/*
		  // setting temperature every rf-cycle - not good for big crystals
		  if (((t-1) % ((int) StepsPrPeriode/2)) == 0) // Dealing with velocity distribution
		  {
			  if (t > 1 && counter == 3)
			  {
				  counter = 0;
				  ensemble.VzSecRMS = 0;
				  for(int N=0; N < ensemble.NumberOfIons; N++)
				  {
					  ensemble.VzSecRMS = ensemble.VzSecRMS + pow(ensemble.ions[N].VzSec,2);
					  ensemble.ions[N].VzSec = 0;
				  }
				  ensemble.VzSecRMS = sqrt(ensemble.VzSecRMS/(ensemble.NumberOfIons*pow(StepsPrPeriode,2)));
				  // rescale velocity distribution
				  ensemble.RescaleVelocityDistribution();
			  }
			  else if (counter == 1)
				  for(int N=0; N < ensemble.NumberOfIons; N++)
					  ensemble.ions[N].VzSec = 0;
			  counter++;
		  }*/

		if (((t-1) % ((int) StepsPrPeriode)) == 0) // Dealing with velocity distribution
		{

			AvgTemp += ensemble.GetCurrentTemperature()*1000/((double) TimeSteps)*StepsPrPeriode ;
			TempSquare += pow(ensemble.GetCurrentTemperature()*1000,2)/((double) TimeSteps)*StepsPrPeriode;
			ensemble.VzSecRMS = 0;
			for(int N=0; N < ensemble.NumberOfIons; N++)
			{
				ensemble.VzSecRMS = ensemble.VzSecRMS + pow(ensemble.ions[N].VzSec,2);
				ensemble.ions[N].VzSec = 0;
			}
			ensemble.VzSecRMS = sqrt(ensemble.VzSecRMS/(ensemble.NumberOfIons*pow(StepsPrPeriode,2)));
			// rescale velocity distribution
			//cout << ensemble.GetCurrentTemperature()*1000 << ";\n";
			//if (((t-1) % ((int) (10*StepsPrPeriode))) == 0)
			//ensemble.RescaleVelocityDistribution();



			for(int N=0; N < ensemble.NumberOfIons; N++)
				ensemble.ions[N].VzSec = 0;
		}

		if (((t-RandomNumber) % ((int) (3*StepsPrPeriode))) == 0)
		{
			RandomNumber = rand() % ((int) StepsPrPeriode-1) + 1;
			while ( ((2 < RandomNumber) && (13 > RandomNumber)) || ((17 < RandomNumber) && (28 > RandomNumber)))
				RandomNumber = rand() % ((int) StepsPrPeriode-1) + 1;
			ensemble.RescaleVelocityDistribution();
		}


	}

	ensemble.SetActualTemperature(AvgTemp);
	ensemble.SetActualTemperatureSTD(sqrt(TempSquare-pow(AvgTemp,2)));
	TotalTemp = TotalTemp/((double) TimeSteps)*1000;
	TotalTempStd = sqrt(TotalTempStd/((double) TimeSteps)*1000*1000-pow(TotalTemp,2));
	cout << "Total z Velocity " << TotalTemp << "+/-" << TotalTempStd << "Secular z Velocity " << ensemble.GetActualTemperature() << '\n';
	cout << ensemble.GetActualTemperature() << '\t' << ensemble.GetActualTemperatureSTD() <<'\n';

}

// Normal Cuda integrator
void CudaLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend)
{
	// allocate memory to the histogram
	ensemble.InitialiseHistogram();
	ensemble.InitialiseVelocityHistogram();

	// make Cuda Pos and Force arrays - its float because we dont want to overflow the GPU
	float *PosX, *PosY, *PosZ, *ForceX, *ForceY, *ForceZ;

	PosX = new float [ensemble.GetNumberOfIons()];
	PosY = new float [ensemble.GetNumberOfIons()];
	PosZ = new float [ensemble.GetNumberOfIons()];
	ForceX = new float [ensemble.GetNumberOfIons()];
	ForceY = new float [ensemble.GetNumberOfIons()];
	ForceZ = new float [ensemble.GetNumberOfIons()];

	// declare and allocate memory on GPU
	float *PosX_d, *PosY_d, *PosZ_d, *ForceX_d, *ForceY_d, *ForceZ_d;
	CudaCoulombAlloc(&PosX_d, &PosY_d, &PosZ_d, &ForceX_d, &ForceY_d, &ForceZ_d, ensemble.GetNumberOfIons());


	// making tempeary position and velocities
	double **TempPos;
	double **TempVel;


	TempPos = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempPos[i] = new double [3];

	TempVel = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempVel[i] = new double [3];


	int counter = 1;

	double AvgTemp=0;
	double TempSquare=0;
	int RandomNumber = 1;
	double TotalVzSquared = 0;
	double TotalVrSquared = 0;
	double TotalTemp = 0;
	double TotalTempStd = 0;

	double TotalRadialTemp = 0;
	double TotalVxSquared = 0; 
	double TotalVySquared = 0; 
	double TotalTempx = 0;
	double TotalTempy = 0;
	double TotalV = 0;
	ofstream Temperaturefile ("TemperatureData.txt");
	Temperaturefile << "First line of output" << endl ;
	
	// writing out to textfile
	FILE *f;

	char FilNavn[14];
	_snprintf(FilNavn, 14, "%s%d%s", "MDPos", ensemble.NumberOfIons, ".xyz"); //Change to _snprintf

	f=fopen(FilNavn, "w");
	// writing start format to xyz file
	fprintf(f, "%d\n%s\n", ensemble.NumberOfIons, "Ca");

	//
	for(int t=1; t < TimeSteps ; t++)
	{

		// Prepare Data for CUDA transfer
		for(int N = 0; N < ensemble.NumberOfIons;N++)
		{
			PosX[N] = ((float) ensemble.ions[N].Pos[0]);
			PosY[N] = ((float) ensemble.ions[N].Pos[1]);
			PosZ[N] = ((float) ensemble.ions[N].Pos[2]);
		}

		// updating position, looping through ions and x-y-z
		//CoulombWrapper( ForceX, ForceY, ForceZ, PosX, PosY, PosZ, ensemble.GetNumberOfIons());
		FastCoulombWrapper( ForceX, ForceY, ForceZ, PosX, PosY, PosZ, ForceX_d, ForceY_d, ForceZ_d, PosX_d, PosY_d, PosZ_d, ensemble.GetNumberOfIons()); // calculating coulomb force on GPU with cuda


		TotalVzSquared = 0; 
		TotalVrSquared = 0;

		TotalVxSquared = 0; 
		TotalVySquared = 0; 
		TotalV = 0;
		// Burde der ikke være en TotalVrSquared = 0  her ?


		// Updating position, looping through ions and x-y-z
		for(int N=0; N < ensemble.NumberOfIons; N++)
		{


			// this is dim x
			double Ftot = Ftrap(ensemble, N, t-1, 0, Vrf, Vend) + ((double) ForceX[N]);
			//double Ftot = Ftrap(ensemble, N, t-1, 0, Vrf, Vend) + ((double) ForceX[N]) + Ffriction(ensemble, N, 0); //Mads: added friction
			TempPos[N][0]  = ensemble.ions[N].Pos[0] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[0];
			TempVel[N][0]  = ensemble.ions[N].Vel[0] + dt/ensemble.ions[N].m*Ftot;

			// for the velocity histogram
			ensemble.ions[N].Vavg[0] += TempVel[N][0];

			// this is dim y
			//Ftot = Ftrap(ensemble, N, t-1, 0, Vrf, Vend) + ((double) ForceY[N]) + Ffriction(ensemble, N, 1); //Mads: added friction
			Ftot = Ftrap(ensemble, N, t-1, 1, Vrf, Vend) + ((double) ForceY[N]);
			TempPos[N][1]  = ensemble.ions[N].Pos[1] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[1];
			TempVel[N][1]  = ensemble.ions[N].Vel[1] + dt/ensemble.ions[N].m*Ftot;
			
			// for the velocity histogram
			ensemble.ions[N].Vavg[1] += TempVel[N][1];

			
			// this is dim z
			//Ftot = Ftrap(ensemble, N, t-1, 0, Vrf, Vend) + ((double) ForceZ[N]) + Ffriction(ensemble, N, 2); //Mads: added friction
			Ftot = Ftrap(ensemble, N, t-1, 2, Vrf, Vend) + ((double) ForceZ[N]);
			TempPos[N][2]  = ensemble.ions[N].Pos[2] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[2];
			TempVel[N][2]  = ensemble.ions[N].Vel[2] + dt/ensemble.ions[N].m*Ftot;
			
			// for the velocity histogram
			ensemble.ions[N].Vavg[2] += TempVel[N][2];


			// Calculating V_secular in with the axis (z) and radial (r)
			ensemble.ions[N].VzSec = ensemble.ions[N].VzSec + TempVel[N][2];
			ensemble.ions[N].VrSec = ensemble.ions[N].VrSec + sqrt(pow(TempVel[N][0],2)+pow(TempVel[N][1],2));
			
			ensemble.ions[N].VxSec = ensemble.ions[N].VxSec + TempVel[N][0];
			ensemble.ions[N].VySec = ensemble.ions[N].VySec + TempVel[N][1];

			// The total secular velocities squared.
			TotalVzSquared += pow(TempVel[N][2],2);
			TotalVrSquared += pow(TempVel[N][1],2) + pow(TempVel[N][0],2); // This is not used ?
			
			TotalVxSquared += pow(TempVel[N][0],2);
			TotalVySquared += pow(TempVel[N][1],2);

			TotalV += pow(TempVel[N][1],2) + pow(TempVel[N][0],2) + pow(TempVel[N][2],2);
		}
		
		// Calculation of the total tempeture and temperature std. (in kelvin)
		TotalRadialTemp += ensemble.ions[0].m/Kb*TotalVrSquared/2/ensemble.NumberOfIons;
		TotalTempx 		+= ensemble.ions[0].m/Kb*TotalVxSquared/ensemble.NumberOfIons;
		TotalTempy 		+= ensemble.ions[0].m/Kb*TotalVySquared/ensemble.NumberOfIons;

		TotalTemp += ensemble.ions[0].m/Kb*TotalVzSquared/ensemble.NumberOfIons;
		TotalTempStd += pow(ensemble.ions[0].m/Kb*TotalVzSquared/ensemble.NumberOfIons,2);

		//cout << TotalTempx <<" "<<TotalTempy <<" "<<TotalTemp << endl; 


		Temperaturefile << t <<  ",  "
						<< ((TotalV/ensemble.NumberOfIons)* ensemble.ions[0].m)/(3*Kb) <<  ",  "
						<< ((TotalVzSquared/ensemble.NumberOfIons)* ensemble.ions[0].m)/(Kb) <<  ",  "
						<< ((TotalVrSquared/ensemble.NumberOfIons)* ensemble.ions[0].m)/(2*Kb) <<  ",  " 
						<< endl;


		// If a periode has been completed : time to calculated Vavg of ions and Vsec.
		if (((t-1) % ((int) StepsPrPeriode)) == 0 && t != 0)
		{
			for(int N=0; N < ensemble.NumberOfIons; N++)
			{
				// calculating secular velocity
				for(int dim=0; dim < 3; dim++)
					ensemble.ions[N].Vavg[dim] = ensemble.ions[N].Vavg[dim]/StepsPrPeriode;
				// setting ions secular velocity
				ensemble.ions[N].Vsec = sqrt(pow(ensemble.ions[N].Vavg[0],2) + pow(ensemble.ions[N].Vavg[1],2) + pow(ensemble.ions[N].Vavg[2],2));

				// THIS IS FOR SINGLE ION RESCALE
				ensemble.Vz_IonArray[N] = ensemble.ions[N].Vavg[2];
				ensemble.Vx_IonArray[N] = ensemble.ions[N].Vavg[0];
				ensemble.Vy_IonArray[N] = ensemble.ions[N].Vavg[1];


				// setting to zero
				for(int dim=0; dim < 3; dim++)
					ensemble.ions[N].Vavg[dim] = 0;
			}
		}

		// saving the time step
		for(int N=0; N < ensemble.NumberOfIons; N++)
		{
			for(int dim=0; dim < 3; dim++)
			{
				ensemble.ions[N].Pos[dim] = TempPos[N][dim];
				ensemble.ions[N].Vel[dim] = TempVel[N][dim];
			}
		}


		// only save last 30 time steps for the xyz file.
		if(t + StepsPrPeriode + 1 >= TimeSteps)
		{
			for(int N=1; N <= ensemble.NumberOfIons; N++)
				fprintf(f, "%s%d\t%f\t%f\t%f\n ", "Ca", N, ensemble.Position(0, N-1)*1e6, ensemble.Position(1, N-1)*1e6, ensemble.Position(2, N-1)*1e6);
		}

		// Only update the histograms after a certain time has passed.
		if (t > StartRecordingHistogram)
		{
			ensemble.UpdateHistogram();
			ensemble.UpdateVelocityHistogram();
		}

		// If a periode has been completed : Time to deal with velocity distribution.
		if (((t-1) % ((int) StepsPrPeriode)) == 0 && t != 0) // Dealing with velocity distribution
		{

			// Maybe find away to removed this, kinda useless.
			// AvgTemp is the 
			AvgTemp += ensemble.GetCurrentTemperature()/((double) TimeSteps)*StepsPrPeriode ;
			TempSquare += pow(ensemble.GetCurrentTemperature(),2)/((double) TimeSteps)*StepsPrPeriode;

			// Setting things to zero before we calculate.
			for (int i = 0; i < ensemble.NumberOfRings ; i ++ ) 
			{
				ensemble.VzSecRMS_in_ring[i]			= 0.0;
				ensemble.VrSecRMS_in_ring[i]			= 0.0;
				ensemble.NumberOfIonsInRadius[i]		= 0;
			}


			// THIS IS FOR RESCALE DONE BEFORE I TRIED TO FIX IT
			ensemble.VzSecRMS = 0;
			ensemble.VrSecRMS = 0;

			
			for(int N=0; N < ensemble.NumberOfIons; N++)
			{
				ensemble.VzSecRMS = ensemble.VzSecRMS +pow(ensemble.ions[N].VzSec,2);
				ensemble.VrSecRMS = ensemble.VrSecRMS +pow(ensemble.ions[N].VrSec,2);
			}


			ensemble.VzSecRMS = sqrt(ensemble.VzSecRMS/(ensemble.NumberOfIons*pow(StepsPrPeriode,2)));
			ensemble.VrSecRMS = sqrt(ensemble.VrSecRMS/(ensemble.NumberOfIons*pow(StepsPrPeriode,2)));




			// THIS PART IS FOR RADIALRESCAl
			// The main loop for finding which "ring" an ion is located in.
			for(int n = 0; n < ensemble.NumberOfIons; n++)
			{
				// Postion as a point in the histogram
				int Nx = ((int) ((ensemble.ions[n].Position(0))/PixelToDistance+((double) HistNx)/2)  - (HistNx/2));
				int Ny = ((int) ((ensemble.ions[n].Position(1))/PixelToDistance+((double) HistNy)/2)  - (HistNy/2));
				int Nz = ((int) ((ensemble.ions[n].Position(2))/PixelToDistance+((double) HistNz)/2)  - (HistNz/2));
		
				// The square of the position
				double Nxsquare = Nx*Nx;
				double Nysquare = Ny*Ny;
				double Nzsquare = Nz*Nz;

				// Length from center of histogram to the ion.
				double r = sqrt(Nxsquare + Nysquare + Nzsquare); 

				double bin_d = 150.0/ ensemble.NumberOfRings;  

				// Finding what ring the ions is in.
				for (int i = 1; i <= ensemble.NumberOfRings ; i ++ ) 
				{
					if ( r <= i*bin_d && r > (i-1)*bin_d)  // Svarer til histx/2 / 100 
					{
						ensemble.RadiiOfIons[n] = i;		  // What ring the ion was found in.
						ensemble.NumberOfIonsInRadius[i-1]++; //increse the number of ions found.
					}
				}
			}

			for(int N = 0; N < ensemble.NumberOfIons; N++)
			{
				int rOfIon = ensemble.RadiiOfIons[N]; // What radius is the current ion in, ie what ring.

				ensemble.VzSecRMS_in_ring[rOfIon-1] = ensemble.VzSecRMS_in_ring[rOfIon-1] + pow(ensemble.ions[N].VzSec,2); // adding the Temperature of each ion in their goup.
				ensemble.VrSecRMS_in_ring[rOfIon-1] = ensemble.VrSecRMS_in_ring[rOfIon-1] + pow(ensemble.ions[N].VrSec,2); // adding the Temperature of each ion in their goup.
			
				ensemble.ions[N].VzSec = 0;
				ensemble.ions[N].VrSec = 0;
			}

			// Calculated the avg temp.
			for (int i = 0; i < ensemble.NumberOfRings ; i++ ) 
			{
				if ( ensemble.NumberOfIonsInRadius[i] != 0)
				{
				ensemble.VzSecRMS_in_ring[i] = sqrt(ensemble.VzSecRMS_in_ring[i]/(ensemble.NumberOfIonsInRadius[i]*pow(StepsPrPeriode,2))); // Average by the number of ions in the ring.
				ensemble.VrSecRMS_in_ring[i] = sqrt(ensemble.VrSecRMS_in_ring[i]/(ensemble.NumberOfIonsInRadius[i]*pow(StepsPrPeriode,2)));
				}
			}
		}

		/* The version before
		// If a periode has been completed : Time to deal with velocity distribution.
		if (((t-1) % ((int) StepsPrPeriode)) == 0) // Dealing with velocity distribution
		{

			// AvgTemp is the 
			AvgTemp += ensemble.GetCurrentTemperature()*1000/((double) TimeSteps)*StepsPrPeriode ;
			TempSquare += pow(ensemble.GetCurrentTemperature()*1000,2)/((double) TimeSteps)*StepsPrPeriode;
			ensemble.VzSecRMS = 0;
			ensemble.VrSecRMS = 0;

			
			for(int N=0; N < ensemble.NumberOfIons; N++)
			{

				ensemble.VzSecRMS = ensemble.VzSecRMS +pow(ensemble.ions[N].VzSec,2);
				ensemble.VrSecRMS = ensemble.VrSecRMS +pow(ensemble.ions[N].VrSec,2);

				ensemble.ions[N].VzSec = 0;
				ensemble.ions[N].VrSec = 0;
			}


			ensemble.VzSecRMS = sqrt(ensemble.VzSecRMS/(ensemble.NumberOfIons*pow(StepsPrPeriode,2)));
			ensemble.VrSecRMS = sqrt(ensemble.VrSecRMS/(ensemble.NumberOfIons*pow(StepsPrPeriode,2)));

		
		}
		*/

		// If a periode has been completed : Time to deal with velocity distribution.
		/*
		if (((t-1) % ((int) StepsPrPeriode)) == 0 && t != 0) // Dealing with velocity distribution
		{
			ensemble.RescaleVelocityDistributionRadial();
		}
		*/
		
		if (((t-RandomNumber) % ((int) (3*StepsPrPeriode))) == 0)
		{
			RandomNumber = rand() % ((int) StepsPrPeriode-1) + 1;
			while ( ((2 < RandomNumber) && (13 > RandomNumber)) || ((17 < RandomNumber) && (28 > RandomNumber)))
				RandomNumber = rand() % ((int) StepsPrPeriode-1) + 1;
				

				// Time to rescale velocity
				//ensemble.RescaleVelocityDistributionRadial();
				ensemble.RescaleVelocityDistribution(); //So this might be the best
				//ensemble.RescaleVelocityDistributionIndivialIons();//
				
		}
		
		if ( t % 10000 == 0 ) {
			cout << "t = "<< t<< endl;
		}


	}

	// Har prøvet at rette så det er i Kelvin output.
	ensemble.SetActualTemperature(AvgTemp);
	ensemble.SetActualTemperatureSTD(sqrt(TempSquare-pow(AvgTemp,2)));
	TotalTemp = TotalTemp/((double) TimeSteps);
	TotalTempStd = sqrt(TotalTempStd/((double) TimeSteps)-pow(TotalTemp,2));
	
	cout << "Total z Velocity " << TotalTemp << "+/-" << TotalTempStd << endl;
	
	cout << "Total radial T " << TotalTempx  / ((double) TimeSteps) << "  "<< TotalTempy  / ((double) TimeSteps)<< endl;

	cout << "Secular z Velocity " << ensemble.GetActualTemperature() << '\n';

	cout << "Actual temp : "<< ensemble.GetActualTemperature() << " K  with STD: " << ensemble.GetActualTemperatureSTD() <<" K \n";


	// closing textfile
	Temperaturefile.close();
	fclose(f);
	// cleaning up on GPU
	CudaCoulombFree(PosX_d, PosY_d, PosZ_d, ForceX_d, ForceY_d, ForceZ_d);
	// cleaning up in general
	delete PosX;
	delete PosY;
	delete PosZ;
	delete ForceX;
	delete ForceY;
	delete ForceZ;
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		delete [] TempPos[i];
	delete TempPos;

	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		delete [] TempVel[i];
	delete TempVel;

}

// Testing Cuda Dynamic Temperature control integrator
void DynamicTemperatureLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend)
{


	// allocate memory to the histogram
	ensemble.InitialiseHistogram();

	ensemble.InitialiseVelocityHistogram();
	

	// allocate memory to Temperature Arrays
	ensemble.InitialiseTemperatureArrays(TimeSteps);
	

	// make Cuda Pos and Force arrays - its float because we dont want to overflow the GPU
	float *PosX, *PosY, *PosZ, *ForceX, *ForceY, *ForceZ;

	PosX = new float [ensemble.GetNumberOfIons()];
	PosY = new float [ensemble.GetNumberOfIons()];
	PosZ = new float [ensemble.GetNumberOfIons()];
	ForceX = new float [ensemble.GetNumberOfIons()];
	ForceY = new float [ensemble.GetNumberOfIons()];
	ForceZ = new float [ensemble.GetNumberOfIons()];

	double ***SecVel;
	SecVel = new double **[((int) StepsPrPeriode)];
	for(int i = 0; i < ((int) StepsPrPeriode);i++)
	{
		SecVel[i] = new double *[ensemble.NumberOfIons];
		for(int N = 0; N < ensemble.NumberOfIons;N++)
			SecVel[i][N] = new double [3];
	}

	double ***AvgSpeed;
	AvgSpeed = new double **[((int) StepsPrPeriode)];
	for(int i = 0; i < ((int) StepsPrPeriode);i++)
	{
		AvgSpeed[i] = new double *[ensemble.NumberOfIons];
		for(int N = 0; N < ensemble.NumberOfIons;N++)
			AvgSpeed[i][N] = new double [2];
	}

	double ***rmsVel;
	rmsVel = new double **[((int) StepsPrPeriode)];
	for(int i = 0; i < ((int) StepsPrPeriode);i++)
	{
		rmsVel[i] = new double *[ensemble.NumberOfIons];
		for(int N = 0; N < ensemble.NumberOfIons;N++)
			rmsVel[i][N] = new double [3];
	}

	// declare and allocate memory on GPU
	float *PosX_d, *PosY_d, *PosZ_d, *ForceX_d, *ForceY_d, *ForceZ_d;
	CudaCoulombAlloc(&PosX_d, &PosY_d, &PosZ_d, &ForceX_d, &ForceY_d, &ForceZ_d, ensemble.GetNumberOfIons());


	// making tempeary position and velocities
	double **TempPos;
	double **TempVel;


	TempPos = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempPos[i] = new double [3];

	TempVel = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempVel[i] = new double [3];



	for(int t=1; t < TimeSteps ; t++)
	{

		// Prepare Data for CUDA transfer
		for(int N = 0; N < ensemble.NumberOfIons;N++)
		{
			PosX[N] = ((float) ensemble.ions[N].Pos[0]);
			PosY[N] = ((float) ensemble.ions[N].Pos[1]);
			PosZ[N] = ((float) ensemble.ions[N].Pos[2]);

		}


		// updating position, looping through ions and x-y-z
		//CoulombWrapper( ForceX, ForceY, ForceZ, PosX, PosY, PosZ, ensemble.GetNumberOfIons());
		FastCoulombWrapper( ForceX, ForceY, ForceZ, PosX, PosY, PosZ, ForceX_d, ForceY_d, ForceZ_d, PosX_d, PosY_d, PosZ_d, ensemble.GetNumberOfIons()); // calculating coulomb force on GPU with cuda



		// updating position, looping through ions and x-y-z
		for(int N=0; N < ensemble.NumberOfIons; N++)
		{

			// this is dim x
			double Ftot = Ftrap(ensemble, N, t-1, 0, Vrf, Vend) + ((double) ForceX[N]) + Ffriction(ensemble, N, 0);

			//double Ftot = Fpseudo(ensemble, N, 0, Vrf, Vend) + ((double) ForceX[N]);
			TempPos[N][0]  = ensemble.ions[N].Pos[0] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[0];
			TempVel[N][0]  = ensemble.ions[N].Vel[0] + dt/ensemble.ions[N].m*Ftot;

			// this is dim y
			//Ftot = Fpseudo(ensemble, N, 2, Vrf, Vend) + ((double) ForceY[N]);
			Ftot = Ftrap(ensemble, N, t-1, 1, Vrf, Vend) + ((double) ForceY[N]) + Ffriction(ensemble, N, 1);
				
			TempPos[N][1]  = ensemble.ions[N].Pos[1] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[1];
			TempVel[N][1]  = ensemble.ions[N].Vel[1] + dt/ensemble.ions[N].m*Ftot;


			// this is dim z
			//Ftot = Fpseudo(ensemble, N, 3, Vrf, Vend) + ((double) ForceZ[N]);
			Ftot = Ftrap(ensemble, N, t-1, 2, Vrf, Vend) + ((double) ForceZ[N]) + Ffriction(ensemble, N, 2);
			
			TempPos[N][2]  = ensemble.ions[N].Pos[2] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[2];
			TempVel[N][2]  = ensemble.ions[N].Vel[2] + dt/ensemble.ions[N].m*Ftot;
		}

		// saving time step
		for(int N=0; N < ensemble.NumberOfIons; N++)
			for(int dim=0; dim < 3; dim++)
			{
				ensemble.ions[N].Pos[dim] = TempPos[N][dim];
				ensemble.ions[N].Vel[dim] = TempVel[N][dim];
			}

		// Calculate and Save Temperature
		for(int N=0; N < ensemble.NumberOfIons; N++)
		{
			// Calculating total temperature
			//ensemble.Ttotal[t] += (pow(ensemble.ions[N].Vel[0],2)+pow(ensemble.ions[N].Vel[1],2)+pow(ensemble.ions[N].Vel[2],2))*ensemble.Mass(0)/(3*Kb*ensemble.NumberOfIons);


			// Saving velocities in rf-periode
			for(int RfPhase = 0; RfPhase < ((int) StepsPrPeriode) ; RfPhase++)
			{
				for(int dim = 0; dim < 3; dim++)
				{
					SecVel[RfPhase][N][dim] += ensemble.ions[N].Vel[dim];
					rmsVel[RfPhase][N][dim] += pow(ensemble.ions[N].Vel[dim],2);
				}
				AvgSpeed[RfPhase][N][0] += sqrt(pow(ensemble.ions[N].Vel[0], 2) + pow(ensemble.ions[N].Vel[1], 2)); //avg radial speed
				AvgSpeed[RfPhase][N][1] += fabs(ensemble.ions[N].Vel[2]); //avg axial speed
			}
			// Calculating secular temperature
			//ensemble.Tsec[t] += (pow(SecVel[t % ((int) StepsPrPeriode)][0] ,2)+pow(SecVel[t % ((int) StepsPrPeriode)][1],2)+pow(SecVel[t % ((int) StepsPrPeriode)][2],2))*ensemble.Mass(0)/(3*Kb*ensemble.NumberOfIons*pow(StepsPrPeriode, 2));

		}




		// Calculating temperatures
		ensemble.Tsec[t] = 0;
		for(int N=0; N < ensemble.NumberOfIons; N++)
		{
			ensemble.Tsec[t] +=	(	pow( SecVel[t % ((int) StepsPrPeriode)][N][0] ,2) +
									pow( SecVel[t % ((int) StepsPrPeriode)][N][1] ,2) +
									pow( SecVel[t % ((int) StepsPrPeriode)][N][2] ,2) ) *
									ensemble.Mass(0)/(3*Kb*ensemble.NumberOfIons * pow(StepsPrPeriode, 2));
			
			// Disse 2 er de vigtige for rescale 
			ensemble.Tsecrad[t] += (pow( SecVel[t % ((int) StepsPrPeriode)][N][0] ,2) + 
									pow( SecVel[t % ((int) StepsPrPeriode)][N][1] ,2) ) * 
									ensemble.Mass(0)/(2*Kb*ensemble.NumberOfIons * pow(StepsPrPeriode, 2));

			ensemble.Tsecz[t] +=	pow( SecVel[t % ((int) StepsPrPeriode)][N][2] ,2) * 
									ensemble.Mass(0)/(Kb*ensemble.NumberOfIons * pow(StepsPrPeriode, 2));

			//

			ensemble.Trms[t] +=	(rmsVel[t % ((int) StepsPrPeriode)][N][0] +
								 rmsVel[t % ((int) StepsPrPeriode)][N][1] +
								 rmsVel[t % ((int) StepsPrPeriode)][N][2])* 
								 ensemble.Mass(0)/(3*Kb*ensemble.NumberOfIons*pow(StepsPrPeriode, 1));

			//ensemble.Trms[t] += ensemble.Mass(0)*pow(AvgSpeed[t % ((int) StepsPrPeriode)][N][0]/StepsPrPeriode, 2)/(2*ensemble.NumberOfIons*Kb); // temperature based on radial velocity
			ensemble.Trmsz[t] +=	(rmsVel[t % ((int) StepsPrPeriode)][N][2]) * 
									ensemble.Mass(0)/(Kb*ensemble.NumberOfIons*pow(StepsPrPeriode, 1));

		}

		// rescaling velocity distribution
		if(t % ((int) StepsPrPeriode) == 0 || ((int) StepsPrPeriode) == 15)
			ensemble.RescaleVelocityDistribution(t);

		// setting cell to zero again
		for(int N=0; N < ensemble.NumberOfIons; N++)
		{
			SecVel[t % ((int) StepsPrPeriode)][N][0] = 0;
			SecVel[t % ((int) StepsPrPeriode)][N][1] = 0;
			SecVel[t % ((int) StepsPrPeriode)][N][2] = 0;

			rmsVel[t % ((int) StepsPrPeriode)][N][0] = 0;
			rmsVel[t % ((int) StepsPrPeriode)][N][1] = 0;
			rmsVel[t % ((int) StepsPrPeriode)][N][2] = 0;

			AvgSpeed[t % ((int) StepsPrPeriode)][N][0] = 0;
			AvgSpeed[t % ((int) StepsPrPeriode)][N][1] = 0;

		}


		//cout << "total temperature "<< ensemble.GetTtotal(t) << '\n';
		//cout << "secular temperature " << ensemble.GetTsecular(t) << " rf phase at " << t % ((int) StepsPrPeriode) << " and time at " << t << '\n';

		// update 3d histogram
		if (t > StartRecordingHistogram)
		{
			ensemble.UpdateHistogram();
			ensemble.UpdateVelocityHistogram();
		}

		if ( t % 10000 == 0 ) {
			cout << "t = "<< t<< endl;
		}

	}


	// cleaning up on GPU
	CudaCoulombFree(PosX_d, PosY_d, PosZ_d, ForceX_d, ForceY_d, ForceZ_d);
	// cleaning up in general
	delete PosX;
	delete PosY;
	delete PosZ;
	delete ForceX;
	delete ForceY;
	delete ForceZ;
	// delete SecVel; remember to fix this
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		delete [] TempPos[i];
	delete TempPos;

	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		delete [] TempVel[i];
	delete TempVel;

}

void MADSDynamicTemperatureLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend)
{
	// allocate memory to the histogram
	ensemble.InitialiseHistogram();
	ensemble.InitialiseVelocityHistogram();
	ensemble.InitialiseCountHistogram();
	
	// make Cuda Pos and Force arrays - its float because we dont want to overflow the GPU
	float *PosX, *PosY, *PosZ, *ForceX, *ForceY, *ForceZ;

	PosX = new float [ensemble.GetNumberOfIons()];
	PosY = new float [ensemble.GetNumberOfIons()];
	PosZ = new float [ensemble.GetNumberOfIons()];
	ForceX = new float [ensemble.GetNumberOfIons()];
	ForceY = new float [ensemble.GetNumberOfIons()];
	ForceZ = new float [ensemble.GetNumberOfIons()];

	// Declare and allocate memory on GPU
	float *PosX_d, *PosY_d, *PosZ_d, *ForceX_d, *ForceY_d, *ForceZ_d;
	CudaCoulombAlloc(&PosX_d, &PosY_d, &PosZ_d, &ForceX_d, &ForceY_d, &ForceZ_d, ensemble.GetNumberOfIons());


	// Making tempeary position and velocities ( for use in calculation of new pos and vel after calculation on CUDA)
	double **TempPos;
	double **TempVel;

	TempPos = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempPos[i] = new double [3];

	TempVel = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempVel[i] = new double [3];


	// Array for storing the last position.
	double **PosOneTauAgo;
	PosOneTauAgo = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
	{
		PosOneTauAgo[i] = new double [3];
		for(int dim=0; dim < 3; dim++)
			{
				PosOneTauAgo[i][dim] = ensemble.ions[i].Pos[dim];
			}
	}

	double **PosTwoTauAgo;
	PosTwoTauAgo = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
	{
		PosTwoTauAgo[i] = new double [3];
		for(int dim=0; dim < 3; dim++)
			{
				PosTwoTauAgo[i][dim] = ensemble.ions[i].Pos[dim];
			}
	}

	// Total Vrms of the terminal speed
	double Total_V_r_rms = 0.0;
	double Total_V_z_rms = 0.0;
	double Total_V_x_rms = 0.0;
	double Total_V_y_rms = 0.0;
	// 
	ofstream Temperaturefile ("TemperatureData.txt");
	Temperaturefile << "First line of output" << endl ;
	ofstream Periodefile ("PeriodeData.txt");
	Periodefile << "First line of output" << endl ;

	for(int t=1; t <= TimeSteps ; t++)
	{

		// Prepare Data for CUDA transfer
		for(int N = 0; N < ensemble.GetNumberOfIons();N++)
		{
			PosX[N] = ((float) ensemble.ions[N].Pos[0]);
			PosY[N] = ((float) ensemble.ions[N].Pos[1]);
			PosZ[N] = ((float) ensemble.ions[N].Pos[2]);

		}
		// updating position, looping through ions and x-y-z
		//CoulombWrapper( ForceX, ForceY, ForceZ, PosX, PosY, PosZ, ensemble.GetNumberOfIons());
		FastCoulombWrapper( ForceX, ForceY, ForceZ, PosX, PosY, PosZ, ForceX_d, ForceY_d, ForceZ_d, PosX_d, PosY_d, PosZ_d, ensemble.GetNumberOfIons());
		// calculating coulomb force on GPU with cuda

		//  For output of temperatures.
		double TotalV = 0.0;
		double TotalVz = 0.0;
		double TotalVr = 0.0;
		double TotalVx = 0.0;
		double TotalVy = 0.0;

		// updating position, looping through ions and x-y-z
		for(int N=0; N < ensemble.GetNumberOfIons(); N++)
		{
			// This is dim X
			double Ftot = Ftrap(ensemble, N, t-1, 0, Vrf, Vend) + ((double) ForceX[N]);
			TempPos[N][0]  = ensemble.ions[N].Pos[0] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[0];
			TempVel[N][0]  = ensemble.ions[N].Vel[0] + dt/ensemble.ions[N].m*Ftot;

			// This is dim Y
			Ftot = Ftrap(ensemble, N, t-1, 1, Vrf, Vend) + ((double) ForceY[N]);
				
			TempPos[N][1]  = ensemble.ions[N].Pos[1] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[1];
			TempVel[N][1]  = ensemble.ions[N].Vel[1] + dt/ensemble.ions[N].m*Ftot;

			// This is dim Z
			Ftot = Ftrap(ensemble, N, t-1, 2, Vrf, Vend) + ((double) ForceZ[N]);
			
			TempPos[N][2]  = ensemble.ions[N].Pos[2] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[2];
			TempVel[N][2]  = ensemble.ions[N].Vel[2] + dt/ensemble.ions[N].m*Ftot;

			TotalV	+= pow(TempVel[N][0],2) + pow(TempVel[N][1],2) + pow(TempVel[N][2],2);
			TotalVz += pow(TempVel[N][2],2);
			TotalVr += pow(TempVel[N][0],2) + pow(TempVel[N][1],2);
			TotalVx += pow(TempVel[N][0],2);
			TotalVy += pow(TempVel[N][1],2);
		}

		TotalV  = sqrt(TotalV / ensemble.GetNumberOfIons()); // Now TotalV is the root mean squared of the total speed. (So avg speed of 1 ion)
		TotalVz = sqrt(TotalVz / ensemble.GetNumberOfIons());
		TotalVr = sqrt(TotalVr / ensemble.GetNumberOfIons());
		TotalVx = sqrt(TotalVx / ensemble.GetNumberOfIons());
		TotalVy = sqrt(TotalVy / ensemble.GetNumberOfIons());
		
		// Saving to temperaturefile.
		Temperaturefile << t <<  ",  "
						<< (pow(TotalV,2) * ensemble.ions[0].m)/(3*Kb)	<<  ",  "
						<< (pow(TotalVz,2)* ensemble.ions[0].m)/(Kb)	<<  ",  "
						<< (pow(TotalVx,2)* ensemble.ions[0].m)/(Kb)	<<  ",  " 
						<< (pow(TotalVy,2)* ensemble.ions[0].m)/(Kb)	<<  ",  " 
						<< endl;

		// saving time step
		for(int N=0; N < ensemble.GetNumberOfIons(); N++)
		{
			for(int dim=0; dim < 3; dim++)
			{
				ensemble.ions[N].Pos[dim] = TempPos[N][dim];
				ensemble.ions[N].Vel[dim] = TempVel[N][dim];
			}
		}
		

		// This is when 2*PI/omegaRF has passed. Which is when we calculate termalspeed.
		/*
		if ((t) % ((int) (StepsPrPeriode + 28))) == 0)
		{
			// Calculating the root mean square of termal vel's.
			double tau = (2*PI) / OmegaRF; 

			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				double TermalVx = (PosOneTauAgo[N][0] - ensemble.ions[N].Pos[0] ) / tau;
				double TermalVy	= (PosOneTauAgo[N][1] - ensemble.ions[N].Pos[1] ) / tau;
				
				double TermalVr	= sqrt( pow(TermalVx,2) + pow(TermalVy,2));
				Total_V_r_rms	+= pow(TermalVr,2);

				Total_V_z_rms	+= pow((PosOneTauAgo[N][2] - ensemble.ions[N].Pos[2] ) / tau,2); 
				Total_V_x_rms	+= pow((PosOneTauAgo[N][0] - ensemble.ions[N].Pos[0] ) / tau,2); 
				Total_V_y_rms	+= pow((PosOneTauAgo[N][1] - ensemble.ions[N].Pos[1] ) / tau,2); 
			}

			Total_V_z_rms = Total_V_z_rms / ensemble.GetNumberOfIons();
			Total_V_z_rms = sqrt(Total_V_z_rms);

			Total_V_r_rms = Total_V_r_rms / ensemble.GetNumberOfIons();
			Total_V_r_rms = sqrt(Total_V_r_rms);

			// Storing the postion, for next calculation
			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				for(int dim=0; dim < 3; dim++)
				{
					PosOneTauAgo[N][dim] = ensemble.ions[N].Pos[dim];
				}
			}
		}

		if (((t) % ((int) (StepsPrPeriode))) == 0)
		{
			ensemble.RescaleVelocityXYZ(Total_V_x_rms,Total_V_y_rms,Total_V_z_rms);
		}
		*/

		double LowestPartOfPeriode = floor((StepsPrPeriode/2) + 0.5);
		//if (((int)(t+ LowestPartOfPeriode) % ((int) (StepsPrPeriode))) == 0 && ((t) % ((int) (StepsPrPeriode))) != 0) 
		//if (((t) % ((int) (LowestPartOfPeriode))) == 0 && ((t) % ((int) (StepsPrPeriode))) != 0)
	    if (((int)(t+ LowestPartOfPeriode-1) % ((int) (StepsPrPeriode))) == 0 && ((t) % ((int) (StepsPrPeriode))) != 0)  // for 105	
		{
			// Storing the postion, for next calculation.
			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				for(int dim=0; dim < 3; dim++)
				{
					PosTwoTauAgo[N][dim] = ensemble.ions[N].Pos[dim];
				}
			}

			if (t >= StartRecordingHistogram)
			{
				ensemble.MyUpdateVelocityHistogram();
				ensemble.UpdateCountHistogram();
			


			for(int N = 0; N < ensemble.GetNumberOfIons(); N++)
			{
				if (Periodefile.is_open())
				{
				Periodefile << t								<< ",  " <<
								N + 1							<< ",  " <<
								ensemble.ions[N].GetMass()		<< ",  " <<
								ensemble.ions[N].Pos[0]			<< ",  " <<
								ensemble.ions[N].Pos[1]			<< ",  " <<
								ensemble.ions[N].Pos[2]			<< ",  " <<
								ensemble.ions[N].Vel[0]			<< ",  " <<
								ensemble.ions[N].Vel[1]			<< ",  " <<
								ensemble.ions[N].Vel[2]			<< ",  " <<
								ensemble.ions[N].Velocity()		<< ",  " <<
								(pow(TotalV,2) * ensemble.ions[0].m)/(3*Kb)	<< ",  " <<
								endl;
				}   
			}
			}

		}

		if (((t) % ((int) (StepsPrPeriode))) == 0)
		{

			/*
			// Calculating the root mean square of termal vel's.
			double tau = (2*PI) / OmegaRF; 

			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				double TermalVx = (PosOneTauAgo[N][0] - ensemble.ions[N].Pos[0] ) / tau;
				double TermalVy	= (PosOneTauAgo[N][1] - ensemble.ions[N].Pos[1] ) / tau;
				
				double TermalVr	= sqrt( pow(TermalVx,2) + pow(TermalVy,2));
				Total_V_r_rms	+= pow(TermalVr,2);

				Total_V_z_rms	+= pow((PosOneTauAgo[N][2] - ensemble.ions[N].Pos[2] ) / tau,2); 
				Total_V_x_rms	+= pow((PosOneTauAgo[N][0] - ensemble.ions[N].Pos[0] ) / tau,2); 
				Total_V_y_rms	+= pow((PosOneTauAgo[N][1] - ensemble.ions[N].Pos[1] ) / tau,2); 
			}

			Total_V_z_rms = Total_V_z_rms / ensemble.GetNumberOfIons();
			Total_V_z_rms = sqrt(Total_V_z_rms);

			Total_V_r_rms = Total_V_r_rms / ensemble.GetNumberOfIons();
			Total_V_r_rms = sqrt(Total_V_r_rms);

			// Storing the postion, for next calculation
			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				for(int dim=0; dim < 3; dim++)
				{
					PosOneTauAgo[N][dim] = ensemble.ions[N].Pos[dim];
				}
			}
			*/

						// Calculating the root mean square of termal vel's.
			double tau = (2*PI) / OmegaRF; 

			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				double TermalVx = (PosTwoTauAgo[N][0] - PosOneTauAgo[N][0] ) / tau;
				double TermalVy	= (PosTwoTauAgo[N][1] - PosOneTauAgo[N][1] ) / tau;
				
				double TermalVr	= sqrt( pow(TermalVx,2) + pow(TermalVy,2));
				Total_V_r_rms	+= pow(TermalVr,2);

				Total_V_z_rms	+= pow((PosTwoTauAgo[N][2] - PosOneTauAgo[N][2] ) / tau,2); 
				Total_V_x_rms	+= pow((PosTwoTauAgo[N][0] - PosOneTauAgo[N][0] ) / tau,2); 
				Total_V_y_rms	+= pow((PosTwoTauAgo[N][1] - PosOneTauAgo[N][1] ) / tau,2); 
			}

			Total_V_z_rms = Total_V_z_rms / ensemble.GetNumberOfIons();
			Total_V_z_rms = sqrt(Total_V_z_rms);

			Total_V_x_rms = Total_V_x_rms / ensemble.GetNumberOfIons();
			Total_V_x_rms = sqrt(Total_V_x_rms);

			Total_V_y_rms = Total_V_y_rms / ensemble.GetNumberOfIons();
			Total_V_y_rms = sqrt(Total_V_y_rms);

			Total_V_r_rms = Total_V_r_rms / ensemble.GetNumberOfIons();
			Total_V_r_rms = sqrt(Total_V_r_rms);

			// Storing the postion, for next calculation
			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				for(int dim=0; dim < 3; dim++)
				{
					PosOneTauAgo[N][dim] = PosTwoTauAgo[N][dim];
				}
			}

	

			// Rescale velocities.
			ensemble.RescaleVelocityXYZ(Total_V_x_rms,Total_V_y_rms,Total_V_z_rms);

		}	
		// update 3d histogram
		
		if (t >= StartRecordingHistogram)
		{
			if(t == StartRecordingHistogram)
			{
				cout << "Record started " << endl; 
			}
			ensemble.UpdateHistogram();
		}
		
		if ( t % 10000 == 0 ) {
			cout << "****** t = "<< t<< " *****"<< endl;
		}

	}

	Periodefile.close();
	Temperaturefile.close();

	// cleaning up on GPU
	CudaCoulombFree(PosX_d, PosY_d, PosZ_d, ForceX_d, ForceY_d, ForceZ_d);
	// cleaning up in general
	delete PosX;
	delete PosY;
	delete PosZ;
	delete ForceX;
	delete ForceY;
	delete ForceZ;
	// delete SecVel; remember to fix this
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		delete [] TempPos[i];
	delete TempPos;

	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		delete [] TempVel[i];
	delete TempVel;

}

void OLDCudaLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend)
{
	// allocate memory to the histogram
	ensemble.InitialiseHistogram();
	ensemble.InitialiseVelocityHistogram();

	// make Cuda Pos and Force arrays - its float because we dont want to overflow the GPU
	float *PosX, *PosY, *PosZ, *ForceX, *ForceY, *ForceZ;

	PosX = new float [ensemble.GetNumberOfIons()];
	PosY = new float [ensemble.GetNumberOfIons()];
	PosZ = new float [ensemble.GetNumberOfIons()];
	ForceX = new float [ensemble.GetNumberOfIons()];
	ForceY = new float [ensemble.GetNumberOfIons()];
	ForceZ = new float [ensemble.GetNumberOfIons()];

	// declare and allocate memory on GPU
	float *PosX_d, *PosY_d, *PosZ_d, *ForceX_d, *ForceY_d, *ForceZ_d;
	CudaCoulombAlloc(&PosX_d, &PosY_d, &PosZ_d, &ForceX_d, &ForceY_d, &ForceZ_d, ensemble.GetNumberOfIons());


	// making tempeary position and velocities
	double **TempPos;
	double **TempVel;


	TempPos = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempPos[i] = new double [3];

	TempVel = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempVel[i] = new double [3];


	int counter = 1;

	double AvgTemp=0;
	double TempSquare=0;
	int RandomNumber = 1;
	double TotalVzSquared = 0;
	double TotalVrSquared = 0;
	double TotalTemp = 0;
	double TotalTempStd = 0;

	// writing out to textfile
	FILE *f;

	char FilNavn[14];
	_snprintf(FilNavn, 14, "%s%d%s", "MDPos", ensemble.NumberOfIons, ".xyz"); //Change to _snprintf

	f=fopen(FilNavn, "w");
	// writing start format to xyz file
	fprintf(f, "%d\n%s\n", ensemble.NumberOfIons, "Ca");

	for(int t=1; t < TimeSteps ; t++)
	{

		// Prepare Data for CUDA transfer
		for(int N = 0; N < ensemble.NumberOfIons;N++)
		{
			PosX[N] = ((float) ensemble.ions[N].Pos[0]);
			PosY[N] = ((float) ensemble.ions[N].Pos[1]);
			PosZ[N] = ((float) ensemble.ions[N].Pos[2]);
		}

		// updating position, looping through ions and x-y-z
		//CoulombWrapper( ForceX, ForceY, ForceZ, PosX, PosY, PosZ, ensemble.GetNumberOfIons());
		FastCoulombWrapper( ForceX, ForceY, ForceZ, PosX, PosY, PosZ, ForceX_d, ForceY_d, ForceZ_d, PosX_d, PosY_d, PosZ_d, ensemble.GetNumberOfIons()); // calculating coulomb force on GPU with cuda


		TotalVzSquared = 0;
		// updating position, looping through ions and x-y-z
		for(int N=0; N < ensemble.NumberOfIons; N++)
		{


			// this is dim x
			double Ftot = Ftrap(ensemble, N, t-1, 0, Vrf, Vend) + ((double) ForceX[N]);
			//double Ftot = Fpseudo(ensemble, N, 0, Vrf, Vend) + ((double) ForceX[N]);
			TempPos[N][0]  = ensemble.ions[N].Pos[0] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[0];
			TempVel[N][0]  = ensemble.ions[N].Vel[0] + dt/ensemble.ions[N].m*Ftot;
			// for the velocity histogram
			ensemble.ions[N].Vavg[0] += TempVel[N][0];

			// this is dim y
			//Ftot = Fpseudo(ensemble, N, 2, Vrf, Vend) + ((double) ForceY[N]);
			Ftot = Ftrap(ensemble, N, t-1, 1, Vrf, Vend) + ((double) ForceY[N]);
			TempPos[N][1]  = ensemble.ions[N].Pos[1] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[1];
			TempVel[N][1]  = ensemble.ions[N].Vel[1] + dt/ensemble.ions[N].m*Ftot;
			// for the velocity histogram
			ensemble.ions[N].Vavg[1] += TempVel[N][1];

			// this is dim z
			//Ftot = Fpseudo(ensemble, N, 3, Vrf, Vend) + ((double) ForceZ[N]);
			Ftot = Ftrap(ensemble, N, t-1, 2, Vrf, Vend) + ((double) ForceZ[N]);
			TempPos[N][2]  = ensemble.ions[N].Pos[2] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[2];
			TempVel[N][2]  = ensemble.ions[N].Vel[2] + dt/ensemble.ions[N].m*Ftot;
			ensemble.ions[N].Vavg[2] += TempVel[N][2];


			ensemble.ions[N].VzSec = ensemble.ions[N].VzSec + TempVel[N][2];
			ensemble.ions[N].VrSec = ensemble.ions[N].VrSec + sqrt(pow(TempVel[N][0],2)+pow(TempVel[N][1],2));
			TotalVzSquared += pow(TempVel[N][2],2);
			TotalVrSquared += pow(TempVel[N][1],2)+pow(TempVel[N][0],2);

		}

		TotalTemp += ensemble.ions[0].m/Kb*TotalVzSquared/ensemble.NumberOfIons;
		TotalTempStd += pow(ensemble.ions[0].m/Kb*TotalVzSquared/ensemble.NumberOfIons,2);

		if (((t-1) % ((int) StepsPrPeriode)) == 0)
			for(int N=0; N < ensemble.NumberOfIons; N++)
			{
				// calculating secular velocity
				for(int dim=0; dim < 3; dim++)
					ensemble.ions[N].Vavg[dim] = ensemble.ions[N].Vavg[dim]/StepsPrPeriode;
				// setting ions secular velocity
				ensemble.ions[N].Vsec = sqrt(pow(ensemble.ions[N].Vavg[0],2) + pow(ensemble.ions[N].Vavg[1],2) + pow(ensemble.ions[N].Vavg[2],2));

				// setting to zero
				for(int dim=0; dim < 3; dim++)
					ensemble.ions[N].Vavg[dim] = 0;
			}


		// saving time step
		for(int N=0; N < ensemble.NumberOfIons; N++)
			for(int dim=0; dim < 3; dim++)
			{
				ensemble.ions[N].Pos[dim] = TempPos[N][dim];
				ensemble.ions[N].Vel[dim] = TempVel[N][dim];
			}

		// only save last 30 time steps...
		if(t + StepsPrPeriode + 1 >= TimeSteps)
			for(int N=1; N <= ensemble.NumberOfIons; N++)
				fprintf(f, "%s%d\t%f\t%f\t%f\n ", "Ca", N, ensemble.Position(0, N-1)*1e6, ensemble.Position(1, N-1)*1e6, ensemble.Position(2, N-1)*1e6);


		if (t > StartRecordingHistogram)
		{
			ensemble.UpdateHistogram();
			ensemble.UpdateVelocityHistogram();
		}

		if (((t-1) % ((int) StepsPrPeriode)) == 0) // Dealing with velocity distribution
		{

			AvgTemp += ensemble.GetCurrentTemperature()/((double) TimeSteps)*StepsPrPeriode ;
			TempSquare += pow(ensemble.GetCurrentTemperature(),2)/((double) TimeSteps)*StepsPrPeriode;
			ensemble.VzSecRMS = 0;
			ensemble.VrSecRMS = 0;
			for(int N=0; N < ensemble.NumberOfIons; N++)
			{
				ensemble.VzSecRMS = ensemble.VzSecRMS + pow(ensemble.ions[N].VzSec,2);
				ensemble.VrSecRMS = ensemble.VrSecRMS + pow(ensemble.ions[N].VrSec,2);
				ensemble.ions[N].VzSec = 0;
				ensemble.ions[N].VrSec = 0;
			}
			ensemble.VzSecRMS = sqrt(ensemble.VzSecRMS/(ensemble.NumberOfIons*pow(StepsPrPeriode,2)));
			ensemble.VrSecRMS = sqrt(ensemble.VrSecRMS/(ensemble.NumberOfIons*pow(StepsPrPeriode,2)));


			for(int N=0; N < ensemble.NumberOfIons; N++)
				ensemble.ions[N].VzSec = 0;
			for(int N=0; N < ensemble.NumberOfIons; N++)
				ensemble.ions[N].VrSec = 0;
		}

		if (((t-RandomNumber) % ((int) (3*StepsPrPeriode))) == 0)
		{
			RandomNumber = rand() % ((int) StepsPrPeriode-1) + 1;
			while ( ((2 < RandomNumber) && (13 > RandomNumber)) || ((17 < RandomNumber) && (28 > RandomNumber)))
				RandomNumber = rand() % ((int) StepsPrPeriode-1) + 1;
			ensemble.RescaleVelocityDistribution();
		}


	}

	ensemble.SetActualTemperature(AvgTemp);
	ensemble.SetActualTemperatureSTD(sqrt(TempSquare-pow(AvgTemp,2)));
	TotalTemp = TotalTemp/((double) TimeSteps);
	TotalTempStd = sqrt(TotalTempStd/((double) TimeSteps)-pow(TotalTemp,2));
	cout << "Total z Velocity " << TotalTemp << "+/-" << TotalTempStd << "Secular z Velocity " << ensemble.GetActualTemperature() << '\n';
	cout << ensemble.GetActualTemperature() << '\t' << ensemble.GetActualTemperatureSTD() <<'\n';
	// closing textfile
	fclose(f);
	// cleaning up on GPU
	CudaCoulombFree(PosX_d, PosY_d, PosZ_d, ForceX_d, ForceY_d, ForceZ_d);
	// cleaning up in general
	delete PosX;
	delete PosY;
	delete PosZ;
	delete ForceX;
	delete ForceY;
	delete ForceZ;
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		delete [] TempPos[i];
	delete TempPos;

	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		delete [] TempVel[i];
	delete TempVel;

}

void TauPeriodeCudaLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend)
{
	// allocate memory to the histogram
	ensemble.InitialiseHistogram();
	ensemble.InitialiseVelocityHistogram();
	ensemble.InitialiseCountHistogram();
	
	// make Cuda Pos and Force arrays - its float because we dont want to overflow the GPU
	float *PosX, *PosY, *PosZ, *ForceX, *ForceY, *ForceZ;

	PosX = new float [ensemble.GetNumberOfIons()];
	PosY = new float [ensemble.GetNumberOfIons()];
	PosZ = new float [ensemble.GetNumberOfIons()];
	ForceX = new float [ensemble.GetNumberOfIons()];
	ForceY = new float [ensemble.GetNumberOfIons()];
	ForceZ = new float [ensemble.GetNumberOfIons()];

	// Declare and allocate memory on GPU
	float *PosX_d, *PosY_d, *PosZ_d, *ForceX_d, *ForceY_d, *ForceZ_d;
	CudaCoulombAlloc(&PosX_d, &PosY_d, &PosZ_d, &ForceX_d, &ForceY_d, &ForceZ_d, ensemble.GetNumberOfIons());


	// Making tempeary position and velocities ( for use in calculation of new pos and vel after calculation on CUDA)
	double **TempPos;
	double **TempVel;

	TempPos = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempPos[i] = new double [3];

	TempVel = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempVel[i] = new double [3];


	// Array for storing the last position.
	double **PosOneTauAgo;
	PosOneTauAgo = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
	{
		PosOneTauAgo[i] = new double [3];
		for(int dim=0; dim < 3; dim++)
			{
				PosOneTauAgo[i][dim] = ensemble.ions[i].Pos[dim];
			}
	}

	double **PosTwoTauAgo;
	PosTwoTauAgo = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
	{
		PosTwoTauAgo[i] = new double [3];
		for(int dim=0; dim < 3; dim++)
			{
				PosTwoTauAgo[i][dim] = ensemble.ions[i].Pos[dim];
			}
	}

	// Total Vrms of the terminal speed
	double Total_V_r_rms = 0.0;
	double Total_V_z_rms = 0.0;
	double Total_V_x_rms = 0.0;
	double Total_V_y_rms = 0.0;
	// 
	ofstream Temperaturefile ("TemperatureData.txt");
	Temperaturefile << "First line of output" << endl ;
	ofstream Periodefile ("PeriodeData.txt");
	Periodefile << "First line of output" << endl ;

	for(int t=1; t <= TimeSteps ; t++)
	{

		// Prepare Data for CUDA transfer
		for(int N = 0; N < ensemble.GetNumberOfIons();N++)
		{
			PosX[N] = ((float) ensemble.ions[N].Pos[0]);
			PosY[N] = ((float) ensemble.ions[N].Pos[1]);
			PosZ[N] = ((float) ensemble.ions[N].Pos[2]);

		}
		// updating position, looping through ions and x-y-z
		//CoulombWrapper( ForceX, ForceY, ForceZ, PosX, PosY, PosZ, ensemble.GetNumberOfIons());
		FastCoulombWrapper( ForceX, ForceY, ForceZ, PosX, PosY, PosZ, ForceX_d, ForceY_d, ForceZ_d, PosX_d, PosY_d, PosZ_d, ensemble.GetNumberOfIons());
		// calculating coulomb force on GPU with cuda

		//  For output of temperatures.
		double TotalV_squared = 0.0;
		double TotalVx_squared = 0.0;
		double TotalVy_squared = 0.0;
		double TotalVz_squared = 0.0;

		// updating position, looping through ions and x-y-z
		for(int N=0; N < ensemble.GetNumberOfIons(); N++)
		{
			// This is dim X
			double Ftot = Ftrap(ensemble, N, t-1, 0, Vrf, Vend) + ((double) ForceX[N]);
			TempPos[N][0]  = ensemble.ions[N].Pos[0] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[0];
			TempVel[N][0]  = ensemble.ions[N].Vel[0] + dt/ensemble.ions[N].m*Ftot;

			// This is dim Y
			Ftot = Ftrap(ensemble, N, t-1, 1, Vrf, Vend) + ((double) ForceY[N]);
				
			TempPos[N][1]  = ensemble.ions[N].Pos[1] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[1];
			TempVel[N][1]  = ensemble.ions[N].Vel[1] + dt/ensemble.ions[N].m*Ftot;

			// This is dim Z
			Ftot = Ftrap(ensemble, N, t-1, 2, Vrf, Vend) + ((double) ForceZ[N]);
			
			TempPos[N][2]  = ensemble.ions[N].Pos[2] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[2];
			TempVel[N][2]  = ensemble.ions[N].Vel[2] + dt/ensemble.ions[N].m*Ftot;

			TotalV_squared	+= pow(TempVel[N][0],2) + pow(TempVel[N][1],2) + pow(TempVel[N][2],2);
			TotalVz_squared += pow(TempVel[N][2],2);
			TotalVx_squared += pow(TempVel[N][0],2);
			TotalVx_squared += pow(TempVel[N][1],2);
		}
		
		// Saving to temperaturefile.
		Temperaturefile << t																<<  ",  "
						<< (TotalV_squared  * ensemble.ions[0].m)/(3*Kb*ensemble.GetNumberOfIons())	<<  ",  "
						<< (TotalVz_squared * ensemble.ions[0].m)/(Kb* ensemble.GetNumberOfIons())	<<  ",  "
						<< (TotalVx_squared * ensemble.ions[0].m)/(Kb* ensemble.GetNumberOfIons())	<<  ",  " 
						<< (TotalVz_squared * ensemble.ions[0].m)/(Kb* ensemble.GetNumberOfIons())	<<  ",  " 
						<< endl;
		


		// saving time step
		for(int N=0; N < ensemble.GetNumberOfIons(); N++)
		{
			for(int dim=0; dim < 3; dim++)
			{
				ensemble.ions[N].Pos[dim] = TempPos[N][dim];
				ensemble.ions[N].Vel[dim] = TempVel[N][dim];
			}
		}
		
		double LowestPartOfPeriode = floor((StepsPrPeriode/2) + 0.5);

		if (((t) % ((int) (LowestPartOfPeriode))) == 0 && ((t) % ((int) (StepsPrPeriode))) != 0)
		{
			// Storing the postion, for next calculation.
			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				for(int dim=0; dim < 3; dim++)
				{
					PosTwoTauAgo[N][dim] = ensemble.ions[N].Pos[dim];
				}
			}

			if (t >= StartRecordingHistogram)
			{
				ensemble.MyUpdateVelocityHistogram();
				ensemble.UpdateCountHistogram();
			}

			/*
			for(int N = 0; N < ensemble.GetNumberOfIons(); N++)
			{
				if (Periodefile.is_open())
				{
				Periodefile <<			N + 1					<< ",  " <<
								ensemble.ions[N].GetMass()		<< ",  " <<
								ensemble.ions[N].Pos[0]			<< ",  " <<
								ensemble.ions[N].Pos[1]			<< ",  " <<
								ensemble.ions[N].Pos[2]			<< ",  " <<
								ensemble.ions[N].Vel[0]			<< ",  " <<
								ensemble.ions[N].Vel[1]			<< ",  " <<
								ensemble.ions[N].Vel[2]			<< ",  " <<
								ensemble.ions[N].Velocity()		<< ",  " <<
								(pow(TotalV,2) * ensemble.ions[0].m)/(3*Kb)	<< ",  " <<
								endl;
				}   
			}
			*/
		}

		if (((t) % ((int) (StepsPrPeriode))) == 0)
		{
			double tau = (2*PI) / OmegaRF; 

			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{			
				Total_V_z_rms	+= pow((PosTwoTauAgo[N][2] - PosOneTauAgo[N][2] ) / tau,2); 
				Total_V_x_rms	+= pow((PosTwoTauAgo[N][0] - PosOneTauAgo[N][0] ) / tau,2); 
				Total_V_y_rms	+= pow((PosTwoTauAgo[N][1] - PosOneTauAgo[N][1] ) / tau,2); 
			}

			// Storing the postion, for next calculation
			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				for(int dim=0; dim < 3; dim++)
				{
					PosOneTauAgo[N][dim] = PosTwoTauAgo[N][dim];
				}
			}

			// Rescale velocities.
			ensemble.RescaleVelocityXYZ(Total_V_x_rms,Total_V_y_rms,Total_V_z_rms);

		}	
		// update 3d histogram
		
		if (t >= StartRecordingHistogram)
		{
			if(t == StartRecordingHistogram)
			{
				cout << "Record started " << endl; 
			}
			ensemble.UpdateHistogram();
		}
		
		if ( t % 10000 == 0 ) {
			cout << "***** t = "<< t<< " *********"<<endl;
		}

	}

	Periodefile.close();
	Temperaturefile.close();

	// cleaning up on GPU
	CudaCoulombFree(PosX_d, PosY_d, PosZ_d, ForceX_d, ForceY_d, ForceZ_d);
	// cleaning up in general
	delete PosX;
	delete PosY;
	delete PosZ;
	delete ForceX;
	delete ForceY;
	delete ForceZ;
	// delete SecVel; remember to fix this
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		delete [] TempPos[i];
	delete TempPos;

	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		delete [] TempVel[i];
	delete TempVel;

}