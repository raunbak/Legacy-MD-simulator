// MolecularDynamics.cpp : Defines the entry point for the console application.
//

// include from making project in VS-10
#include "stdafx.h"
#include <stdio.h>

// From the original MD.cpp
#include <cstdlib>
#include <cstring>
#include <iostream>
#include "FastEnsemble.h"
#include "integrator.h"
#include "constants.h"

// To save files.
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <stdarg.h>
using namespace std;


// unix main

int main()
{
clock_t tStart = clock();

	double Vrf = 220;
	double Vend = 3.1; //Kugle = 18. Oval = 3.1
	// HUSK AT TJEKKE starthist når timesteps ændres! (LAV DET AUTO))
	int TimeSteps = 550000;//100000; //  Scale : 10000 = 0.1ms 
	double Temperature = 0.00671;	 // Scale is 0.01 = 10mK eller 0.001 = 1mK	    


	cout << "Allocating memory...\n";
	//int m1, int n1, int m2, int n2
	//FastEnsemble crystal(40,750,40,750);
	FastEnsemble crystal(40,500,40,500);
	
	cout << "Generating crystal...\n";
	crystal.CrystalGenerator(Vrf,Vend);
	crystal.SetSteadyStateTemperature(Temperature);

	cout << crystal.GetNumberOfIons();
	// running sim
	cout << " ions using steps with length = " << dt << "s\n";
	
	MADSDynamicTemperatureLeFrogintegrator(crystal, TimeSteps, Vrf, Vend);
	//DynamicTemperatureLeFrogintegrator(crystal, TimeSteps, Vrf, Vend);


	//OLDCudaLeFrogintegrator(crystal,TimeSteps,Vrf,Vend);
	//CudaLeFrogintegrator(crystal,TimeSteps,Vrf,Vend);

	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	cout << "DONE with simulation\n";
	cout << "Now saving to datafiles\n";
	
	// Saves histogram data to file
	ofstream Histogramfile ("HistogramData.txt");
	Histogramfile << "Bin  i j k" << endl ;

	ofstream VHistfile ("VelocityHistogramData.txt");
	VHistfile << "Sum of all velocity in this bin  i \t j \t k" << endl ;

	ofstream CHistfile ("CountHistogramData.txt");
	CHistfile << "bin  i j k" << endl ;

	for(int i=0; i < HistNx; i++)
	{
		for(int j=0; j < HistNy; j++)
		{
			for(int k=0; k < HistNz; k++)
			{
				if (Histogramfile.is_open())
				{
					Histogramfile	<< crystal.ReturnHist(i, j, k)				<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}
				if (VHistfile.is_open())
				{
					VHistfile		<< crystal.ReturnVelHist(i, j, k)			<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}
				if (CHistfile.is_open())
				{
					CHistfile		<< crystal.ReturnCountHist(i, j, k)			<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}
			}
		}
	}

	Histogramfile.close();
	VHistfile.close();
	CHistfile.close();
    
	cout << "Histograms done" << endl;

	crystal.SaveIonDataToFile();


	cout << "Ion data done" << endl;
	// cleaning up
    //crystal.FreeTemperatureArrays(); // May haveing this not in helps.
	
	cout << "Arrays freed" << endl;

	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	cout << "Quest completed! I mean data saved!\n";
	cout << "Press a key to end the program\n";


	cin.get();
	
	


	return 0;
}
