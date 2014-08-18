#include "stdafx.h" // Pre-compiled header.

#include <cstdlib>
#include <cmath>
#include <cstring>
//#include <mex.h>      // SOMETHING IS WROOOONG
#include <iostream>
#include <complex>
#include "Ensemble.h"
#include "integrator.h"
#include "constants.h"

using namespace std;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	// testing new temperature control
	if (nrhs == 0 && nlhs == 6)
	{
		int N = 800;
		double Vrf = 210;
		double Vend = 2.8;
		int TimeSteps = 400000;
		double Temperature = 0.035;

		mexPrintf("Allocating memory...\n");
		FastEnsemble crystal(N);
		mexPrintf("Generating crystal...\n");
		crystal.CrystalGenerator(Vrf,Vend);
		crystal.SetSteadyStateTemperature(Temperature);

		// running sim
		DynamicTemperatureLeFrogintegrator(crystal, TimeSteps, Vrf, Vend);

		mexPrintf("Copying data to matlab...\n");
		const mwSize ImDims[]={HistNx,HistNy,HistNz};
		const mwSize Tempdims[]={1,TimeSteps};

		double *HistOut;
		double *Ttotalout;
		double *Trmszout;
		double *Tsecularout;
		double *Tsecularradout;
		double *Tsecularzout;

		plhs[0] = mxCreateNumericArray(3,ImDims,mxDOUBLE_CLASS,mxREAL);
		HistOut = mxGetPr(plhs[0]);
		plhs[1] = mxCreateNumericArray(2,Tempdims,mxDOUBLE_CLASS,mxREAL);
		Ttotalout = mxGetPr(plhs[1]);
		plhs[2] = mxCreateNumericArray(2,Tempdims,mxDOUBLE_CLASS,mxREAL);
		Trmszout = mxGetPr(plhs[2]);
		plhs[3] = mxCreateNumericArray(2,Tempdims,mxDOUBLE_CLASS,mxREAL);
		Tsecularout = mxGetPr(plhs[3]);
		plhs[4] = mxCreateNumericArray(2,Tempdims,mxDOUBLE_CLASS,mxREAL);
		Tsecularradout = mxGetPr(plhs[4]);
		plhs[5] = mxCreateNumericArray(2,Tempdims,mxDOUBLE_CLASS,mxREAL);
		Tsecularzout = mxGetPr(plhs[5]);


		for(int i=0; i < HistNx; i++)
			for(int j=0; j < HistNy; j++)
				for(int k=0; k < HistNz; k++)
					HistOut[k * HistNx*HistNy + j*HistNx + i] = ((double) crystal.ReturnHist(i, j, k));

		for(int t = 0; t < TimeSteps; t++)
		{
			//cout << "Total temperature " << crystal.GetTtotal(t) << '\n';
			Ttotalout[t] = crystal.GetTrms(t);
		}

		for(int t = 0; t < TimeSteps; t++)
		{
			//cout << "RMS z temperature " << crystal.GetTtotal(t) << '\n';
			Trmszout[t] = crystal.GetTrmsz(t);
		}

		for(int t = 0; t < TimeSteps; t++)
		{
			//cout << "Secular temperature " << crystal.GetTsecular(t) << '\n';
			Tsecularout[t] = crystal.GetTsecular(t);
		}

		for(int t = 0; t < TimeSteps; t++)
		{
			//cout << "Secular temperature " << crystal.GetTsecular(t) << '\n';
			Tsecularradout[t] = crystal.GetTsecularrad(t);
		}

		for(int t = 0; t < TimeSteps; t++)
		{
			//cout << "Secular temperature " << crystal.GetTsecular(t) << '\n';
			Tsecularzout[t] = crystal.GetTsecularz(t);
		}

		// cleaning up
		crystal.FreeTemperatureArrays();
		//mxDestroyArray(plhs[1]);
		//mxDestroyArray(plhs[2]);
	}

	if (nlhs == 6 && nrhs == 5)
	{
		// input variables...
		int N = mxGetScalar(prhs[0]);
		double Vrf = mxGetScalar(prhs[1]);
		double Vend = mxGetScalar(prhs[2]);
		int TimeSteps = mxGetScalar(prhs[3]);
		double Temperature = mxGetScalar(prhs[4]);

		// initializing stuff...
		mexPrintf("Allocating memory...\n");
		FastEnsemble crystal(N);
		mexPrintf("Generating crystal...\n");
		crystal.CrystalGenerator(Vrf,Vend);
		crystal.SetSteadyStateTemperature(Temperature);

		// running sim...
		DynamicTemperatureLeFrogintegrator(crystal, TimeSteps, Vrf, Vend);


		mexPrintf("Copying data to matlab...\n");
		const mwSize ImDims[]={HistNx,HistNy,HistNz};
		const mwSize Tempdims[]={1,TimeSteps};

		double *HistOut;
		double *Ttotalout;
		double *Trmszout;
		double *Tsecularout;
		double *Tsecularradout;
		double *Tsecularzout;

		plhs[0] = mxCreateNumericArray(3,ImDims,mxDOUBLE_CLASS,mxREAL);
		HistOut = mxGetPr(plhs[0]);
		plhs[1] = mxCreateNumericArray(2,Tempdims,mxDOUBLE_CLASS,mxREAL);
		Ttotalout = mxGetPr(plhs[1]);
		plhs[2] = mxCreateNumericArray(2,Tempdims,mxDOUBLE_CLASS,mxREAL);
		Trmszout = mxGetPr(plhs[2]);
		plhs[3] = mxCreateNumericArray(2,Tempdims,mxDOUBLE_CLASS,mxREAL);
		Tsecularout = mxGetPr(plhs[3]);
		plhs[4] = mxCreateNumericArray(2,Tempdims,mxDOUBLE_CLASS,mxREAL);
		Tsecularradout = mxGetPr(plhs[4]);
		plhs[5] = mxCreateNumericArray(2,Tempdims,mxDOUBLE_CLASS,mxREAL);
		Tsecularzout = mxGetPr(plhs[5]);


		for(int i=0; i < HistNx; i++)
			for(int j=0; j < HistNy; j++)
				for(int k=0; k < HistNz; k++)
					HistOut[k * HistNx*HistNy + j*HistNx + i] = ((double) crystal.ReturnHist(i, j, k));

		for(int t = 0; t < TimeSteps; t++)
		{
			//cout << "Total temperature " << crystal.GetTtotal(t) << '\n';
			Ttotalout[t] = crystal.GetTrms(t);
		}

		for(int t = 0; t < TimeSteps; t++)
		{
			//cout << "RMS z temperature " << crystal.GetTtotal(t) << '\n';
			Trmszout[t] = crystal.GetTrmsz(t);
		}

		for(int t = 0; t < TimeSteps; t++)
		{
			//cout << "Secular temperature " << crystal.GetTsecular(t) << '\n';
			Tsecularout[t] = crystal.GetTsecular(t);
		}

		for(int t = 0; t < TimeSteps; t++)
		{
			//cout << "Secular temperature " << crystal.GetTsecular(t) << '\n';
			Tsecularradout[t] = crystal.GetTsecularrad(t);
		}

		for(int t = 0; t < TimeSteps; t++)
		{
			//cout << "Secular temperature " << crystal.GetTsecular(t) << '\n';
			Tsecularzout[t] = crystal.GetTsecularz(t);
		}

		// cleaning up
		crystal.FreeTemperatureArrays();
		//mxDestroyArray(plhs[1]);
		//mxDestroyArray(plhs[2]);
	}


	/*

	// checking correct number of inputs
	if (nrhs == 4)
	{
		int N = mxGetScalar(prhs[0]);
		double Vrf = mxGetScalar(prhs[1]);
		double Vend = mxGetScalar(prhs[2]);
		int TimeSteps = mxGetScalar(prhs[3]);

		mexPrintf("Allocating memory...\n");
		Ensemble crystal(N, TimeSteps);
		crystal.SetSteadyStateTemperature(0.05);

		mexPrintf("Generating crystal...\n");
		//crystal.CrystalGenerator();
		crystal.CrystalGenerator(Vrf,Vend);

		mexPrintf("Simulating...\n");
		LeFrogintegrator(crystal, TimeSteps, Vrf, Vend);
		//crystal.SavePositionToFile(TimeSteps);


		//now write some code that copies position and velocities to mxArray plhs
		if (nlhs == 2)
		{
			mexPrintf("Copying data to matlab...\n");
			const mwSize dims[]={3,N,TimeSteps};
			double *PosOut;
			double *VelOut;
			plhs[0] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
			PosOut = mxGetPr(plhs[0]);
			plhs[1] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
			VelOut = mxGetPr(plhs[1]);

			for(int dim = 0; dim < 3; dim++)
				for(int t = 0; t < TimeSteps; t++)
					for(int n = 0; n < N;++n)
					{
						// array[Z * dimX*dimY + Y*dimX + X] stupid 1d array conversion
						PosOut[t*3*N + n*3 + dim] = crystal.Position(dim, n, t); // Position and Velocity not finished
						VelOut[t*3*N + n*3 + dim] = crystal.Velocity(dim, n, t);
					}
		}
		else if(nlhs == 4)
		{
			mexPrintf("Copying data to matlab...\n");
			const mwSize dims[]={3,N,TimeSteps};
			const mwSize Tempdims[]={1,TimeSteps};
			double *PosOut;
			double *VelOut;
			double *Tout;
			double *Tsecout;

			plhs[0] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
			PosOut = mxGetPr(plhs[0]);
			plhs[1] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
			VelOut = mxGetPr(plhs[1]);
			plhs[2] = mxCreateNumericArray(2,Tempdims,mxDOUBLE_CLASS,mxREAL);
			Tout = mxGetPr(plhs[2]);
			plhs[3] = mxCreateNumericArray(2,Tempdims,mxDOUBLE_CLASS,mxREAL);
			Tsecout = mxGetPr(plhs[3]);

			for(int dim = 0; dim < 3; dim++)
				for(int t = 0; t < TimeSteps; t++)
					for(int n = 0; n < N;++n)
					{
						// array[Z * dimX*dimY + Y*dimX + X] stupid 1d array conversion
						PosOut[t*3*N + n*3 + dim] = crystal.Position(dim, n, t); // Position and Velocity not finished
						VelOut[t*3*N + n*3 + dim] = crystal.Velocity(dim, n, t);
					}

			for(int t = 0; t < TimeSteps; t++)
			{
				Tout[t] = crystal.Ttot(t);
				Tsecout[t] = crystal.Tsec(t);
			}

		}
	}
	else if (nrhs == 5)
	{
		int N = mxGetScalar(prhs[0]);
		double Vrf = mxGetScalar(prhs[1]);
		double Vend = mxGetScalar(prhs[2]);
		int TimeSteps = mxGetScalar(prhs[3]);
		double Temperature = mxGetScalar(prhs[4]);

		mexPrintf("Allocating memory...\n");
		FastEnsemble crystal(N);


		mexPrintf("Generating crystal...\n");
		//crystal.CrystalGenerator();
		crystal.CrystalGeneratorForPseudoPotential(Vrf,Vend);
		crystal.SetSteadyStateTemperature(Temperature);

		mexPrintf("Simulating...\n");
		time_t start,end;
		double dif;
		time (&start);
		//LeFrogintegrator(crystal, TimeSteps, Vrf, Vend);
		//CudaLeFrogintegrator(crystal, TimeSteps, Vrf, Vend);
		TestCudaLeFrogintegrator(crystal, TimeSteps, Vrf, Vend);
		//GridIntegrator(crystal, TimeSteps, Vrf, Vend);
		//crystal.SavePositionToFile(TimeSteps);
		//LeFrogMGintegrator(crystal, TimeSteps, Vrf, Vend);
		//GridIntegrator(crystal, TimeSteps, Vrf, Vend);
		time (&end);
		dif = difftime (end,start);
		cout << "Simulation took " << dif << " second\n";

		//mexPrintf("Saving positions to VMD file...\n");
		//crystal.SavePositionToFile();



		//now write some code that copies position and velocities to mxArray plhs
		if (nlhs == 2)
		{
			mexPrintf("Copying data to matlab...\n");
			const mwSize dims[]={3,N};
			double *PosOut;
			double *VelOut;
			plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
			PosOut = mxGetPr(plhs[0]);
			plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
			VelOut = mxGetPr(plhs[1]);

			for(int dim = 0; dim < 3; dim++)
				for(int n = 0; n < N;++n)
				{
					// array[Z * dimX*dimY + Y*dimX + X] stupid 1d array conversion
					PosOut[n*3 + dim] = crystal.Position(dim, n); // Position and Velocity not finished
					VelOut[n*3 + dim] = crystal.Velocity(dim, n);
				}

			// cleaning up
			mxDestroyArray(plhs[0]);
			mxDestroyArray(plhs[1]);



		}
		else if(nlhs == 3)
		{
			mexPrintf("Copying data to matlab...\n");
			const mwSize dims[]={3,N};
			const mwSize ImDims[]={HistNx,HistNy,HistNz};

			double *PosOut;
			double *VelOut;
			double *HistOut;
			plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
			PosOut = mxGetPr(plhs[0]);
			plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
			VelOut = mxGetPr(plhs[1]);
			plhs[2] = mxCreateNumericArray(3,ImDims,mxDOUBLE_CLASS,mxREAL);
			HistOut = mxGetPr(plhs[2]);

			for(int dim = 0; dim < 3; dim++)
				for(int n = 0; n < N;++n)
				{
					// array[Z * dimX*dimY + Y*dimX + X] stupid 1d array conversion
					PosOut[n*3 + dim] = crystal.Position(dim, n); // Position and Velocity not finished
					VelOut[n*3 + dim] = crystal.Velocity(dim, n);
				}

			for(int i=0; i < HistNx; i++)
				for(int j=0; j < HistNy; j++)
					for(int k=0; k < HistNz; k++)
						HistOut[k * HistNx*HistNy + j*HistNx + i] = ((double) crystal.ReturnHist(i, j, k));

			// cleaning up
			mxDestroyArray(plhs[0]);
			mxDestroyArray(plhs[1]);
			// mxDestroyArray(plhs[2]);

		}
		else if(nlhs == 4)
		{
			mexPrintf("Copying data to matlab...\n");
			const mwSize dims[]={3,N};
			const mwSize ImDims[]={HistNx,HistNy,HistNz};
			const mwSize VelDims[]={HistNx,HistNy,HistNz};

			double *PosOut;
			double *VelOut;
			double *HistOut;
			double *VelHistOut;
			plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
			PosOut = mxGetPr(plhs[0]);
			plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
			VelOut = mxGetPr(plhs[1]);
			plhs[2] = mxCreateNumericArray(3,ImDims,mxDOUBLE_CLASS,mxREAL);
			HistOut = mxGetPr(plhs[2]);
			plhs[3] = mxCreateNumericArray(3,VelDims,mxDOUBLE_CLASS,mxREAL);
			VelHistOut = mxGetPr(plhs[3]);

			for(int dim = 0; dim < 3; dim++)
				for(int n = 0; n < N;++n)
				{
					// array[Z * dimX*dimY + Y*dimX + X] stupid 1d array conversion
					PosOut[n*3 + dim] = crystal.Position(dim, n); // Position and Velocity not finished
					VelOut[n*3 + dim] = crystal.Velocity(dim, n);
				}

			for(int i=0; i < HistNx; i++)
				for(int j=0; j < HistNy; j++)
					for(int k=0; k < HistNz; k++)
					{
						HistOut[k * HistNx*HistNy + j*HistNx + i] = ((double) crystal.ReturnHist(i, j, k));
						VelHistOut[k * HistNx*HistNy + j*HistNx + i] = ((double) crystal.ReturnVelHist(i, j, k));
					}

			// cleaning up
			mxDestroyArray(plhs[0]);
			mxDestroyArray(plhs[1]);
			//	mxDestroyArray(plhs[2]);
			//	mxDestroyArray(plhs[3]);



		}
		else if(nlhs == 6)
		{
			mexPrintf("Copying data to matlab...\n");
			const mwSize dims[]={3,N};
			const mwSize ImDims[]={HistNx,HistNy,HistNz};
			const mwSize VelDims[]={HistNx,HistNy,HistNz};
			const mwSize ScalarDims[] = {1};

			double *PosOut;
			double *VelOut;
			double *HistOut;
			double *VelHistOut;
			double *AvgTempOut;
			double *StdTempOut;

			plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
			PosOut = mxGetPr(plhs[0]);
			plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
			VelOut = mxGetPr(plhs[1]);
			plhs[2] = mxCreateNumericArray(3,ImDims,mxDOUBLE_CLASS,mxREAL);
			HistOut = mxGetPr(plhs[2]);
			plhs[3] = mxCreateNumericArray(3,VelDims,mxDOUBLE_CLASS,mxREAL);
			VelHistOut = mxGetPr(plhs[3]);
			plhs[4] = mxCreateNumericArray(1,ScalarDims,mxDOUBLE_CLASS,mxREAL);
			AvgTempOut = mxGetPr(plhs[4]);
			plhs[5] = mxCreateNumericArray(1,ScalarDims,mxDOUBLE_CLASS,mxREAL);
			StdTempOut = mxGetPr(plhs[5]);



			for(int dim = 0; dim < 3; dim++)
				for(int n = 0; n < N;++n)
				{
					// array[Z * dimX*dimY + Y*dimX + X] stupid 1d array conversion
					PosOut[n*3 + dim] = crystal.Position(dim, n); // Position and Velocity not finished
					VelOut[n*3 + dim] = crystal.Velocity(dim, n);
				}

			for(int i=0; i < HistNx; i++)
				for(int j=0; j < HistNy; j++)
					for(int k=0; k < HistNz; k++)
					{
						HistOut[k * HistNx*HistNy + j*HistNx + i] = ((double) crystal.ReturnHist(i, j, k));
						VelHistOut[k * HistNx*HistNy + j*HistNx + i] = ((double) crystal.ReturnVelHist(i, j, k));
					}

			AvgTempOut[0] = crystal.GetActualTemperature();
			StdTempOut[0] = crystal.GetActualTemperatureSTD();

			// cleaning up
			mxDestroyArray(plhs[0]);
			mxDestroyArray(plhs[1]);
			//	mxDestroyArray(plhs[2]);
			//	mxDestroyArray(plhs[3]);



		}
		// cleaning up
		cout << "cleaning up...\n";
		crystal.CleanUpEnsemble();
	}
	 */

}
