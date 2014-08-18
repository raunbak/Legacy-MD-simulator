//integrator.h

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include "FastIon.h"
#include "FastEnsemble.h"

#include <cstring>
#include <iostream>

void LeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend); // FastEnsemble simple integrator
void CudaLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend); // FastEnsemble simple integrator
void DynamicTemperatureLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend); // Testing temperature control as in articles

void MADSDynamicTemperatureLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend); // Testing temperature control as in articles
void OLDCudaLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend);

void TauPeriodeCudaLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend);
#endif
