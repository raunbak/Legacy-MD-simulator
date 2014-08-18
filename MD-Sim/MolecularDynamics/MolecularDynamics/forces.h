// forces.h

#ifndef FORCES_H_
#define FORCES_H_


#include "FastIon.h"
#include "FastEnsemble.h"
#include <cstring>

// Direct calculation forces

double Ftot(FastEnsemble & ensemble, int N, int TimeStep, int dim, double Vrf, double Vend);
double Fcoulumb(FastEnsemble & ensemble, int N, int dim);
double Ftrap(FastEnsemble & ensemble, int N, int TimeStep, int dim, double Vrf, double Vend);
double Ffriction(FastEnsemble & ensemble, int N, int dim);
double Fpseudo(FastEnsemble & ensemble, int N, int dim, double Vrf, double Vend);




#endif
