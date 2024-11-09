/* ---------------------------------------------------------------------
   LES module header file
   --------------------------------------------------------------------- */

#include "pluto.h"

#ifndef LES_FMIN
#define LES_FMIN 0.5 /* blah blah  */
#endif

void LES_ViscousFlux(const Data *d, Data_Arr flux, Grid *grid);

void LES_ViscousFluxNew(const Data *d, double **ViF, double **ViS,
                        double *dcoeff, int beg, int end, Grid *grid);

void LES_ViscousRHS(const Data *d, Data_Arr dU, double *dcoeff, double **aflux,
                    double dt, int beg, int end, Grid *grid);

double ***LES_GetFilter();
