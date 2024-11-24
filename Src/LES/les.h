/* ---------------------------------------------------------------------
   LES module header file
   --------------------------------------------------------------------- */

#ifndef LES_FILE_H
#define LES_FILE_H

#include "les_alpha.h"

#ifndef LES_FMIN
#define LES_FMIN 0.5 /* blah blah  */
#endif

double get_LES_Cs();
void set_LES_Cs(double cs);

double get_LES_Ys();
void set_LES_Ys(double ys);

// void LES_ViscousFlux(const Data *d, Data_Arr flux, Grid *grid);

void LES_ViscousFluxNew(const Data *d, double **ViF, double **ViS,
                        double *dcoeff, int beg, int end, Grid *grid);

void LES_ViscousRHS(const Data *d, Data_Arr dU, double *dcoeff, double **aflux,
                    double dt, int beg, int end, Grid *grid);

double ***LES_GetFilter();

#endif
