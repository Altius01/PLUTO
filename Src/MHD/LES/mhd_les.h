/* ---------------------------------------------------------------------
   LES module header file
   --------------------------------------------------------------------- */

#include "pluto.h"

void LES_ViscousFlux_MHD(const Data *d, double **ViF, double **ViS,
                         double *dcoeff, int beg, int end, Grid *grid);

void LES_ViscousRHS_MHD(const Data *d, Data_Arr dU, double *dcoeff,
                        double **aflux, double dt, int beg, int end,
                        Grid *grid);
