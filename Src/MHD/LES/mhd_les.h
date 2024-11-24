/* ---------------------------------------------------------------------
   LES module header file
   --------------------------------------------------------------------- */
#pragma once

#include "pluto.h"

#include "les_phi.h"

double get_MHD_LES_Ds();
void set_MHD_LES_Ds(double ds);

double LES_PHI_MHD(const Data *d, int i, int j, int k);

void LES_ViscousFlux_MHD(const Data *d, double **ViF, double **ViS,
                         double *dcoeff, int beg, int end, Grid *grid);

void LES_ViscousRHS_MHD(const Data *d, Data_Arr dU, double *dcoeff,
                        double **aflux, double dt, int beg, int end,
                        Grid *grid);
