/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Lax-Friedrechs (Rusanov) Riemann solver for RHD.

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    Dec 03, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"

/* ********************************************************************* */
void LF_Solver (const Sweep *sweep, int beg, int end, 
                double *cmax, Grid *grid)
/*!
 * Solve Riemann problem for the relativistic Hydro (RHD) equations
 * using the Lax-Friedrichs (Rusanov) Riemann solver.
 * 
 * \param[in,out] sweep   pointer to Sweep structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */

{
  int    nv, i;
  const State   *stateL = &(sweep->stateL);
  const State   *stateR = &(sweep->stateR);
  double *uL, *uR, cL, cR;
  double **fL = stateL->flux, **fR = stateR->flux;
  double  *pL = stateL->prs,   *pR = stateR->prs;
  static double *cminL, *cmaxL;
  static double *cminR, *cmaxR;

  if (cminL == NULL){
    cminL = ARRAY_1D(NMAX_POINT, double);
    cmaxL = ARRAY_1D(NMAX_POINT, double);
    cminR = ARRAY_1D(NMAX_POINT, double);
    cmaxR = ARRAY_1D(NMAX_POINT, double);
  }

#if TIME_STEPPING == CHARACTERISTIC_TRACING
{
  static int first_call = 1;
  if (first_call){
    print ("! LF_Solver(): employment of this solver with ");
    print ("CHARACTERISTIC_TRACING may degrade order of accuracy to 1.\n");
    first_call = 0;
  }
}
#endif
  
/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (stateL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (stateR, beg, end, FACE_CENTER, grid);

  Flux (stateL, beg, end);
  Flux (stateR, beg, end);

  MaxSignalSpeed (stateL, cminL, cmaxL, beg, end);
  MaxSignalSpeed (stateR, cminR, cmaxR, beg, end);

  for (i = beg; i <= end; i++) {
    cL = MAX(fabs(cminL[i]), fabs(cmaxL[i]));
    cR = MAX(fabs(cminR[i]), fabs(cmaxR[i]));
    cmax[i] = MAX(cL, cR);

    uL = stateL->u[i];
    uR = stateR->u[i];
    NFLX_LOOP(nv) {
      sweep->flux[i][nv] = 0.5*(fL[i][nv] + fR[i][nv] - cmax[i]*(uR[nv] - uL[nv]));
    }
    sweep->press[i] = 0.5*(pL[i] + pR[i]);
  }
}
