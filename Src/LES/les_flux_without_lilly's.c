#include "pluto.h"

#ifndef LES_FMIN
#define LES_FMIN 0.5 /* blah blah  */
#endif

/* ********************************************************************* */
void LES_ViscousFluxNew(const Data *d, double **ViF, double **ViS,
                        double *dcoeff, int beg, int end, Grid *grid)
/*!
 *
 *  \param [in]      V  data array containing cell-centered quantities
 *  \param [in,out]  ViF pointer to viscous fluxes
 *  \param [in,out]  ViS pointer to viscous source terms
 *  \param [in,out]  dcoeff  pointer to diffusion coefficient for dt
 *calculation \param [in]      beg     integer, index for loop beg \param [in]
 *end     integer, index for loop end \param [in]      grid  pointer to array
 *of Grid structures
 *
 *  \return This function has no return value.
 ************************************************************************* */
{
  int i, j, k, n, nv;
  const double c13 = (1.0) / (3.0);
  double nu_t;
  double dx, dy, dz;
  double inv_dx1, inv_dx2, inv_dx3;
  double *x1 = grid->x[IDIR], *dx1 = grid->dx[IDIR];
  double *x2 = grid->x[JDIR], *dx2 = grid->dx[JDIR];
  double *x3 = grid->x[KDIR], *dx3 = grid->dx[KDIR];
  double *x1r = grid->xr[IDIR];
  double *x2r = grid->xr[JDIR];
  double *x3r = grid->xr[KDIR];
  double dxVx, dxVy, dxVz, dyVx, dyVy, dyVz, dzVx, dzVy, dzVz;
  double S, Sxx, Syy, Szz, Sxy, Sxz, Syz;
  double divV, divr, scrh;
  double wx, wy, wz, qx, qy, qz;
  double fmx, fmy, fmz, fe;
  // double dtdx, dtdy, dtdz;
  double nu_max;
  double Cs = 0.1, Csdl2, dH;
  double Prandtl = 0.71;

  static double ***H;

  double ***Vx = d->Vc[VX1];
  double ***Vy = d->Vc[VX2];
  double ***Vz = d->Vc[VX3];
  double ***rho, ***pr;
  double *r, *th, r_1, dr, s_1, tan_1;
  double vc[NVAR], vi[NVAR];        /* Center and interface values */
  static double *one_dVr, *one_dmu; /*auxillary volume components for r_1
                                       singularity @ cylindrical and spherical*/

  /* --------------------------------------------------------
     0. Set pointers to coordinates and grid indices,
        allocate memory
     -------------------------------------------------------- */

  i = g_i;
  j = g_j;
  k = g_k;

  dxVx = dxVy = dxVz = 0.0;
  dyVx = dyVy = dyVz = 0.0;
  dzVx = dzVy = dzVz = 0.0;

  if (H == NULL) {
    // H = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  }

  dx = grid->dx[IDIR][IBEG];
  dy = grid->dx[JDIR][JBEG];
  dz = grid->dx[KDIR][KBEG];

  rho = d->Vc[RHO];

  // pr = d->Vc[PRS];

  Csdl2 = Cs * Cs * pow(dx * dy * dz, 2.0 / 3.0);

  nu_max = 1.e-12;

  /* ----------------------------------------------
       5. Compute Specific enthalpy
       ---------------------------------------------- */

  int i_h, j_h, k_h;
  KTOT_LOOP(k_h) {
    JTOT_LOOP(j_h) {
      ITOT_LOOP(i_h) {
        S = Vx[k_h][j_h][i_h] * Vx[k_h][j_h][i_h] +
            Vy[k_h][j_h][i_h] * Vy[k_h][j_h][i_h] +
            Vz[k_h][j_h][i_h] * Vz[k_h][j_h][i_h];

        // H[k_h][j_h][i_h] =
        //     pr[k_h][j_h][i_h] * g_gamma / (g_gamma - 1.0) /
        //     rho[k_h][j_h][i_h] + 0.5 * S;
      }
    }
  }

  if (g_dir == IDIR) {

    /* ------------------------------------------------------
       1. Compute derivatives of velocity at x1-zone
          interfaces
       ------------------------------------------------------ */

    inv_dx2 = 1.0 / dx2[j];
    inv_dx3 = 1.0 / dx3[k];
    for (i = beg; i <= end; i++) {
      inv_dx1 = 1.0 / dx1[i];

      /* -- 1a. Compute face- and cell-centered values  -- */

      NVAR_LOOP(nv) {
        vi[nv] = 0.5 * (d->Vc[nv][k][j][i] + d->Vc[nv][k][j][i + 1]);
        vc[nv] = d->Vc[nv][k][j][i];
      }
      /* -- 1b. Compute viscosity and velocity derivatives -- */

#if INCLUDE_IDIR
      dxVx = FDIFF_X1(Vx, k, j, i) * inv_dx1;
      dxVy = FDIFF_X1(Vy, k, j, i) * inv_dx1;
      dxVz = FDIFF_X1(Vz, k, j, i) * inv_dx1;
#endif
#if INCLUDE_JDIR
      dyVx =
          0.5 * (CDIFF_X2(Vx, k, j, i) + CDIFF_X2(Vx, k, j, i + 1)) * inv_dx2;
      dyVy =
          0.5 * (CDIFF_X2(Vy, k, j, i) + CDIFF_X2(Vy, k, j, i + 1)) * inv_dx2;
      dyVz =
          0.5 * (CDIFF_X2(Vz, k, j, i) + CDIFF_X2(Vz, k, j, i + 1)) * inv_dx2;
#endif
#if INCLUDE_KDIR
      dzVx =
          0.5 * (CDIFF_X3(Vx, k, j, i) + CDIFF_X3(Vx, k, j, i + 1)) * inv_dx3;
      dzVy =
          0.5 * (CDIFF_X3(Vy, k, j, i) + CDIFF_X3(Vy, k, j, i + 1)) * inv_dx3;
      dzVz =
          0.5 * (CDIFF_X3(Vz, k, j, i) + CDIFF_X3(Vz, k, j, i + 1)) * inv_dx3;
#endif

      /* ------------------------------------------
         1c. Compute stress tensor components and
             geometrical source terms in different
             coordinate systems
         ------------------------------------------ */

#if GEOMETRY == CARTESIAN

      Sxx = 0.5 * (dxVx + dxVx);
      Sxy = 0.5 * (dyVx + dxVy);
      Sxz = 0.5 * (dzVx + dxVz);

      Syy = 0.5 * (dyVy + dyVy);
      Syz = 0.5 * (dzVy + dyVz);
      Szz = 0.5 * (dzVz + dzVz);

      S = Sxx * Sxx + Syy * Syy + Szz * Szz +
          2.0 * (Sxy * Sxy + Sxz * Sxz + Syz * Syz);
      S = sqrt(2.0 * S);

      divV = DIM_EXPAND(dxVx, +dyVy, +dzVz);

      nu_t = Csdl2 * S;
      nu_max = MAX(nu_max, nu_t);
      dcoeff[i] = nu_t;

      scrh = 2.0 * vi[RHO] * nu_t;

      fmx = scrh * (Sxx - divV / 3.0);
      fmy = scrh * Sxy;
      fmz = scrh * Sxz;

      // dH = (H[k][j][i + 1] - H[k][j][i]);
      // fe = vi[RHO] * nu_t / Prandtl * dH / dx;

      ViS[i][MX1] = 0.0;
      ViS[i][MX2] = 0.0;
      ViS[i][MX3] = 0.0;
      // ViS[i][ENG] = 0.0;
#endif /* -- end #if GEOMETRY -- */

      /* ------------------------------------------
         1d. Compute fluxes at x1-faces
         ------------------------------------------ */

      ViF[i][MX1] = fmx;
      ViF[i][MX2] = fmy;
      ViF[i][MX3] = fmz;

      // ViF[i][ENG] = 0;
#if HAVE_ENERGY
      ViF[i][ENG] = fe;
#endif
      // printf("x sweep: VIS[MX1] = %g\n", ViS[i][MX1]);
      // printf("x sweep: VIS[MX3] = %g\n", ViS[i][MX3]);
      // printf("x sweep: VIS[ENG] = %g\n", ViS[i][ENG]);
      // printf("x sweep: VIS[MX2] = %g\n", ViS[i][MX2]);
      // printf("x sweep: VIF[MX1] = %g\n", ViF[i][MX1]);
      // printf("x sweep: VIF[MX2] = %g\n", ViF[i][MX2]);
      // printf("x sweep: VIF[MX3] = %g\n", ViF[i][MX3]);
      // printf("x sweep: VIF[ENG] = %g\n", ViF[i][ENG]);
    }
  } else if (g_dir == JDIR) {

    /* ------------------------------------------------------
       2. Compute derivatives of velocity at x2-zone
          interfaces
       ------------------------------------------------------ */

    inv_dx1 = 1.0 / dx1[i];
    inv_dx3 = 1.0 / dx3[k];
    for (j = beg; j <= end; j++) {
      inv_dx2 = 1.0 / dx2[j];

      /* -- 2a. Compute face- and cell-centered values  -- */

      NVAR_LOOP(nv) {
        vi[nv] = 0.5 * (d->Vc[nv][k][j + 1][i] + d->Vc[nv][k][j][i]);
        vc[nv] = d->Vc[nv][k][j][i];
      }
      /* -- 2b. Compute viscosity and velocity derivatives -- */

#if INCLUDE_IDIR
      dxVx =
          0.5 * (CDIFF_X1(Vx, k, j, i) + CDIFF_X1(Vx, k, j + 1, i)) * inv_dx1;
      dxVy =
          0.5 * (CDIFF_X1(Vy, k, j, i) + CDIFF_X1(Vy, k, j + 1, i)) * inv_dx1;
      dxVz =
          0.5 * (CDIFF_X1(Vz, k, j, i) + CDIFF_X1(Vz, k, j + 1, i)) * inv_dx1;
#endif
#if INCLUDE_JDIR
      dyVx = FDIFF_X2(Vx, k, j, i) * inv_dx2;
      dyVy = FDIFF_X2(Vy, k, j, i) * inv_dx2;
      dyVz = FDIFF_X2(Vz, k, j, i) * inv_dx2;
#endif
#if INCLUDE_KDIR
      dzVx =
          0.5 * (CDIFF_X3(Vx, k, j, i) + CDIFF_X3(Vx, k, j + 1, i)) * inv_dx3;
      dzVy =
          0.5 * (CDIFF_X3(Vy, k, j, i) + CDIFF_X3(Vy, k, j + 1, i)) * inv_dx3;
      dzVz =
          0.5 * (CDIFF_X3(Vz, k, j, i) + CDIFF_X3(Vz, k, j + 1, i)) * inv_dx3;
#endif

      /* ------------------------------------------
         2c. Compute stress tensor components and
             geometrical source terms in different
             coordinate systems
         ------------------------------------------ */

#if GEOMETRY == CARTESIAN
      divV = DIM_EXPAND(dxVx, +dyVy, +dzVz);

      Sxx = 0.5 * (dxVx + dxVx);
      Sxy = 0.5 * (dyVx + dxVy);
      Sxz = 0.5 * (dzVx + dxVz);

      Syy = 0.5 * (dyVy + dyVy);
      Syz = 0.5 * (dzVy + dyVz);
      Szz = 0.5 * (dzVz + dzVz);

      S = Sxx * Sxx + Syy * Syy + Szz * Szz +
          2.0 * (Sxy * Sxy + Sxz * Sxz + Syz * Syz);
      S = sqrt(2.0 * S);
      nu_t = Csdl2 * S;
      nu_max = MAX(nu_max, nu_t);
      dcoeff[j] = nu_t;

      scrh = 2.0 * vi[RHO] * nu_t;

      fmx = scrh * Sxy;
      fmy = scrh * (Syy - divV / 3.0);
      fmz = scrh * Syz;

      // dH = (H[k][j + 1][i] - H[k][j][i]);
      // fe = vi[RHO] * nu_t / Prandtl * dH / dy;

      ViS[j][MX1] = 0.0;
      ViS[j][MX2] = 0.0;
      ViS[j][MX3] = 0.0;
      // ViS[j][ENG] = 0.0;
#endif /* -- GEOMETRY -- */

      /* ------------------------------------------
         2d. Compute fluxes at x2-faces
         ------------------------------------------ */

      ViF[j][MX1] = fmx;
      ViF[j][MX2] = fmy;
      ViF[j][MX3] = fmz;

      // ViF[j][ENG] = 0;
#if HAVE_ENERGY
      ViF[j][ENG] = fe;
#endif
    }

    // printf("y sweep: VIS[MX1] = %g\n", ViS[i][MX1]);
    // printf("y sweep: VIS[MX2] = %g\n", ViS[i][MX2]);
    // printf("y sweep: VIS[MX3] = %g\n", ViS[i][MX3]);
    // printf("y sweep: VIS[ENG] = %g\n", ViS[i][ENG]);
    // printf("y sweep: VIF[MX1] = %g\n", ViF[i][MX1]);
    // printf("y sweep: VIF[MX2] = %g\n", ViF[i][MX2]);
    // printf("y sweep: VIF[MX3] = %g\n", ViF[i][MX3]);
    // printf("y sweep: VIF[ENG] = %g\n", ViF[i][ENG]);

  } else if (g_dir == KDIR) {

    /* ------------------------------------------------------
       3. Compute derivatives of velocity at x3-zone
          interfaces
       ------------------------------------------------------ */

    inv_dx1 = 1.0 / dx1[i];
    inv_dx2 = 1.0 / dx2[j];
    for (k = beg; k <= end; k++) {
      inv_dx3 = 1.0 / grid->dx[KDIR][k];

      /* -- 3a. Compute face- and cell-centered values  -- */

      NVAR_LOOP(nv) {
        vi[nv] = 0.5 * (d->Vc[nv][k][j][i] + d->Vc[nv][k + 1][j][i]);
        vc[nv] = d->Vc[nv][k][j][i];
      }
      /* -- 3b. Compute viscosity and velocity derivatives -- */

#if INCLUDE_IDIR
      dxVx =
          0.5 * (CDIFF_X1(Vx, k, j, i) + CDIFF_X1(Vx, k + 1, j, i)) * inv_dx1;
      dxVy =
          0.5 * (CDIFF_X1(Vy, k, j, i) + CDIFF_X1(Vy, k + 1, j, i)) * inv_dx1;
      dxVz =
          0.5 * (CDIFF_X1(Vz, k, j, i) + CDIFF_X1(Vz, k + 1, j, i)) * inv_dx1;
#endif
#if INCLUDE_JDIR
      dyVx =
          0.5 * (CDIFF_X2(Vx, k, j, i) + CDIFF_X2(Vx, k + 1, j, i)) * inv_dx2;
      dyVy =
          0.5 * (CDIFF_X2(Vy, k, j, i) + CDIFF_X2(Vy, k + 1, j, i)) * inv_dx2;
      dyVz =
          0.5 * (CDIFF_X2(Vz, k, j, i) + CDIFF_X2(Vz, k + 1, j, i)) * inv_dx2;
#endif
#if INCLUDE_KDIR
      dzVx = FDIFF_X3(Vx, k, j, i) * inv_dx3;
      dzVy = FDIFF_X3(Vy, k, j, i) * inv_dx3;
      dzVz = FDIFF_X3(Vz, k, j, i) * inv_dx3;
#endif

      /* ------------------------------------------
         3c. Compute stress tensor components and
             geometrical source terms in different
             coordinate systems
         ------------------------------------------ */

#if GEOMETRY == CARTESIAN

      divV = dxVx + dyVy + dzVz;

      Sxx = 0.5 * (dxVx + dxVx);
      Sxy = 0.5 * (dyVx + dxVy);
      Sxz = 0.5 * (dzVx + dxVz);

      Syy = 0.5 * (dyVy + dyVy);
      Syz = 0.5 * (dzVy + dyVz);
      Szz = 0.5 * (dzVz + dzVz);

      S = Sxx * Sxx + Syy * Syy + Szz * Szz +
          2.0 * (Sxy * Sxy + Sxz * Sxz + Syz * Syz);
      S = sqrt(2.0 * S);
      nu_t = Csdl2 * S;
      nu_max = MAX(nu_max, nu_t);
      dcoeff[k] = nu_t;

      scrh = 2.0 * vi[RHO] * nu_t;

      fmx = scrh * Sxz;
      fmy = scrh * Syz;
      fmz = scrh * (Szz - divV / 3.0);

      // dH = (H[k + 1][j][i] - H[k][j][i]);
      // fe = nu_t / Prandtl * dH / dz;

      ViS[k][MX1] = 0.0;
      ViS[k][MX2] = 0.0;
      ViS[k][MX3] = 0.0;
      // ViS[k][ENG] = 0.0;
#endif /* -- end #if GEOMETRY -- */

      /* ------------------------------------------
         3d. Compute fluxes at x3-faces
         ------------------------------------------ */

      ViF[k][MX1] = fmx;
      ViF[k][MX2] = fmy;
      ViF[k][MX3] = fmz;

      // ViF[k][ENG] = 0;

#if HAVE_ENERGY
      ViF[k][ENG] = fe;
#endif
    } /*loop*/
    // printf("z sweep: VIS[MX1] = %g\n", ViS[i][MX1]);
    // printf("z sweep: VIS[MX2] = %g\n", ViS[i][MX2]);
    // printf("z sweep: VIS[MX3] = %g\n", ViS[i][MX3]);
    // printf("z sweep: VIS[ENG] = %g\n", ViS[i][ENG]);
    // printf("z sweep: VIF[MX1] = %g\n", ViF[i][MX1]);
    // printf("z sweep: VIF[MX2] = %g\n", ViF[i][MX2]);
    // printf("z sweep: VIF[MX3] = %g\n", ViF[i][MX3]);
    // printf("z sweep: VIF[ENG] = %g\n", ViF[i][ENG]);

  } /*sweep*/
}
