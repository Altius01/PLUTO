#include "pluto.h"

#ifndef LES_FMIN
#define LES_FMIN 0.5 /* blah blah  */
#endif

static double ***f;

/* ********************************************************************* */
void ViscousFlux(const Data *d, double **ViF, double **ViS, double *dcoeff,
                 int beg, int end, Grid *grid) {
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
#define D_DX_I(q) (q[k][j][i + 1] - q[k][j][i])
#define D_DY_J(q) (q[k][j + 1][i] - q[k][j][i])
#define D_DZ_K(q) (q[k + 1][j][i] - q[k][j][i])

#define D_DY_I(q)                                                              \
  (0.25 * (q[k][j + 1][i] + q[k][j + 1][i + 1]) -                              \
   0.25 * (q[k][j - 1][i] + q[k][j - 1][i + 1]))

#define D_DZ_I(q)                                                              \
  (0.25 * (q[k + 1][j][i] + q[k + 1][j][i + 1]) -                              \
   0.25 * (q[k - 1][j][i] + q[k - 1][j][i + 1]))

#define D_DX_J(q)                                                              \
  (0.25 * (q[k][j][i + 1] + q[k][j + 1][i + 1]) -                              \
   0.25 * (q[k][j][i - 1] + q[k][j + 1][i - 1]))

#define D_DZ_J(q)                                                              \
  (0.25 * (q[k + 1][j][i] + q[k + 1][j + 1][i]) -                              \
   0.25 * (q[k - 1][j][i] + q[k - 1][j + 1][i]))

#define D_DX_K(q)                                                              \
  (0.25 * (q[k][j][i + 1] + q[k + 1][j][i + 1]) -                              \
   0.25 * (q[k][j][i - 1] + q[k + 1][j][i - 1]))

#define D_DY_K(q)                                                              \
  (0.25 * (q[k][j + 1][i] + q[k + 1][j + 1][i]) -                              \
   0.25 * (q[k][j - 1][i] + q[k + 1][j - 1][i]))

  {
    int i, j, k, n, nv;
    double nu1, nu2;
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
    double div, nu_t, rhom;
    double wx, wy, wz, qx, qy, qz;
    double fmx, fmy, fmz, fe;
    double dtdx, dtdy, dtdz;
    double nu_max;
    double Cs = 0.1, Csdl2, dH;
    double Prandtl = 0.71;

    static double *tau_xx, *tau_xy, *tau_xz, *tau_yx, *tau_yy, *tau_yz, *tau_zx,
        *tau_zy, *tau_zz;
    static double **Ux_ave, **Uy_ave, **Uz_ave;
    static double **dUx_ave_dx, **dUx_ave_dz;
    static double **dUy_ave_dx, **dUy_ave_dz;
    static double **dUz_ave_dx, **dUz_ave_dz;
    static double ***H;

    double ***Vx = d->Vc[VX1];
    double ***Vy = d->Vc[VX2];
    double ***Vz = d->Vc[VX3];
    double ***vx, ***vy, ***vz, ***rho, ***pr;
    double *r, *th, r_1, dr, s_1, tan_1;
    double vc[NVAR], vi[NVAR]; /* Center and interface values */
    static double *one_dVr,
        *one_dmu; /*auxillary volume components for r_1 singularity @
                     cylindrical and spherical*/

#if INCLUDE_LES == NO
    return; /* do nothing, just return */
#endif

    /* --------------------------------------------------------
       0. Set pointers to coordinates and grid indices,
          allocate memory
       -------------------------------------------------------- */

    if (tau_xx == NULL) {
      tau_xx = ARRAY_1D(NMAX_POINT, double);
      tau_xy = ARRAY_1D(NMAX_POINT, double);
      tau_xz = ARRAY_1D(NMAX_POINT, double);
      tau_yx = ARRAY_1D(NMAX_POINT, double);
      tau_yy = ARRAY_1D(NMAX_POINT, double);
      tau_yz = ARRAY_1D(NMAX_POINT, double);
      tau_zx = ARRAY_1D(NMAX_POINT, double);
      tau_zy = ARRAY_1D(NMAX_POINT, double);
      tau_zz = ARRAY_1D(NMAX_POINT, double);
    }

    i = g_i;
    j = g_j;
    k = g_k;

    dxVx = dxVy = dxVz = 0.0;
    dyVx = dyVy = dyVz = 0.0;
    dzVx = dzVy = dzVz = 0.0;

    if (f == NULL) {
      f = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
      Ux_ave = ARRAY_2D(NX3_TOT, NX1_TOT, double);
      Uy_ave = ARRAY_2D(NX3_TOT, NX1_TOT, double);
      Uz_ave = ARRAY_2D(NX3_TOT, NX1_TOT, double);

      dUx_ave_dx = ARRAY_2D(NX3_TOT, NX1_TOT, double);
      dUx_ave_dz = ARRAY_2D(NX3_TOT, NX1_TOT, double);

      dUy_ave_dx = ARRAY_2D(NX3_TOT, NX1_TOT, double);
      dUy_ave_dz = ARRAY_2D(NX3_TOT, NX1_TOT, double);

      dUz_ave_dx = ARRAY_2D(NX3_TOT, NX1_TOT, double);
      dUz_ave_dz = ARRAY_2D(NX3_TOT, NX1_TOT, double);

      H = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    }

    dx = grid->dx[IDIR][IBEG];
    dy = grid->dx[JDIR][JBEG];
    dz = grid->dx[KDIR][KBEG];

    dtdx = g_dt / dx;
    dtdy = g_dt / dy;
    dtdz = g_dt / dz;

    rho = d->Vc[RHO];
    vx = d->Vc[VX1];
    vy = d->Vc[VX2];
    vz = d->Vc[VX3];
    pr = d->Vc[PRS];

    Csdl2 = Cs * Cs * pow(dx * dy * dz, 2.0 / 3.0);

    nu_max = 1.e-12;

    /* ----------------------------------------------
         1. Compute average velocity along y
      .     [remember to run with the -no-x2par switch]
         ---------------------------------------------- */

    for (k = KBEG - 2; k <= KEND + 2; k++) {
      for (i = IBEG - 2; i <= IEND + 2; i++) {

        Ux_ave[k][i] = Uy_ave[k][i] = Uz_ave[k][i] = 0.0;
        for (j = JBEG; j <= JEND; j++) {
          Ux_ave[k][i] += vx[k][j][i];
          Uy_ave[k][i] += vy[k][j][i];
          Uz_ave[k][i] += vz[k][j][i];
        }
        Ux_ave[k][i] /= (double)NX2;
        Uy_ave[k][i] /= (double)NX2;
        Uz_ave[k][i] /= (double)NX2;
      }
    }

    /* ----------------------------------------------
       2. Compute average <x> and <z> derivatives
       ---------------------------------------------- */

    for (k = KBEG - 1; k <= KEND + 1; k++) {
      for (i = IBEG - 1; i <= IEND + 1; i++) {
        dUx_ave_dx[k][i] = 0.5 * (Ux_ave[k][i + 1] - Ux_ave[k][i - 1]) / dx;
        dUx_ave_dz[k][i] = 0.5 * (Ux_ave[k + 1][i] - Ux_ave[k - 1][i]) / dz;

        dUy_ave_dx[k][i] = 0.5 * (Uy_ave[k][i + 1] - Uy_ave[k][i - 1]) / dx;
        dUy_ave_dz[k][i] = 0.5 * (Uy_ave[k + 1][i] - Uy_ave[k - 1][i]) / dz;

        dUz_ave_dx[k][i] = 0.5 * (Uz_ave[k][i + 1] - Uz_ave[k][i - 1]) / dx;
        dUz_ave_dz[k][i] = 0.5 * (Uz_ave[k + 1][i] - Uz_ave[k - 1][i]) / dz;
      }
    }

    /* ----------------------------------------------
       3. Compute "f" function
       ---------------------------------------------- */

    for (k = KBEG - 1; k <= KEND + 1; k++) {
      for (j = JBEG - 1; j <= JEND + 1; j++) {
        for (i = IBEG - 1; i <= IEND + 1; i++) {

          dxVx =
              0.5 * (vx[k][j][i + 1] - vx[k][j][i - 1]) / dx - dUx_ave_dx[k][i];
          dxVy =
              0.5 * (vy[k][j][i + 1] - vy[k][j][i - 1]) / dx - dUy_ave_dx[k][i];
          dxVz =
              0.5 * (vz[k][j][i + 1] - vz[k][j][i - 1]) / dx - dUz_ave_dx[k][i];

          dyVx = 0.5 * (vx[k][j + 1][i] - vx[k][j - 1][i]) / dy;
          dyVy = 0.5 * (vy[k][j + 1][i] - vy[k][j - 1][i]) / dy;
          dyVz = 0.5 * (vz[k][j + 1][i] - vz[k][j - 1][i]) / dy;

          dzVx =
              0.5 * (vx[k + 1][j][i] - vx[k - 1][j][i]) / dz - dUx_ave_dz[k][i];
          dxVy =
              0.5 * (vy[k + 1][j][i] - vy[k - 1][j][i]) / dz - dUy_ave_dz[k][i];
          dzVz =
              0.5 * (vz[k + 1][j][i] - vz[k - 1][j][i]) / dz - dUz_ave_dz[k][i];

          wx = dyVz - dxVy;
          wy = dzVx - dxVz;
          wz = dxVy - dyVx;

          qx = wx * dxVx + wy * dyVx + wz * dzVx;
          qy = wx * dxVy + wy * dyVy + wz * dxVy;
          qz = wx * dxVz + wy * dyVz + wz * dzVz;

          f[k][j][i] = sqrt(qx * qx + qy * qy + qz * qz);
          f[k][j][i] /= (wx * wx + wy * wy + wz * wz) * LES_FMIN + 1.e-6;
          /*
              if (f[k][j][i] < 1.e-3) f[k][j][i] = -1.0;
              else f[k][j][i] -= (wx*wx + wy*wy + wz*wz)*g_inputParam[FMIN];
          */
          /*
          printf ("%d %d %d   N = %12.6e, D = %12.6e\n",
               i,j,k, f[k][j][i],(wx*wx + wy*wy + wz*wz) );
          */

          // for (nv = 0; nv < NVAR; nv++)
          //   rhs[nv][k][j][i] = 0.0;
        }
      }
    }

    /* ----------------------------------------------
       5. Compute Specific enthalpy
       ---------------------------------------------- */

    KTOT_LOOP(k) {
      JTOT_LOOP(j) {
        ITOT_LOOP(i) {
          S = vx[k][j][i] * vx[k][j][i] + vy[k][j][i] * vy[k][j][i] +
              vz[k][j][i] * vz[k][j][i];

          H[k][j][i] =
              pr[k][j][i] * g_gamma / (g_gamma - 1.0) / rho[k][j][i] + 0.5 * S;
        }
      }
    }

    /* ----------------------------------------------
       6. Compute Viscous Fluxes
       ---------------------------------------------- */

    if (g_dir == IDIR) {

      /* ------------------------
         6.1  X-Sweep
         ------------------------ */

      for (k = KBEG; k <= KEND; k++) {
        for (j = JBEG; j <= JEND; j++) {
          for (i = IBEG - 1; i <= IEND; i++) {

            rhom = 0.5 * (rho[k][j][i] + rho[k][j][i + 1]);

            dxVx = (D_DX_I(vx)) / dx;
            dxVy = (D_DX_I(vy)) / dx;
            dxVz = (D_DX_I(vz)) / dx;

            dyVx = (D_DY_I(vx)) / dy;
            dyVy = (D_DY_I(vy)) / dy;
            dyVz = (D_DY_I(vz)) / dy;

            dzVx = (D_DZ_I(vx)) / dz;
            dzVy = (D_DZ_I(vy)) / dz;
            dzVz = (D_DZ_I(vz)) / dz;

            Sxx = 0.5 * (dxVx + dxVx);
            Sxy = 0.5 * (dyVx + dxVy);
            Sxz = 0.5 * (dzVx + dxVz);

            Syy = 0.5 * (dyVy + dyVy);
            Syz = 0.5 * (dzVy + dyVz);
            Szz = 0.5 * (dzVz + dzVz);

            div = dxVx + dyVy + dzVz;

            S = Sxx * Sxx + Syy * Syy + Szz * Szz +
                2.0 * (Sxy * Sxy + Sxz * Sxz + Syz * Syz);
            S = sqrt(2.0 * S);
            nu_t = Csdl2 * S;
            nu_max = MAX(nu_max, nu_t);

            scrh = dtdx * 2.0 * rhom * nu_t;

            fmx = scrh * (Sxx - div / 3.0);
            fmy = scrh * Sxy;
            fmz = scrh * Sxz;

            dH = dtdx * (H[k][j][i + 1] - H[k][j][i]);
            fe = rhom * nu_t / Prandtl * dH / dx;

            ViF[i][MX1] = fmx;
            ViF[i][MX2] = fmy;
            ViF[i][MX3] = fmz;

#if HAVE_ENERGY
            ViF[i][ENG] = fe;
#endif
          }
        }
      }

    } else if (g_dir == JDIR) {

      /* ------------------------
         6.2  Y-Sweep
         ------------------------ */

      for (k = KBEG; k <= KEND; k++) {
        for (j = JBEG - 1; j <= JEND; j++) {
          for (i = IBEG; i <= IEND; i++) {

            rhom = 0.5 * (rho[k][j][i] + rho[k][j + 1][i]);

            dxVx = (D_DX_J(vx)) / dx;
            dxVy = (D_DX_J(vy)) / dx;
            dxVz = (D_DX_J(vz)) / dx;

            dyVx = (D_DY_J(vx)) / dy;
            dyVy = (D_DY_J(vy)) / dy;
            dyVz = (D_DY_J(vz)) / dy;

            dzVx = (D_DZ_J(vx)) / dz;
            dzVy = (D_DZ_J(vy)) / dz;
            dzVz = (D_DZ_J(vz)) / dz;

            Sxx = 0.5 * (dxVx + dxVx);
            Sxy = 0.5 * (dyVx + dxVy);
            Sxz = 0.5 * (dzVx + dxVz);

            Syy = 0.5 * (dyVy + dyVy);
            Syz = 0.5 * (dzVy + dyVz);
            Szz = 0.5 * (dzVz + dzVz);

            div = dxVx + dyVy + dzVz;

            S = Sxx * Sxx + Syy * Syy + Szz * Szz +
                2.0 * (Sxy * Sxy + Sxz * Sxz + Syz * Syz);
            S = sqrt(2.0 * S);
            nu_t = Csdl2 * S;
            nu_max = MAX(nu_max, nu_t);

            scrh = dtdy * 2.0 * rhom * nu_t;

            fmx = scrh * Sxy;
            fmy = scrh * (Syy - div / 3.0);
            fmz = scrh * Syz;

            dH = dtdy * (H[k][j + 1][i] - H[k][j][i]);
            fe = rhom * nu_t / Prandtl * dH / dy;

            ViF[j][MX1] = fmx;
            ViF[j][MX2] = fmy;
            ViF[j][MX3] = fmz;

#if HAVE_ENERGY
            ViF[j][ENG] = fe;
#endif
          }
        }
      }

    } else if (g_dir == KDIR) {

      /* ------------------------
         6.3 Z-Sweep
         ------------------------ */

      for (k = KBEG - 1; k <= KEND; k++) {
        for (j = JBEG; j <= JEND; j++) {
          for (i = IBEG; i <= IEND; i++) {

            rhom = 0.5 * (rho[k + 1][j][i] + rho[k][j][i]);

            dxVx = (D_DX_K(vx)) / dx;
            dxVy = (D_DX_K(vy)) / dx;
            dxVz = (D_DX_K(vz)) / dx;

            dyVx = (D_DY_K(vx)) / dy;
            dyVy = (D_DY_K(vy)) / dy;
            dyVz = (D_DY_K(vz)) / dy;

            dzVx = (D_DZ_K(vx)) / dz;
            dzVy = (D_DZ_K(vy)) / dz;
            dzVz = (D_DZ_K(vz)) / dz;

            Sxx = 0.5 * (dxVx + dxVx);
            Sxy = 0.5 * (dyVx + dxVy);
            Sxz = 0.5 * (dzVx + dxVz);

            Syy = 0.5 * (dyVy + dyVy);
            Syz = 0.5 * (dzVy + dyVz);
            Szz = 0.5 * (dzVz + dzVz);

            div = dxVx + dyVy + dzVz;

            S = Sxx * Sxx + Syy * Syy + Szz * Szz +
                2.0 * (Sxy * Sxy + Sxz * Sxz + Syz * Syz);
            S = sqrt(2.0 * S);
            nu_t = Csdl2 * S;
            nu_max = MAX(nu_max, nu_t);

            scrh = dtdz * 2.0 * rhom * nu_t;

            fmx = scrh * Sxz;
            fmy = scrh * Syz;
            fmz = scrh * (Szz - div / 3.0);

            dH = dtdz * (H[k + 1][j][i] - H[k][j][i]);
            fe = rhom * nu_t / Prandtl * dH / dz;

            ViF[k][MX1] = fmx;
            ViF[k][MX2] = fmy;
            ViF[k][MX3] = fmz;

#if HAVE_ENERGY
            ViF[k][ENG] = fe;
#endif
          }
        }
      }
    }
  }
}

#undef D_DX_I
#undef D_DY_I
#undef D_DZ_I
#undef D_DX_J
#undef D_DY_J
#undef D_DZ_J
#undef D_DX_K
#undef D_DY_K
#undef D_DZ_K

double ***LES_GetFilter() { return (f); }
