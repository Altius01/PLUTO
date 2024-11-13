#include "pluto.h"

#define ITOT_LOOP_FILTERED(i) for ((i) = 0; (i) < NX1_TOT; (i) += 2)
#define JTOT_LOOP_FILTERED(j) for ((j) = 0; (j) < NX2_TOT; (j) += 2)
#define KTOT_LOOP_FILTERED(k) for ((k) = 0; (k) < NX3_TOT; (k) += 2)

/* ********************************************************************* */
void LES_ViscousFlux_MHD(const Data *d, double **ViF, double **ViS,
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

  double ***Jx1, ***Jx2, ***Jx3;

  Jx1 = d->J[IDIR];
  Jx2 = d->J[JDIR];
  Jx3 = d->J[KDIR];

  double ***Vx = d->Vc[VX1];
  double ***Vy = d->Vc[VX2];
  double ***Vz = d->Vc[VX3];

  double dxBx, dxBy, dxBz, dyBx, dyBy, dyBz, dzBx, dzBy, dzBz;
  double S;
  double J, Jxx, Jyy, Jzz, Jxy, Jxz, Jyx, Jzx, Jyz, Jzy;
  Jxx = Jyy = Jzz = 0.0;

  double divJ, divr, scrh;
  double wx, wy, wz, qx, qy, qz;
  double fBx, fBy, fBz, fe;
  double nu_max;
  double Ds = 0.1, Dsdl2, dH;
  double Prandtl = 0.71;

  static double ***Lu_xx, ***Lu_yy, ***Lu_zz, ***Lu_xy, ***Lu_xz, ***Lu_yz;
  static double ***Lb_xy, ***Lb_xz, ***Lb_yz;
  static double ***Mu_xx, ***Mu_xy, ***Mu_xz, ***Mu_yx, ***Mu_yy, ***Mu_yz,
      ***Mu_zx, ***Mu_zy, ***Mu_zz;
  static double ***Mb_xx, ***Mb_xy, ***Mb_xz, ***Mb_yx, ***Mb_yy, ***Mb_yz,
      ***Mb_zx, ***Mb_zy, ***Mb_zz;

  static double ***P, ***Q;

  double ***Bx = d->Vc[BX1];
  double ***By = d->Vc[BX2];
  double ***Bz = d->Vc[BX3];
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

  dxBx = dxBy = dxBz = 0.0;
  dyBx = dyBy = dyBz = 0.0;
  dzBx = dzBy = dzBz = 0.0;

  if (Mb_xx == NULL) {
    Lu_xx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Lu_yy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Lu_zz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Lu_xy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Lu_xz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Lu_yz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    Lb_xy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Lb_xz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Lb_yz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    Mb_xx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Mb_xy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Mb_xz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Mb_yx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Mb_yy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Mb_yz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Mb_zx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Mb_zy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Mb_zz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    Mu_xx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Mu_xy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Mu_xz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Mu_yx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Mu_yy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Mu_yz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Mu_zx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Mu_zy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Mu_zz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    P = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Q = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  }

  dx = grid->dx[IDIR][IBEG];
  dy = grid->dx[JDIR][JBEG];
  dz = grid->dx[KDIR][KBEG];

  rho = d->Vc[RHO];
  pr = d->Vc[PRS];

  Dsdl2 = Ds * Ds * pow(dx * dy * dz, 2.0 / 3.0);

  nu_max = 1.e-12;

  /* ----------------------------------------------
       5. Compute Specific enthalpy
       ---------------------------------------------- */

  int i_h, j_h, k_h;
  KTOT_LOOP_FILTERED(k_h) {
    JTOT_LOOP_FILTERED(j_h) {
      ITOT_LOOP_FILTERED(i_h) {
        double divider = 1.0;
        double rho_f, rhovx_f, rhovy_f, rhovz_f, Bx_f, By_f, Bz_f = 0;
        double BxBx_f, ByBy_f, BzBz_f, BxBy_f, BxBz_f, ByBz_f = 0;

#if INCLUDE_IDIR
        divider *= 2;
        rho_f += (rho[k_h][j_h][i_h] + rho[k_h][j_h][i_h + 1]);

        rhovx_f += (rho[k_h][j_h][i_h] * Vx[k_h][j_h][i_h] +
                    rho[k_h][j_h][i_h + 1] * Vx[k_h][j_h][i_h + 1]);
        rhovy_f += (rho[k_h][j_h][i_h] * Vy[k_h][j_h][i_h] +
                    rho[k_h][j_h][i_h + 1] * Vy[k_h][j_h][i_h + 1]);
        rhovz_f += (rho[k_h][j_h][i_h] * Vz[k_h][j_h][i_h] +
                    rho[k_h][j_h][i_h + 1] * Vz[k_h][j_h][i_h + 1]);

        Bx_f += (Bx[k_h][j_h][i_h] + Bx[k_h][j_h][i_h + 1]);
        By_f += (By[k_h][j_h][i_h] + By[k_h][j_h][i_h + 1]);
        Bz_f += (Bz[k_h][j_h][i_h] + Bz[k_h][j_h][i_h + 1]);

        BxBx_f += (Bx[k_h][j_h][i_h] * Bx[k_h][j_h][i_h] +
                   Bx[k_h][j_h][i_h + 1] * Bx[k_h][j_h][i_h + 1]);
        ByBy_f += (By[k_h][j_h][i_h] * By[k_h][j_h][i_h] +
                   By[k_h][j_h][i_h + 1] * By[k_h][j_h][i_h + 1]);
        BzBz_f += (Bz[k_h][j_h][i_h] * Bz[k_h][j_h][i_h] +
                   Bz[k_h][j_h][i_h + 1] * Bz[k_h][j_h][i_h + 1]);
        BxBy_f += (Bx[k_h][j_h][i_h] * By[k_h][j_h][i_h] +
                   Bx[k_h][j_h][i_h + 1] * By[k_h][j_h][i_h + 1]);
        BxBz_f += (Bx[k_h][j_h][i_h] * Bz[k_h][j_h][i_h] +
                   Bx[k_h][j_h][i_h + 1] * Bz[k_h][j_h][i_h + 1]);
        ByBz_f += (By[k_h][j_h][i_h] * Bz[k_h][j_h][i_h] +
                   By[k_h][j_h][i_h + 1] * Bz[k_h][j_h][i_h + 1]);
#endif
#if INCLUDE_JDIR
        divider *= 2;
        rho_f += (rho[k_h][j_h][i_h] + rho[k_h][j_h + 1][i_h]);

        rhovx_f += (rho[k_h][j_h][i_h] * Vx[k_h][j_h][i_h] +
                    rho[k_h][j_h + 1][i_h] * Vx[k_h][j_h + 1][i_h]);
        rhovy_f += (rho[k_h][j_h][i_h] * Vy[k_h][j_h][i_h] +
                    rho[k_h][j_h + 1][i_h] * Vy[k_h][j_h + 1][i_h]);
        rhovz_f += (rho[k_h][j_h][i_h] * Vz[k_h][j_h][i_h] +
                    rho[k_h][j_h + 1][i_h] * Vz[k_h][j_h + 1][i_h]);

        Bx_f += (Bx[k_h][j_h][i_h] + Bx[k_h][j_h + 1][i_h]);
        By_f += (By[k_h][j_h][i_h] + By[k_h][j_h + 1][i_h]);
#endif
#if INCLUDE_KDIR
        divider *= 2;
        rho_f += (rho[k_h][j_h][i_h] + rho[k_h + 1][j_h][i_h]);

        rhovx_f += (rho[k_h][j_h][i_h] * Vx[k_h][j_h][i_h] +
                    rho[k_h + 1][j_h][i_h] * Vx[k_h + 1][j_h][i_h]);
        rhovy_f += (rho[k_h][j_h][i_h] * Vy[k_h][j_h][i_h] +
                    rho[k_h + 1][j_h][i_h] * Vy[k_h + 1][j_h][i_h]);
        rhovz_f += (rho[k_h][j_h][i_h] * Vz[k_h][j_h][i_h] +
                    rho[k_h + 1][j_h][i_h] * Vz[k_h + 1][j_h][i_h]);

        Bx_f += (Bx[k_h][j_h][i_h] + Bx[k_h + 1][j_h][i_h]);
        By_f += (By[k_h][j_h][i_h] + By[k_h + 1][j_h][i_h]);
        Bz_f += (Bz[k_h][j_h][i_h] + Bz[k_h + 1][j_h][i_h]);

#endif

        Lu_xx[k_h][j_h][i_h] =
            rho[k_h][j_h][i_h] * Vx[k_h][j_h][i_h] * Vx[k_h][j_h][i_h];
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
      dxBx = FDIFF_X1(Bx, k, j, i) * inv_dx1;
      dxBy = FDIFF_X1(By, k, j, i) * inv_dx1;
      dxBz = FDIFF_X1(Bz, k, j, i) * inv_dx1;
#endif
#if INCLUDE_JDIR
      dyBx =
          0.5 * (CDIFF_X2(Bx, k, j, i) + CDIFF_X2(Bx, k, j, i + 1)) * inv_dx2;
      dyBy =
          0.5 * (CDIFF_X2(By, k, j, i) + CDIFF_X2(By, k, j, i + 1)) * inv_dx2;
      dyBz =
          0.5 * (CDIFF_X2(Bz, k, j, i) + CDIFF_X2(Bz, k, j, i + 1)) * inv_dx2;
#endif
#if INCLUDE_KDIR
      dzBx =
          0.5 * (CDIFF_X3(Bx, k, j, i) + CDIFF_X3(Bx, k, j, i + 1)) * inv_dx3;
      dzBy =
          0.5 * (CDIFF_X3(By, k, j, i) + CDIFF_X3(By, k, j, i + 1)) * inv_dx3;
      dzBz =
          0.5 * (CDIFF_X3(Bz, k, j, i) + CDIFF_X3(Bz, k, j, i + 1)) * inv_dx3;
#endif

      /* ------------------------------------------
         1c. Compute stress tensor components and
             geometrical source terms in different
             coordinate systems
         ------------------------------------------ */

#if GEOMETRY == CARTESIAN

      Jxy = 0.5 * (dyBx - dxBy);
      Jyx = -Jxy;
      Jxz = 0.5 * (dzBx - dxBz);
      Jzx = -Jxz;
      Jyz = 0.5 * (dzBy - dyBz);
      Jzy = -Jyz;

      J = Jx1[k][j][i] * Jx1[k][j][i] + Jx2[k][j][i] * Jx2[k][j][i] +
          Jx3[k][j][i] * Jx3[k][j][i];

      J = sqrt(J);

      nu_t = Dsdl2 * J;
      nu_max = MAX(nu_max, nu_t);
      dcoeff[i] = nu_t;

      scrh = 2.0 * nu_t;

      fBx = scrh * Jxx;
      fBy = scrh * Jxy;
      fBz = scrh * Jxz;

      // dH = (H[k][j][i + 1] - H[k][j][i]);
      // fe = vi[RHO] * nu_t / Prandtl * dH / dx;

      ViS[i][BX1] = 0.0;
      ViS[i][BX2] = 0.0;
      ViS[i][BX3] = 0.0;
      ViS[i][ENG] = 0.0;
#endif /* -- end #if GEOMETRY -- */

      /* ------------------------------------------
         1d. Compute fluxes at x1-faces
         ------------------------------------------ */

      ViF[i][BX1] = fBx;
      ViF[i][BX2] = fBy;
      ViF[i][BX3] = fBz;

      ViF[i][ENG] = 0;
#if HAVE_ENERGY
      // ViF[i][ENG] = fe;
#endif
      // printf("x sweep: VIS[BX1] = %g\n", ViS[i][BX1]);
      // printf("x sweep: VIS[BX3] = %g\n", ViS[i][BX3]);
      // printf("x sweep: VIS[BX2] = %g\n", ViS[i][BX2]);
      // printf("x sweep: VIS[ENG] = %g\n", ViS[i][ENG]);
      // printf("x sweep: VIF[BX1] = %g\n", ViF[i][BX1]);
      // printf("x sweep: VIF[BX2] = %g\n", ViF[i][BX2]);
      // printf("x sweep: VIF[BX3] = %g\n", ViF[i][BX3]);
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
      dxBx =
          0.5 * (CDIFF_X1(Bx, k, j, i) + CDIFF_X1(Bx, k, j + 1, i)) * inv_dx1;
      dxBy =
          0.5 * (CDIFF_X1(By, k, j, i) + CDIFF_X1(By, k, j + 1, i)) * inv_dx1;
      dxBz =
          0.5 * (CDIFF_X1(Bz, k, j, i) + CDIFF_X1(Bz, k, j + 1, i)) * inv_dx1;
#endif
#if INCLUDE_JDIR
      dyBx = FDIFF_X2(Bx, k, j, i) * inv_dx2;
      dyBy = FDIFF_X2(By, k, j, i) * inv_dx2;
      dyBz = FDIFF_X2(Bz, k, j, i) * inv_dx2;
#endif
#if INCLUDE_KDIR
      dzBx =
          0.5 * (CDIFF_X3(Bx, k, j, i) + CDIFF_X3(Bx, k, j + 1, i)) * inv_dx3;
      dzBy =
          0.5 * (CDIFF_X3(By, k, j, i) + CDIFF_X3(By, k, j + 1, i)) * inv_dx3;
      dzBz =
          0.5 * (CDIFF_X3(Bz, k, j, i) + CDIFF_X3(Bz, k, j + 1, i)) * inv_dx3;
#endif

      /* ------------------------------------------
         2c. Compute stress tensor components and
             geometrical source terms in different
             coordinate systems
         ------------------------------------------ */

#if GEOMETRY == CARTESIAN
      Jxy = 0.5 * (dyBx - dxBy);
      Jyx = -Jxy;
      Jxz = 0.5 * (dzBx - dxBz);
      Jzx = -Jxz;
      Jyz = 0.5 * (dzBy - dyBz);
      Jzy = -Jyz;

      J = Jx1[k][j][i] * Jx1[k][j][i] + Jx2[k][j][i] * Jx2[k][j][i] +
          Jx3[k][j][i] * Jx3[k][j][i];

      J = sqrt(J);

      nu_t = Dsdl2 * J;

      nu_max = MAX(nu_max, nu_t);
      dcoeff[j] = nu_t;

      scrh = 2.0 * nu_t;

      fBx = scrh * Jxy;
      fBy = scrh * Jyy;
      fBz = scrh * Jyz;

      // dH = (H[k][j + 1][i] - H[k][j][i]);
      // fe = vi[RHO] * nu_t / Prandtl * dH / dy;

      ViS[j][BX1] = 0.0;
      ViS[j][BX2] = 0.0;
      ViS[j][BX3] = 0.0;
      ViS[j][ENG] = 0.0;
#endif /* -- GEOMETRY -- */

      /* ------------------------------------------
         2d. Compute fluxes at x2-faces
         ------------------------------------------ */

      ViF[j][BX1] = fBx;
      ViF[j][BX2] = fBy;
      ViF[j][BX3] = fBz;

      ViF[j][ENG] = 0;
#if HAVE_ENERGY
      // ViF[j][ENG] = fe;
#endif
    }

    // printf("y sweep: VIS[BX1] = %g\n", ViS[i][BX1]);
    // printf("y sweep: VIS[BX2] = %g\n", ViS[i][BX2]);
    // printf("y sweep: VIS[BX3] = %g\n", ViS[i][BX3]);
    // printf("y sweep: VIS[ENG] = %g\n", ViS[i][ENG]);
    // printf("y sweep: VIF[BX1] = %g\n", ViF[i][BX1]);
    // printf("y sweep: VIF[BX2] = %g\n", ViF[i][BX2]);
    // printf("y sweep: VIF[BX3] = %g\n", ViF[i][BX3]);
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
      dxBx =
          0.5 * (CDIFF_X1(Bx, k, j, i) + CDIFF_X1(Bx, k + 1, j, i)) * inv_dx1;
      dxBy =
          0.5 * (CDIFF_X1(By, k, j, i) + CDIFF_X1(By, k + 1, j, i)) * inv_dx1;
      dxBz =
          0.5 * (CDIFF_X1(Bz, k, j, i) + CDIFF_X1(Bz, k + 1, j, i)) * inv_dx1;
#endif
#if INCLUDE_JDIR
      dyBx =
          0.5 * (CDIFF_X2(Bx, k, j, i) + CDIFF_X2(Bx, k + 1, j, i)) * inv_dx2;
      dyBy =
          0.5 * (CDIFF_X2(By, k, j, i) + CDIFF_X2(By, k + 1, j, i)) * inv_dx2;
      dyBz =
          0.5 * (CDIFF_X2(Bz, k, j, i) + CDIFF_X2(Bz, k + 1, j, i)) * inv_dx2;
#endif
#if INCLUDE_KDIR
      dzBx = FDIFF_X3(Bx, k, j, i) * inv_dx3;
      dzBy = FDIFF_X3(By, k, j, i) * inv_dx3;
      dzBz = FDIFF_X3(Bz, k, j, i) * inv_dx3;
#endif

      /* ------------------------------------------
         3c. Compute stress tensor components and
             geometrical source terms in different
             coordinate systems
         ------------------------------------------ */

#if GEOMETRY == CARTESIAN

      Jxy = 0.5 * (dyBx - dxBy);
      Jyx = -Jxy;
      Jxz = 0.5 * (dzBx - dxBz);
      Jzx = -Jxz;
      Jyz = 0.5 * (dzBy - dyBz);
      Jzy = -Jyz;

      J = Jx1[k][j][i] * Jx1[k][j][i] + Jx2[k][j][i] * Jx2[k][j][i] +
          Jx3[k][j][i] * Jx3[k][j][i];

      J = sqrt(J);

      nu_t = Dsdl2 * J;
      nu_max = MAX(nu_max, nu_t);
      dcoeff[k] = nu_t;

      scrh = 2.0 * vi[RHO] * nu_t;

      fBx = scrh * Jxz;
      fBy = scrh * Jyz;
      fBz = scrh * Jzz;

      // dH = (H[k + 1][j][i] - H[k][j][i]);
      // fe = nu_t / Prandtl * dH / dz;

      ViS[k][BX1] = 0.0;
      ViS[k][BX2] = 0.0;
      ViS[k][BX3] = 0.0;
      ViS[k][ENG] = 0.0;
#endif /* -- end #if GEOMETRY -- */

      /* ------------------------------------------
         3d. Compute fluxes at x3-faces
         ------------------------------------------ */

      ViF[k][BX1] = fBx;
      ViF[k][BX2] = fBy;
      ViF[k][BX3] = fBz;

      ViF[k][ENG] = 0;

#if HAVE_ENERGY
      // ViF[k][ENG] = fe;
#endif
    } /*loop*/
    // printf("z sweep: VIS[BX1] = %g\n", ViS[i][BX1]);
    // printf("z sweep: VIS[BX2] = %g\n", ViS[i][BX2]);
    // printf("z sweep: VIS[BX3] = %g\n", ViS[i][BX3]);
    // printf("z sweep: VIS[ENG] = %g\n", ViS[i][ENG]);
    // printf("z sweep: VIF[BX1] = %g\n", ViF[i][BX1]);
    // printf("z sweep: VIF[BX2] = %g\n", ViF[i][BX2]);
    // printf("z sweep: VIF[BX3] = %g\n", ViF[i][BX3]);
    // printf("z sweep: VIF[ENG] = %g\n", ViF[i][ENG]);

  } /*sweep*/
}
