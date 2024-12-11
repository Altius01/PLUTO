#include "pluto.h"

#ifndef LES_FMIN
#define LES_FMIN 0.5 /* blah blah  */
#endif

static double Cs = INIT_CS, Ys = INIT_DS;

double get_LES_Cs() { return Cs; }

double get_LES_Ys() { return Ys; }

void set_LES_Cs(double cs) {
  Cs = cs > 0 ? MIN(cs, 1) : 0;
  // printf("CS: %lf | cs: %lf\n", Cs, cs);
  g_maxCs = MAX(g_maxCs, Cs);
}

void set_LES_Ys(double ys) {
  Ys = ys;
  // printf("YS: %lf | ys: %lf\n", Ys, ys);
  g_maxYs = MAX(g_maxYs, Ys);
}

#define ITOT_LOOP_FILTERED(i) for ((i) = 0; (i) < NX1_TOT; (i) += 2)
#define JTOT_LOOP_FILTERED(j) for ((j) = 0; (j) < NX2_TOT; (j) += 2)
#define KTOT_LOOP_FILTERED(k) for ((k) = 0; (k) < NX3_TOT; (k) += 2)

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
  double nu_t, nu_t_trace;
  double dx, dy, dz;
  double inv_dx1, inv_dx2, inv_dx3;
  double *x1 = grid->x[IDIR], *dx1 = grid->dx[IDIR];
  double *x2 = grid->x[JDIR], *dx2 = grid->dx[JDIR];
  double *x3 = grid->x[KDIR], *dx3 = grid->dx[KDIR];
  double *x1r = grid->xr[IDIR];
  double *x2r = grid->xr[JDIR];
  double *x3r = grid->xr[KDIR];
  double dxVx, dxVy, dxVz, dyVx, dyVy, dyVz, dzVx, dzVy, dzVz;
  double dxBx, dxBy, dxBz, dyBx, dyBy, dyBz, dzBx, dzBy, dzBz;

  double S, Sxx, Syy, Szz, Sxy, Sxz, Syz;
  double divr;
  double wx, wy, wz, qx, qy, qz;
  double fmx, fmy, fmz, fe;
  double nu_max;
  double Cs_, Ys_, Csdl2, dH;
  double Prandtl = 0.71;

  // static double ***H;

#if DYNAMIC_PROCEDURE
  static double Lu_xx, Lu_yy, Lu_zz, Lu_xy, Lu_xz, Lu_yz;
  static double Mu_xx, Mu_xy, Mu_xz, Mu_yx, Mu_yy, Mu_yz, Mu_zx, Mu_zy, Mu_zz;

  static double ALPHA_Y;

#endif
  LES_alpha_ctx *ctx = les_alpha_create();

  double ***Vx = d->Vc[VX1];
  double ***Vy = d->Vc[VX2];
  double ***Vz = d->Vc[VX3];

#if PHYSICS == MHD || PHYSICS == RMHD
  double ***Bx = d->Vc[BX1];
  double ***By = d->Vc[BX2];
  double ***Bz = d->Vc[BX3];
#endif

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

  // if (H == NULL) {
  // H = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  // }

  dx = grid->dx[IDIR][IBEG];
  dy = grid->dx[JDIR][JBEG];
  dz = grid->dx[KDIR][KBEG];

  rho = d->Vc[RHO];

  // pr = d->Vc[PRS];

  nu_max = 1.e-12;

  /* ----------------------------------------------
       5. Compute Specific enthalpy
       ---------------------------------------------- */

  init_les_alpha_ctx(ctx, d, grid);
#if DYNAMIC_PROCEDURE

  int i_h, j_h, k_h;

  double Cs_averaged_numerator = 0, Cs_averaged_denominator = 0;
  double Ys_averaged_numerator = 0, Ys_averaged_denominator = 0;

  KTOT_LOOP_FILTERED(k_h) {
    JTOT_LOOP_FILTERED(j_h) {
      ITOT_LOOP_FILTERED(i_h) {
        double divider = 0.0;

        double rho_f = 0, rhovx_f = 0, rhovy_f = 0, rhovz_f = 0, Bx_f = 0,
               By_f = 0, Bz_f = 0;
        double BxBx_f = 0, ByBy_f = 0, BzBz_f = 0, BxBy_f = 0, BxBz_f = 0,
               ByBz_f = 0;

        double divV_f = 0, S_f = 0, Sxx_f = 0, Sxy_f = 0, Sxz_f = 0, Syy_f = 0,
               Syz_f = 0, Szz_f = 0;

        double Lu_Vxx_f = 0, Lu_Vxy_f = 0, Lu_Vxz_f = 0, Lu_Vyy_f = 0,
               Lu_Vyz_f = 0, Lu_Vzz_f = 0;
        double Mu_Sxx_f = 0, Mu_Sxy_f = 0, Mu_Sxz_f = 0, Mu_Syy_f = 0,
               Mu_Syx_f = 0, Mu_Syz_f = 0, Mu_Szx_f = 0, Mu_Szy_f = 0,
               Mu_Szz_f = 0;

        double A_f = 0, Axx_f = 0, Axy_f = 0, Axz_f = 0, Ayy_f = 0, Ayx_f = 0,
               Ayz_f = 0, Azx_f = 0, Azy_f = 0, Azz_f = 0;

        double Y_f = 0;

        int delta_i = 0, delta_j = 0, delta_k = 0;

#if INCLUDE_IDIR
        delta_i = 1;
#endif
#if INCLUDE_JDIR
        delta_j = 1;
#endif
#if INCLUDE_KDIR
        delta_k = 1;
#endif

        for (int i_shift = i_h; i_shift <= i_h + delta_i; i_shift++) {
          for (int j_shift = j_h; j_shift <= j_h + delta_j; j_shift++) {
            for (int k_shift = k_h; k_shift <= k_h + delta_k; k_shift++) {
              double rho_ = rho[k_shift][j_shift][i_shift];
              double Vx_ = Vx[k_shift][j_shift][i_shift];
              double Vy_ = Vy[k_shift][j_shift][i_shift];
              double Vz_ = Vz[k_shift][j_shift][i_shift];

              double local_volume = dx1[i_shift] * dx2[j_shift] * dx3[k_shift];

              rho_f += (rho_)*local_volume;

              rhovx_f += (rho_ * Vx_) * local_volume;
              rhovy_f += (rho_ * Vy_) * local_volume;
              rhovz_f += (rho_ * Vz_) * local_volume;

#if PHYSICS == MHD || PHYSICS == RMHD
              double Bx_ = Bx[k_shift][j_shift][i_shift];
              double By_ = By[k_shift][j_shift][i_shift];
              double Bz_ = Bz[k_shift][j_shift][i_shift];

              Bx_f += (Bx_)*local_volume;
              By_f += (By_)*local_volume;
              Bz_f += (Bz_)*local_volume;

              BxBx_f += (Bx_ * Bx_) * local_volume;
              ByBy_f += (By_ * By_) * local_volume;
              BzBz_f += (Bz_ * Bz_) * local_volume;
              BxBy_f += (Bx_ * By_) * local_volume;
              BxBz_f += (Bx_ * Bz_) * local_volume;
              ByBz_f += (By_ * Bz_) * local_volume;
#endif

              Lu_Vxx_f += (rho_ * Vx_ * Vx_) * local_volume;
              Lu_Vxy_f += (rho_ * Vx_ * Vy_) * local_volume;
              Lu_Vxz_f += (rho_ * Vx_ * Vz_) * local_volume;
              Lu_Vyy_f += (rho_ * Vy_ * Vy_) * local_volume;
              Lu_Vyz_f += (rho_ * Vy_ * Vz_) * local_volume;
              Lu_Vzz_f += (rho_ * Vz_ * Vz_) * local_volume;

              NVAR_LOOP(nv) {
                vi[nv] = 0.5 * (d->Vc[nv][k][j][i] + d->Vc[nv][k][j][i + 1]);
                vc[nv] = d->Vc[nv][k][j][i];
              }

              update_les_alpha_ctx(ctx, vi, grid, i_shift, j_shift, k_shift);

              Sxx_f += rho_ * les_alpha_get_Sxx(ctx) * local_volume;
              Sxy_f += rho_ * les_alpha_get_Sxy(ctx) * local_volume;
              Sxz_f += rho_ * les_alpha_get_Sxz(ctx) * local_volume;
              Syy_f += rho_ * les_alpha_get_Syy(ctx) * local_volume;
              Syz_f += rho_ * les_alpha_get_Syz(ctx) * local_volume;
              Szz_f += rho_ * les_alpha_get_Szz(ctx) * local_volume;
              S_f += rho_f * les_alpha_get_S(ctx) * local_volume;
              divV_f += rho_f * les_alpha_get_divV(ctx) * local_volume;

              A_f += LES_Alpha(ctx, i_shift, j_shift, k_shift) * local_volume;

              Axx_f +=
                  LES_Alpha_xx(ctx, i_shift, j_shift, k_shift) * local_volume;
              Axy_f +=
                  LES_Alpha_xy(ctx, i_shift, j_shift, k_shift) * local_volume;
              Axz_f +=
                  LES_Alpha_xz(ctx, i_shift, j_shift, k_shift) * local_volume;

              Ayx_f +=
                  LES_Alpha_yx(ctx, i_shift, j_shift, k_shift) * local_volume;
              Ayy_f +=
                  LES_Alpha_yy(ctx, i_shift, j_shift, k_shift) * local_volume;
              Ayz_f +=
                  LES_Alpha_yz(ctx, i_shift, j_shift, k_shift) * local_volume;

              Azx_f +=
                  LES_Alpha_zx(ctx, i_shift, j_shift, k_shift) * local_volume;
              Azy_f +=
                  LES_Alpha_zy(ctx, i_shift, j_shift, k_shift) * local_volume;
              Azz_f +=
                  LES_Alpha_zz(ctx, i_shift, j_shift, k_shift) * local_volume;

              Y_f += rho_ * LES_Alpha(ctx, i_shift, j_shift, k_shift) *
                     (les_alpha_get_S(ctx))*local_volume;

              Mu_Sxx_f +=
                  LES_Alpha_xx(ctx, i_shift, j_shift, k_shift) * rho_ *
                  (les_alpha_get_Sxx(ctx) - c13 * les_alpha_get_divV(ctx)) *
                  local_volume;
              Mu_Sxy_f += LES_Alpha_xy(ctx, i_shift, j_shift, k_shift) * rho_ *
                          (les_alpha_get_Sxy(ctx))*local_volume;
              Mu_Sxz_f += LES_Alpha_xz(ctx, i_shift, j_shift, k_shift) * rho_ *
                          (les_alpha_get_Sxz(ctx))*local_volume;
              Mu_Syx_f += LES_Alpha_yx(ctx, i_shift, j_shift, k_shift) * rho_ *
                          (les_alpha_get_Sxy(ctx))*local_volume;
              Mu_Syy_f +=
                  LES_Alpha_yy(ctx, i_shift, j_shift, k_shift) * rho_ *
                  (les_alpha_get_Syy(ctx) - c13 * les_alpha_get_divV(ctx)) *
                  local_volume;
              Mu_Syz_f += LES_Alpha_yz(ctx, i_shift, j_shift, k_shift) * rho_ *
                          (les_alpha_get_Syz(ctx))*local_volume;
              Mu_Szx_f += LES_Alpha_zx(ctx, i_shift, j_shift, k_shift) * rho_ *
                          (les_alpha_get_Sxz(ctx))*local_volume;
              Mu_Szy_f += LES_Alpha_zy(ctx, i_shift, j_shift, k_shift) * rho_ *
                          (les_alpha_get_Syz(ctx))*local_volume;
              Mu_Szz_f +=
                  LES_Alpha_zz(ctx, i_shift, j_shift, k_shift) * rho_ *
                  (les_alpha_get_Szz(ctx) - c13 * les_alpha_get_divV(ctx)) *
                  local_volume;

              divider += local_volume;
            }
          }
        }

        if (divider == 0.0)
          divider += 1e-20;

        divider = 1.0 / divider;

        rho_f *= divider;

        if (rho_f == 0.0)
          rho_f += 1e-20;

        double inv_rho_f = 1.0 / rho_f;

        rhovx_f *= divider;
        rhovy_f *= divider;
        rhovz_f *= divider;

        Sxx_f *= inv_rho_f * divider;
        Sxy_f *= inv_rho_f * divider;
        Sxz_f *= inv_rho_f * divider;
        Syy_f *= inv_rho_f * divider;
        Syz_f *= inv_rho_f * divider;
        Szz_f *= inv_rho_f * divider;
        S_f *= inv_rho_f * divider;
        divV_f *= inv_rho_f * divider;

        Y_f *= inv_rho_f * divider;

#if PHYSICS == MHD || PHYSICS == RMHD
        Bx_f *= divider;
        By_f *= divider;
        Bz_f *= divider;

        BxBx_f *= divider;
        ByBy_f *= divider;
        BzBz_f *= divider;
        BxBy_f *= divider;
        BxBz_f *= divider;
        ByBz_f *= divider;
#endif

        Lu_Vxx_f *= divider;
        Lu_Vxy_f *= divider;
        Lu_Vxz_f *= divider;
        Lu_Vyy_f *= divider;
        Lu_Vyz_f *= divider;
        Lu_Vzz_f *= divider;

        A_f *= divider;
        Axx_f *= divider;
        Axy_f *= divider;
        Axz_f *= divider;
        Ayx_f *= divider;
        Ayy_f *= divider;
        Ayz_f *= divider;
        Azx_f *= divider;
        Azy_f *= divider;
        Azz_f *= divider;

        Mu_Sxx_f *= inv_rho_f * divider;
        Mu_Sxy_f *= inv_rho_f * divider;
        Mu_Sxz_f *= inv_rho_f * divider;
        Mu_Syx_f *= inv_rho_f * divider;
        Mu_Syy_f *= inv_rho_f * divider;
        Mu_Syz_f *= inv_rho_f * divider;
        Mu_Szx_f *= inv_rho_f * divider;
        Mu_Szy_f *= inv_rho_f * divider;
        Mu_Szz_f *= inv_rho_f * divider;

#if PHYSICS == MHD || PHYSICS == RMHD
        Lu_xx =
            Lu_Vxx_f - rhovx_f * rhovx_f * inv_rho_f - (BxBx_f - Bx_f * Bx_f);

        Lu_yy =
            Lu_Vyy_f - rhovy_f * rhovy_f * inv_rho_f - (ByBy_f - By_f * By_f);

        Lu_zz =
            Lu_Vzz_f - rhovz_f * rhovz_f * inv_rho_f - (BzBz_f - Bz_f * Bz_f);

        Lu_xy =
            Lu_Vxy_f - rhovx_f * rhovy_f * inv_rho_f - (BxBy_f - Bx_f * By_f);

        Lu_xz =
            Lu_Vxz_f - rhovx_f * rhovz_f * inv_rho_f - (BxBz_f - Bx_f * Bz_f);

        Lu_yz =
            Lu_Vyz_f - rhovy_f * rhovz_f * inv_rho_f - (ByBz_f - By_f * Bz_f);

#else
        Lu_xx = Lu_Vxx_f - rhovx_f * rhovx_f * inv_rho_f;
        Lu_yy = Lu_Vyy_f - rhovy_f * rhovy_f * inv_rho_f;
        Lu_zz = Lu_Vzz_f - rhovz_f * rhovz_f * inv_rho_f;
        Lu_xy = Lu_Vxy_f - rhovx_f * rhovy_f * inv_rho_f;
        Lu_xz = Lu_Vxz_f - rhovx_f * rhovz_f * inv_rho_f;
        Lu_yz = Lu_Vyz_f - rhovy_f * rhovz_f * inv_rho_f;
#endif

        Mu_xx = Axx_f * (Sxx_f - c13 * divV_f) - Mu_Sxx_f;
        Mu_xy = Axy_f * (Sxy_f)-Mu_Sxy_f;
        Mu_xz = Axz_f * (Sxz_f)-Mu_Sxz_f;

        Mu_yx = Ayx_f * (Sxy_f)-Mu_Syx_f;
        Mu_yy = Ayy_f * (Syy_f - c13 * divV_f) - Mu_Syy_f;
        Mu_yz = Ayz_f * (Syz_f)-Mu_Syz_f;

        Mu_zx = Azx_f * (Sxz_f)-Mu_Szx_f;
        Mu_zy = Azy_f * (Syz_f)-Mu_Szy_f;
        Mu_zz = Azz_f * (Szz_f - c13 * divV_f) - Mu_Szz_f;

        ALPHA_Y = A_f * S_f - Y_f;

        double local_volume = dx1[i_h] * dx2[j_h] * dx3[k_h];

        Cs_averaged_numerator +=
            (Lu_xx * Mu_xx + Lu_xy * Mu_xy + Lu_xz * Mu_xz + Lu_xy * Mu_yx +
             Lu_yy * Mu_yy + Lu_yz * Mu_yz + Lu_xz * Mu_zx + Lu_yz * Mu_zy +
             Lu_zz * Mu_zz) *
            local_volume;

        Cs_averaged_denominator +=
            (Mu_xx * Mu_xx + Mu_xy * Mu_xy + Mu_xz * Mu_xz + Mu_yx * Mu_yx +
             Mu_yy * Mu_yy + Mu_yz * Mu_yz + Mu_zx * Mu_zx + Mu_zy * Mu_zy +
             Mu_zz * Mu_zz) *
            local_volume;

        Ys_averaged_numerator +=
            (Lu_xx * Lu_xx + Lu_yy * Lu_yy + Lu_zz * Lu_zz) * local_volume;

        Ys_averaged_denominator += (ALPHA_Y)*local_volume;
      }
    }
  }

  set_LES_Cs(Cs_averaged_numerator / (Cs_averaged_denominator + 1e-10));

  // printf(
  //     "Ys: %.10e Ys_averaged_numerator: %.10e Ys_averaged_denominator:
  //     %.10e\n", Ys_averaged_numerator / (Ys_averaged_denominator + 1e-10),
  //     Ys_averaged_numerator, Ys_averaged_denominator);

  set_LES_Ys(Ys_averaged_numerator / (Ys_averaged_denominator + 1e-10));
#endif

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

      update_les_alpha_ctx(ctx, vi, grid, i, j, k);

      /* ------------------------------------------
         1c. Compute stress tensor components and
             geometrical source terms in different
             coordinate systems
         ------------------------------------------ */

#if GEOMETRY == CARTESIAN

      Cs_ = get_LES_Cs();
      Ys_ = get_LES_Ys();

      nu_t = Cs_ * les_alpha_get_nut(ctx);
      nu_t_trace = Ys_ * les_alpha_get_nut_trace(ctx);

      nu_max = MAX(nu_max, MAX(nu_t, nu_t_trace));
      dcoeff[i] = fabs(nu_max);

      fmx = -2.0 * Cs_ * LES_Alpha_xx(ctx, i, j, k) *
                (les_alpha_get_Sxx(ctx) - c13 * les_alpha_get_divV(ctx)) +
            2.0 * c13 * Ys_ * LES_Alpha(ctx, i, j, k) * les_alpha_get_S(ctx);
      fmy = -2.0 * Cs_ * LES_Alpha_xy(ctx, i, j, k) * les_alpha_get_Sxy(ctx);
      fmz = -2.0 * Cs_ * LES_Alpha_xz(ctx, i, j, k) * les_alpha_get_Sxz(ctx);

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

      update_les_alpha_ctx(ctx, vi, grid, i, j, k);

      /* ------------------------------------------
         2c. Compute stress tensor components and
             geometrical source terms in different
             coordinate systems
         ------------------------------------------ */

#if GEOMETRY == CARTESIAN
      Cs_ = get_LES_Cs();
      Ys_ = get_LES_Ys();

      nu_t = Cs_ * les_alpha_get_nut(ctx);
      nu_t_trace = Ys_ * les_alpha_get_nut_trace(ctx);

      nu_max = MAX(nu_max, MAX(nu_t, nu_t_trace));
      dcoeff[j] = fabs(nu_max);

      fmx = -2.0 * Cs_ * LES_Alpha_yz(ctx, i, j, k) * les_alpha_get_Sxy(ctx);
      fmy = -2.0 * Cs_ * LES_Alpha_yy(ctx, i, j, k) *
                (les_alpha_get_Syy(ctx) - c13 * les_alpha_get_divV(ctx)) +
            2.0 * c13 * Ys_ * LES_Alpha(ctx, i, j, k) * les_alpha_get_S(ctx);
      fmz = -2.0 * Cs_ * LES_Alpha_yz(ctx, i, j, k) * les_alpha_get_Syz(ctx);

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

      update_les_alpha_ctx(ctx, vi, grid, i, j, k);

      /* ------------------------------------------
         3c. Compute stress tensor components and
             geometrical source terms in different
             coordinate systems
         ------------------------------------------ */

#if GEOMETRY == CARTESIAN
      Cs_ = get_LES_Cs();
      Ys_ = get_LES_Ys();

      nu_t = Cs_ * les_alpha_get_nut(ctx);
      nu_t_trace = Ys_ * les_alpha_get_nut_trace(ctx);

      nu_max = MAX(nu_max, MAX(nu_t, nu_t_trace));
      dcoeff[k] = fabs(nu_max);

      fmx = -2.0 * Cs_ * LES_Alpha_zx(ctx, i, j, k) * les_alpha_get_Sxz(ctx);
      fmy = -2.0 * Cs_ * LES_Alpha_zy(ctx, i, j, k) * les_alpha_get_Syz(ctx);
      fmz = -2.0 * Cs_ * LES_Alpha_zz(ctx, i, j, k) *
                (les_alpha_get_Szz(ctx) - c13 * les_alpha_get_divV(ctx)) +
            2.0 * c13 * Ys_ * LES_Alpha(ctx, i, j, k) * les_alpha_get_S(ctx);

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
  les_alpha_delete(ctx);
}
