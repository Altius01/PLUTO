#include "pluto.h"

#include "LES/les_phi.h"
#include <stdlib.h>

static double Ds = INIT_DS;

double get_MHD_LES_Ds() { return Ds; }
void set_MHD_LES_Ds(double ds) {
  Ds = ds > 0 ? MIN(ds, 0.3) : 0;
  // // printf("DS: %lf | ds: %lf\n", Ds, ds);
  g_maxDs = MAX(g_maxDs, Ds);
}

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
  double Jxx, Jyy, Jzz, Jxy, Jxz, Jyx, Jzx, Jyz, Jzy;
  Jxx = Jyy = Jzz = 0.0;

  double divJ, divr;
  double wx, wy, wz, qx, qy, qz;
  double fBx, fBy, fBz, fe;
  double nu_max;
  double Ds_ = 0.1, Dsdl2, dH;
  double Prandtl = 0.71;

#if DYNAMIC_PROCEDURE
  static double Lb_xy, Lb_xz, Lb_yz;
  static double Mb_xx, Mb_xy, Mb_xz, Mb_yx, Mb_yy, Mb_yz, Mb_zx, Mb_zy, Mb_zz;
#endif

  LES_phi_ctx *ctx = les_phi_create();

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

  dx = grid->dx[IDIR][IBEG];
  dy = grid->dx[JDIR][JBEG];
  dz = grid->dx[KDIR][KBEG];

  rho = d->Vc[RHO];
  // pr = d->Vc[PRS];

  Dsdl2 = Ds * Ds * pow(dx * dy * dz, 2.0 / 3.0);

  nu_max = 1.e-12;

  /* ----------------------------------------------
       5. Compute Specific enthalpy
       ---------------------------------------------- */

  init_les_phi_ctx(ctx, d, grid);

#if DYNAMIC_PROCEDURE
  double Ds_averaged_numerator = 0, Ds_averaged_denominator = 0;

  int i_h, j_h, k_h;
  KTOT_LOOP_FILTERED(k_h) {
    JTOT_LOOP_FILTERED(j_h) {
      ITOT_LOOP_FILTERED(i_h) {
        double divider = 0.0;
        double rho_f = 0, rhovx_f = 0, rhovy_f = 0, rhovz_f = 0, Bx_f = 0,
               By_f = 0, Bz_f = 0;

        double Phixx_f = 0, Phixy_f = 0, Phixz_f = 0, Phiyy_f = 0, Phiyx_f = 0,
               Phiyz_f = 0, Phizx_f = 0, Phizy_f = 0, Phizz_f = 0;

        double Jxy_f = 0, Jxz_f = 0, Jyz_f = 0;

        double Lb_VxBy_f = 0, Lb_VxBz_f = 0, Lb_VyBx_f = 0, Lb_VyBz_f = 0,
               Lb_VzBy_f = 0, Lb_VzBx_f = 0;

        double Mb_xy_f = 0, Mb_xz_f = 0, Mb_yx_f = 0, Mb_yz_f = 0, Mb_zy_f = 0,
               Mb_zx_f = 0;

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
              double Bx_ = Bx[k_shift][j_shift][i_shift];
              double By_ = By[k_shift][j_shift][i_shift];
              double Bz_ = Bz[k_shift][j_shift][i_shift];

              double local_volume = dx1[i_shift] * dx2[j_shift] * dx3[k_shift];

              rho_f += (rho_)*local_volume;

              rhovx_f += (rho_ * Vx_) * local_volume;
              rhovy_f += (rho_ * Vy_) * local_volume;
              rhovz_f += (rho_ * Vz_) * local_volume;

              Bx_f += (Bx_)*local_volume;
              By_f += (By_)*local_volume;
              Bz_f += (Bz_)*local_volume;

              Lb_VxBy_f += (rho_ * By_ * Vx_) * local_volume;
              Lb_VxBz_f += (rho_ * Bz_ * Vx_) * local_volume;
              Lb_VyBx_f += (rho_ * Bx_ * Vy_) * local_volume;
              Lb_VyBz_f += (rho_ * Bz_ * Vy_) * local_volume;
              Lb_VzBx_f += (rho_ * Bx_ * Vz_) * local_volume;
              Lb_VzBy_f += (rho_ * By_ * Vz_) * local_volume;

              NVAR_LOOP(nv) {
                vi[nv] = 0.5 * (d->Vc[nv][k][j][i] + d->Vc[nv][k][j][i + 1]);
                vc[nv] = d->Vc[nv][k][j][i];
              }

              update_les_phi_ctx(ctx, vi, grid, i_shift, j_shift, k_shift);

              Jxy_f += les_phi_get_Jxy(ctx) * local_volume;
              Jxz_f += les_phi_get_Jxz(ctx) * local_volume;
              Jyz_f += les_phi_get_Jyz(ctx) * local_volume;

              Phixx_f +=
                  LES_Phi_xx(ctx, i_shift, j_shift, k_shift) * local_volume;
              Phixy_f +=
                  LES_Phi_xy(ctx, i_shift, j_shift, k_shift) * local_volume;
              Phixz_f +=
                  LES_Phi_xz(ctx, i_shift, j_shift, k_shift) * local_volume;

              Phiyx_f +=
                  LES_Phi_yx(ctx, i_shift, j_shift, k_shift) * local_volume;
              Phiyy_f +=
                  LES_Phi_yy(ctx, i_shift, j_shift, k_shift) * local_volume;
              Phiyz_f +=
                  LES_Phi_yz(ctx, i_shift, j_shift, k_shift) * local_volume;

              Phizx_f +=
                  LES_Phi_zx(ctx, i_shift, j_shift, k_shift) * local_volume;
              Phizy_f +=
                  LES_Phi_zy(ctx, i_shift, j_shift, k_shift) * local_volume;
              Phizz_f +=
                  LES_Phi_zz(ctx, i_shift, j_shift, k_shift) * local_volume;

              Mb_xy_f += (Phixy_f * Jxy_f) * local_volume;
              Mb_yx_f -= (Phiyx_f * Jxy_f) * local_volume;
              Mb_xz_f += (Phixz_f * Jxz_f) * local_volume;
              Mb_zx_f -= (Phizx_f * Jxz_f) * local_volume;
              Mb_yz_f += (Phiyz_f * Jyz_f) * local_volume;
              Mb_zy_f -= (Phizy_f * Jyz_f) * local_volume;

              divider += local_volume;
            }
          }
        }

        if (divider == 0.0) {
          divider += 1e-20;
        }
        divider = 1.0 / divider;

        rho_f *= divider;

        if (rho_f == 0.0) {
          rho_f += 1e-20;
        }
        double inv_rho_f = 1.0 / rho_f;

        rhovx_f *= inv_rho_f * divider;
        rhovy_f *= inv_rho_f * divider;
        rhovz_f *= inv_rho_f * divider;

        Bx_f *= divider;
        By_f *= divider;
        Bz_f *= divider;

        Lb_VxBy_f *= inv_rho_f * divider;
        Lb_VxBz_f *= inv_rho_f * divider;
        Lb_VyBx_f *= inv_rho_f * divider;
        Lb_VyBz_f *= inv_rho_f * divider;
        Lb_VzBx_f *= inv_rho_f * divider;
        Lb_VzBy_f *= inv_rho_f * divider;

        Lb_xy = (Lb_VxBy_f - rhovx_f * inv_rho_f * By_f) -
                (Lb_VyBx_f - rhovy_f * inv_rho_f * Bx_f);

        Lb_xz = (Lb_VxBz_f - rhovx_f * inv_rho_f * Bz_f) -
                (Lb_VzBx_f - rhovz_f * inv_rho_f * Bx_f);

        Lb_yz = (Lb_VyBz_f - rhovy_f * inv_rho_f * Bz_f) -
                (Lb_VzBy_f - rhovz_f * inv_rho_f * By_f);

        Mb_xx = 0;
        Mb_xy = Phixy_f * Jxy_f - Mb_xy_f;
        Mb_xz = Phixz_f * Jxz_f - Mb_xz_f;

        Mb_yx = -Phiyx_f * Jxy_f - Mb_yx_f;
        Mb_yy = 0;
        Mb_yz = Phiyz_f * Jyz_f - Mb_yz_f;

        Mb_zx = -Phizx_f * Jxz_f - Mb_zx_f;
        Mb_zy = -Phizy_f * Jyz_f - Mb_zy_f;
        Mb_zz = 0;

        double local_volume = dx1[i_h] * dx2[j_h] * dx3[k_h];

        Ds_averaged_numerator +=
            (Lb_xy * Mb_xy + Lb_xz * Mb_xz - Lb_xy * Mb_yx + Lb_yz * Mb_yz -
             Lb_xz * Mb_zx - Lb_yz * Mb_zy) *
            local_volume;

        Ds_averaged_denominator +=
            (Mb_xx * Mb_xx + Mb_xy * Mb_xy + Mb_xz * Mb_xz + Mb_yx * Mb_yx +
             Mb_yy * Mb_yy + Mb_yz * Mb_yz + Mb_zx * Mb_zx + Mb_zy * Mb_zy +
             Mb_zz * Mb_zz) *
            local_volume;

        // if (Ds_averaged_numerator / (Ds_averaged_denominator + 1e-10) > 10) {
        //   printf("local_volume: %.10e\n", local_volume);
        //   printf("Lb_xy: %.10e, Lb_xz: %.10e, Lb_yz: %.10e \n", Lb_xy, Lb_xz,
        //          Lb_yz);
        //   printf("Mb_xx: %.10e, Mb_xy: %.10e, Mb_xz: %.10e, Mb_yx: %.10e, "
        //          "Mb_yy: %.10e, Mb_yz: %.10e, Mb_zx: %.10e, Mb_zy: %.10e\n",
        //          Mb_xx, Mb_xy, Mb_xz, Mb_yx, Mb_yy, Mb_yz, Mb_zx, Mb_zy);
        //   exit(1);
        // }
      }
    }
  }

  set_MHD_LES_Ds(Ds_averaged_numerator / (Ds_averaged_denominator + 1e-10));
#endif

  if (g_dir == IDIR) {

    /* ------------------------------------------------------
       1. Compute derivatives of velocity at x1-zone
          interfaces
       ------------------------------------------------------ */

    for (i = beg; i <= end; i++) {
      /* ------------------------------------------
         1c. Compute stress tensor components and
             geometrical source terms in different
             coordinate systems
         ------------------------------------------ */

#if GEOMETRY == CARTESIAN

      Ds_ = get_MHD_LES_Ds();
      nu_t = Ds_ * les_phi_get_eta(ctx);

      nu_max = MAX(nu_max, nu_t);
      dcoeff[i] = nu_max;

      fBx = 0;
      fBy = -2.0 * get_MHD_LES_Ds() * LES_Phi_xy(ctx, i, j, k) *
            les_phi_get_Jxy(ctx);
      fBz = -2.0 * get_MHD_LES_Ds() * LES_Phi_xz(ctx, i, j, k) *
            les_phi_get_Jxz(ctx);

      // dH = (H[k][j][i + 1] - H[k][j][i]);
      // fe = vi[RHO] * nu_t / Prandtl * dH / dx;

      ViS[i][BX1] = 0.0;
      ViS[i][BX2] = 0.0;
      ViS[i][BX3] = 0.0;
      // ViS[i][ENG] = 0.0;
#endif /* -- end #if GEOMETRY -- */

      /* ------------------------------------------
         1d. Compute fluxes at x1-faces
         ------------------------------------------ */

      ViF[i][BX1] = fBx;
      ViF[i][BX2] = fBy;
      ViF[i][BX3] = fBz;

      // ViF[i][ENG] = 0;
#if HAVE_ENERGY
      // ViF[i][ENG] = fe;
#endif
    }
  } else if (g_dir == JDIR) {

    /* ------------------------------------------------------
       2. Compute derivatives of velocity at x2-zone
          interfaces
       ------------------------------------------------------ */

    for (j = beg; j <= end; j++) {
      /* ------------------------------------------
         2c. Compute stress tensor components and
             geometrical source terms in different
             coordinate systems
         ------------------------------------------ */

#if GEOMETRY == CARTESIAN

      Ds_ = get_MHD_LES_Ds();
      nu_t = Ds_ * les_phi_get_eta(ctx);

      nu_max = MAX(nu_max, nu_t);
      dcoeff[j] = nu_max;

      fBx = -2.0 * get_MHD_LES_Ds() * LES_Phi_yx(ctx, i, j, k) *
            les_phi_get_Jyx(ctx);
      fBy = 0;
      fBz = -2.0 * get_MHD_LES_Ds() * LES_Phi_yz(ctx, i, j, k) *
            les_phi_get_Jyz(ctx);

      // dH = (H[k][j + 1][i] - H[k][j][i]);
      // fe = vi[RHO] * nu_t / Prandtl * dH / dy;

      ViS[j][BX1] = 0.0;
      ViS[j][BX2] = 0.0;
      ViS[j][BX3] = 0.0;
      // ViS[j][ENG] = 0.0;
#endif /* -- GEOMETRY -- */

      /* ------------------------------------------
         2d. Compute fluxes at x2-faces
         ------------------------------------------ */

      ViF[j][BX1] = fBx;
      ViF[j][BX2] = fBy;
      ViF[j][BX3] = fBz;

      // ViF[j][ENG] = 0;
#if HAVE_ENERGY
      // ViF[j][ENG] = fe;
#endif
    }

  } else if (g_dir == KDIR) {

    /* ------------------------------------------------------
       3. Compute derivatives of velocity at x3-zone
          interfaces
       ------------------------------------------------------ */

    for (k = beg; k <= end; k++) {
      /* ------------------------------------------
         3c. Compute stress tensor components and
             geometrical source terms in different
             coordinate systems
         ------------------------------------------ */

#if GEOMETRY == CARTESIAN

      Ds_ = get_MHD_LES_Ds();
      nu_t = Ds_ * les_phi_get_eta(ctx);

      nu_max = MAX(nu_max, nu_t);
      dcoeff[k] = nu_max;

      fBx = -2.0 * get_MHD_LES_Ds() * LES_Phi_zx(ctx, i, j, k) *
            les_phi_get_Jzx(ctx);
      fBy = -2.0 * get_MHD_LES_Ds() * LES_Phi_zy(ctx, i, j, k) *
            les_phi_get_Jzy(ctx);
      fBz = 0;

      // dH = (H[k + 1][j][i] - H[k][j][i]);
      // fe = nu_t / Prandtl * dH / dz;

      ViS[k][BX1] = 0.0;
      ViS[k][BX2] = 0.0;
      ViS[k][BX3] = 0.0;
      // ViS[k][ENG] = 0.0;
#endif /* -- end #if GEOMETRY -- */

      /* ------------------------------------------
         3d. Compute fluxes at x3-faces
         ------------------------------------------ */

      ViF[k][BX1] = fBx;
      ViF[k][BX2] = fBy;
      ViF[k][BX3] = fBz;

      // ViF[k][ENG] = 0;

#if HAVE_ENERGY
      // ViF[k][ENG] = fe;
#endif
    } /*loop*/
  } /*sweep*/
  les_phi_delete(ctx);
}
