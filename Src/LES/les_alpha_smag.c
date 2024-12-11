#include "LES/les.h"
#include "pluto.h"

#include "les_alpha_smag.h"

void init_les_alpha_ctx(LES_alpha_ctx *ctx, const Data *d, Grid *grid) {
  double dx = grid->dx[IDIR][IBEG];
  double dy = grid->dx[JDIR][JBEG];
  double dz = grid->dx[KDIR][KBEG];

  ctx->dl2 = pow(dx * dy * dz, 2.0 / 3.0);

  ctx->x1 = grid->x[IDIR], ctx->dx1 = grid->dx[IDIR];
  ctx->x2 = grid->x[JDIR], ctx->dx2 = grid->dx[JDIR];
  ctx->x3 = grid->x[KDIR], ctx->dx3 = grid->dx[KDIR];

  ctx->Vx = d->Vc[VX1];
  ctx->Vy = d->Vc[VX2];
  ctx->Vz = d->Vc[VX3];
}

void update_les_alpha_ctx(LES_alpha_ctx *ctx, double *vi, Grid *grid, int i,
                          int j, int k) {

  double inv_dx1 = 1.0 / ctx->dx1[i];
  double inv_dx3 = 1.0 / ctx->dx3[k];
  double inv_dx2 = 1.0 / ctx->dx2[j];

  /* -- Compute viscosity and velocity derivatives -- */
#if INCLUDE_IDIR
  ctx->dxVx = FDIFF_X1(ctx->Vx, k, j, i) * inv_dx1;
  ctx->dxVy = FDIFF_X1(ctx->Vy, k, j, i) * inv_dx1;
  ctx->dxVz = FDIFF_X1(ctx->Vz, k, j, i) * inv_dx1;
#endif
#if INCLUDE_JDIR
  ctx->dyVx = 0.5 *
              (CDIFF_X2(ctx->Vx, k, j, i) + CDIFF_X2(ctx->Vx, k, j, i + 1)) *
              inv_dx2;
  ctx->dyVy = 0.5 *
              (CDIFF_X2(ctx->Vy, k, j, i) + CDIFF_X2(ctx->Vy, k, j, i + 1)) *
              inv_dx2;
  ctx->dyVz = 0.5 *
              (CDIFF_X2(ctx->Vz, k, j, i) + CDIFF_X2(ctx->Vz, k, j, i + 1)) *
              inv_dx2;
#endif
#if INCLUDE_KDIR
  ctx->dzVx = 0.5 *
              (CDIFF_X3(ctx->Vx, k, j, i) + CDIFF_X3(ctx->Vx, k, j, i + 1)) *
              inv_dx3;
  ctx->dzVy = 0.5 *
              (CDIFF_X3(ctx->Vy, k, j, i) + CDIFF_X3(ctx->Vy, k, j, i + 1)) *
              inv_dx3;
  ctx->dzVz = 0.5 *
              (CDIFF_X3(ctx->Vz, k, j, i) + CDIFF_X3(ctx->Vz, k, j, i + 1)) *
              inv_dx3;
#endif

  ctx->Sxx = 0.5 * (ctx->dxVx + ctx->dxVx);
  ctx->Sxy = 0.5 * (ctx->dyVx + ctx->dxVy);
  ctx->Sxz = 0.5 * (ctx->dzVx + ctx->dxVz);

  ctx->Syy = 0.5 * (ctx->dyVy + ctx->dyVy);
  ctx->Syz = 0.5 * (ctx->dzVy + ctx->dyVz);
  ctx->Szz = 0.5 * (ctx->dzVz + ctx->dzVz);

  ctx->S =
      ctx->Sxx * ctx->Sxx + ctx->Syy * ctx->Syy + ctx->Szz * ctx->Szz +
      2.0 * (ctx->Sxy * ctx->Sxy + ctx->Sxz * ctx->Sxz + ctx->Syz * ctx->Syz);
  ctx->S = sqrt(2.0 * ctx->S);

  ctx->divV = DIM_EXPAND(ctx->dxVx, +ctx->dyVy, +ctx->dzVz);

  ctx->rho = vi[RHO];

  ctx->nu_t = ctx->dl2 * ctx->S;
  ctx->nu_t_prime = ctx->dl2 * ctx->S;

  ctx->scrh = ctx->rho * ctx->nu_t;

  ctx->scrh_trace = ctx->rho * ctx->nu_t_prime;
}
