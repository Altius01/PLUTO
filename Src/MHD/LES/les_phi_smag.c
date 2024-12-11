#include "LES/les.h"
#include "pluto.h"

#include "les_phi_smag.h"

void init_les_phi_ctx(LES_phi_ctx *ctx, const Data *d, Grid *grid) {
  double dx = grid->dx[IDIR][IBEG];
  double dy = grid->dx[JDIR][JBEG];
  double dz = grid->dx[KDIR][KBEG];

  ctx->dl2 = pow(dx * dy * dz, 2.0 / 3.0);

  ctx->x1 = grid->x[IDIR], ctx->dx1 = grid->dx[IDIR];
  ctx->x2 = grid->x[JDIR], ctx->dx2 = grid->dx[JDIR];
  ctx->x3 = grid->x[KDIR], ctx->dx3 = grid->dx[KDIR];

  ctx->Bx = d->Vc[BX1];
  ctx->By = d->Vc[BX2];
  ctx->Bz = d->Vc[BX3];

  ctx->Jx1 = d->J[IDIR];
  ctx->Jx2 = d->J[JDIR];
  ctx->Jx3 = d->J[KDIR];
}

void update_les_phi_ctx(LES_phi_ctx *ctx, double *vi, Grid *grid, int i, int j,
                        int k) {

  double inv_dx1 = 1.0 / ctx->dx1[i];
  double inv_dx3 = 1.0 / ctx->dx3[k];
  double inv_dx2 = 1.0 / ctx->dx2[j];

  /* -- Compute viscosity and velocity derivatives -- */
#if INCLUDE_IDIR
  ctx->dxBx = FDIFF_X1(ctx->Bx, k, j, i) * inv_dx1;
  ctx->dxBy = FDIFF_X1(ctx->By, k, j, i) * inv_dx1;
  ctx->dxBz = FDIFF_X1(ctx->Bz, k, j, i) * inv_dx1;
#else
  ctx->dxBx = 0;
  ctx->dxBy = 0;
  ctx->dxBz = 0;
#endif
#if INCLUDE_JDIR
  ctx->dyBx = 0.5 *
              (CDIFF_X2(ctx->Bx, k, j, i) + CDIFF_X2(ctx->Bx, k, j, i + 1)) *
              inv_dx2;
  ctx->dyBy = 0.5 *
              (CDIFF_X2(ctx->By, k, j, i) + CDIFF_X2(ctx->By, k, j, i + 1)) *
              inv_dx2;
  ctx->dyBz = 0.5 *
              (CDIFF_X2(ctx->Bz, k, j, i) + CDIFF_X2(ctx->Bz, k, j, i + 1)) *
              inv_dx2;
#else
  ctx->dyBx = 0;
  ctx->dyBy = 0;
  ctx->dyBz = 0;

#endif
#if INCLUDE_KDIR
  ctx->dzBx = 0.5 *
              (CDIFF_X3(ctx->Bx, k, j, i) + CDIFF_X3(ctx->Bx, k, j, i + 1)) *
              inv_dx3;
  ctx->dzBy = 0.5 *
              (CDIFF_X3(ctx->By, k, j, i) + CDIFF_X3(ctx->By, k, j, i + 1)) *
              inv_dx3;
  ctx->dzBz = 0.5 *
              (CDIFF_X3(ctx->Bz, k, j, i) + CDIFF_X3(ctx->Bz, k, j, i + 1)) *
              inv_dx3;
#else
  ctx->dzBx = 0;
  ctx->dzBy = 0;
  ctx->dzBz = 0;

#endif

  ctx->Jp[IDIR] = ctx->dyBz - ctx->dzBy;
  ctx->Jp[JDIR] = ctx->dzBx - ctx->dxBz;
  ctx->Jp[KDIR] = ctx->dxBy - ctx->dyBx;

  ctx->j = ctx->Jp[IDIR] * ctx->Jp[IDIR] + ctx->Jp[JDIR] * ctx->Jp[JDIR] +
           ctx->Jp[KDIR] * ctx->Jp[KDIR];
  ctx->j = sqrt(ctx->j);

  ctx->Jxy = 0.5 * (ctx->dyBx - ctx->dxBy);
  ctx->Jxz = 0.5 * (ctx->dzBx - ctx->dxBz);
  ctx->Jyz = 0.5 * (ctx->dzBy - ctx->dyBz);

  ctx->eta_t = ctx->dl2 * ctx->j;

  ctx->scrh = ctx->eta_t;
}
