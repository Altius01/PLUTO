#ifndef LES_PHI_SMAG_H
#define LES_PHI_SMAG_H

#include "pluto.h"

typedef struct LES_phi_ctx {
  double *x1, *dx1;
  double *x2, *dx2;
  double *x3, *dx3;

  double ***Bx;
  double ***By;
  double ***Bz;

  double ***Jx1;
  double ***Jx2;
  double ***Jx3;

  double vp[3], Jp[3];

  double scrh;
  double eta_t;

  double dl2;
  double j;
  double Jxy, Jxz, Jyz;
  double dxBx, dxBy, dxBz, dyBx, dyBy, dyBz, dzBx, dzBy, dzBz;
} LES_phi_ctx;

LES_phi_ctx *les_phi_create() {
  LES_phi_ctx *result;

  result = (LES_phi_ctx *)malloc(sizeof(LES_phi_ctx));
  if (result == NULL) {
    return NULL;
  }
  return result;
}

void les_phi_delete(LES_phi_ctx *ctx) { free(ctx); }

double les_phi_get_eta(LES_phi_ctx *ctx) { return ctx->eta_t; }

double les_phi_get_Jxy(LES_phi_ctx *ctx) { return ctx->Jxy; }
double les_phi_get_Jxz(LES_phi_ctx *ctx) { return ctx->Jxz; }
double les_phi_get_Jyz(LES_phi_ctx *ctx) { return ctx->Jyz; }

double les_phi_get_Jyx(LES_phi_ctx *ctx) { return -ctx->Jxy; }
double les_phi_get_Jzx(LES_phi_ctx *ctx) { return -ctx->Jxz; }
double les_phi_get_Jzy(LES_phi_ctx *ctx) { return -ctx->Jyz; }

double LES_Phi_xx(LES_phi_ctx *ctx, int i, int j, int k) { return ctx->scrh; };

double LES_Phi_xy(LES_phi_ctx *ctx, int i, int j, int k) { return ctx->scrh; }

double LES_Phi_xz(LES_phi_ctx *ctx, int i, int j, int k) { return ctx->scrh; }

double LES_Phi_yx(LES_phi_ctx *ctx, int i, int j, int k) { return ctx->scrh; }

double LES_Phi_yy(LES_phi_ctx *ctx, int i, int j, int k) { return ctx->scrh; }

double LES_Phi_yz(LES_phi_ctx *ctx, int i, int j, int k) { return ctx->scrh; }

double LES_Phi_zx(LES_phi_ctx *ctx, int i, int j, int k) { return ctx->scrh; }

double LES_Phi_zy(LES_phi_ctx *ctx, int i, int j, int k) { return ctx->scrh; }

double LES_Phi_zz(LES_phi_ctx *ctx, int i, int j, int k) { return ctx->scrh; }
#endif
