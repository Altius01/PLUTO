#ifndef LES_ALPHA_SMAG_H
#define LES_ALPHA_SMAG_H

#include "pluto.h"

typedef struct LES_alpha_ctx {
  double *x1, *dx1;
  double *x2, *dx2;
  double *x3, *dx3;

  double ***Vx;
  double ***Vy;
  double ***Vz;

  double rho;
  double scrh;
  double scrh_trace;
  double nu_t, nu_t_prime;

  double dl2;
  double S, Sxx, Syy, Szz, Sxy, Sxz, Syz;
  double divV, dxVx, dxVy, dxVz, dyVx, dyVy, dyVz, dzVx, dzVy, dzVz;
} LES_alpha_ctx;

LES_alpha_ctx *les_alpha_create() {
  LES_alpha_ctx *result;

  result = (LES_alpha_ctx *)malloc(sizeof(LES_alpha_ctx));
  if (result == NULL) {
    return NULL;
  }
  return result;
}

void les_alpha_delete(LES_alpha_ctx *cstr) { free(cstr); }

double les_alpha_get_nut(LES_alpha_ctx *ctx) { return ctx->nu_t; }
double les_alpha_get_nut_trace(LES_alpha_ctx *ctx) { return ctx->nu_t_prime; }

double les_alpha_get_divV(LES_alpha_ctx *ctx) { return ctx->divV; }
double les_alpha_get_Sxx(LES_alpha_ctx *ctx) { return ctx->Sxx; }
double les_alpha_get_Sxy(LES_alpha_ctx *ctx) { return ctx->Sxy; }
double les_alpha_get_Sxz(LES_alpha_ctx *ctx) { return ctx->Sxz; }
double les_alpha_get_Syy(LES_alpha_ctx *ctx) { return ctx->Syy; }
double les_alpha_get_Syz(LES_alpha_ctx *ctx) { return ctx->Syz; }
double les_alpha_get_Szz(LES_alpha_ctx *ctx) { return ctx->Szz; }

double LES_Alpha(LES_alpha_ctx *ctx, int i, int j, int k) {
  return ctx->scrh_trace;
}

double LES_Alpha_xx(LES_alpha_ctx *ctx, int i, int j, int k) {
  return ctx->scrh;
};

double LES_Alpha_xy(LES_alpha_ctx *ctx, int i, int j, int k) {
  return ctx->scrh;
}

double LES_Alpha_xz(LES_alpha_ctx *ctx, int i, int j, int k) {
  return ctx->scrh;
}

double LES_Alpha_yx(LES_alpha_ctx *ctx, int i, int j, int k) {
  return ctx->scrh;
}

double LES_Alpha_yy(LES_alpha_ctx *ctx, int i, int j, int k) {
  return ctx->scrh;
}

double LES_Alpha_yz(LES_alpha_ctx *ctx, int i, int j, int k) {
  return ctx->scrh;
}

double LES_Alpha_zx(LES_alpha_ctx *ctx, int i, int j, int k) {
  return ctx->scrh;
}

double LES_Alpha_zy(LES_alpha_ctx *ctx, int i, int j, int k) {
  return ctx->scrh;
}

double LES_Alpha_zz(LES_alpha_ctx *ctx, int i, int j, int k) {
  return ctx->scrh;
}
#endif
