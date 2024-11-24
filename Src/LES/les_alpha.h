#ifndef LES_ALPHA_H
#define LES_ALPHA_H

#include "pluto.h"

struct LES_alpha_ctx;

typedef struct LES_alpha_ctx LES_alpha_ctx;

LES_alpha_ctx *les_alpha_create();

void les_alpha_delete(LES_alpha_ctx *cstr);

void init_les_alpha_ctx(LES_alpha_ctx *ctx, const Data *d, Grid *grid);

void update_les_alpha_ctx(LES_alpha_ctx *ctx, double *vi, Grid *grid, int i,
                          int j, int k);

double les_alpha_get_nut(LES_alpha_ctx *ctx);
double les_alpha_get_nut_trace(LES_alpha_ctx *ctx);

double les_alpha_get_divV(LES_alpha_ctx *ctx);
double les_alpha_get_Sxx(LES_alpha_ctx *ctx);
double les_alpha_get_Sxy(LES_alpha_ctx *ctx);
double les_alpha_get_Sxz(LES_alpha_ctx *ctx);
double les_alpha_get_Syy(LES_alpha_ctx *ctx);
double les_alpha_get_Syz(LES_alpha_ctx *ctx);
double les_alpha_get_Szz(LES_alpha_ctx *ctx);

double LES_Alpha(LES_alpha_ctx *ctx, int i, int j, int k);
double LES_Alpha_xx(LES_alpha_ctx *ctx, int i, int j, int k);
double LES_Alpha_xy(LES_alpha_ctx *ctx, int i, int j, int k);
double LES_Alpha_xz(LES_alpha_ctx *ctx, int i, int j, int k);
double LES_Alpha_yy(LES_alpha_ctx *ctx, int i, int j, int k);
double LES_Alpha_yz(LES_alpha_ctx *ctx, int i, int j, int k);
double LES_Alpha_zx(LES_alpha_ctx *ctx, int i, int j, int k);
double LES_Alpha_zz(LES_alpha_ctx *ctx, int i, int j, int k);
double LES_Alpha_zy(LES_alpha_ctx *ctx, int i, int j, int k);
double LES_Alpha_yx(LES_alpha_ctx *ctx, int i, int j, int k);
#endif
