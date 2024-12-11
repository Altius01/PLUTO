#ifndef LES_PHI_H
#define LES_PHI_H

#include "pluto.h"

struct LES_phi_ctx;

typedef struct LES_phi_ctx LES_phi_ctx;

LES_phi_ctx *les_phi_create();

void les_phi_delete(LES_phi_ctx *cstr);

void init_les_phi_ctx(LES_phi_ctx *ctx, const Data *d, Grid *grid);

void update_les_phi_ctx(LES_phi_ctx *ctx, double *vi, Grid *grid, int i, int j,
                        int k);

double les_phi_get_eta(LES_phi_ctx *ctx);

double les_phi_get_Jxx(LES_phi_ctx *ctx);
double les_phi_get_Jxy(LES_phi_ctx *ctx);
double les_phi_get_Jxz(LES_phi_ctx *ctx);
double les_phi_get_Jyy(LES_phi_ctx *ctx);
double les_phi_get_Jyz(LES_phi_ctx *ctx);
double les_phi_get_Jzz(LES_phi_ctx *ctx);
double les_phi_get_Jyx(LES_phi_ctx *ctx);
double les_phi_get_Jzx(LES_phi_ctx *ctx);
double les_phi_get_Jzy(LES_phi_ctx *ctx);

double LES_Phi_xx(LES_phi_ctx *ctx, int i, int j, int k);
double LES_Phi_xy(LES_phi_ctx *ctx, int i, int j, int k);
double LES_Phi_xz(LES_phi_ctx *ctx, int i, int j, int k);
double LES_Phi_yy(LES_phi_ctx *ctx, int i, int j, int k);
double LES_Phi_yz(LES_phi_ctx *ctx, int i, int j, int k);
double LES_Phi_zx(LES_phi_ctx *ctx, int i, int j, int k);
double LES_Phi_zz(LES_phi_ctx *ctx, int i, int j, int k);
double LES_Phi_zy(LES_phi_ctx *ctx, int i, int j, int k);
double LES_Phi_yx(LES_phi_ctx *ctx, int i, int j, int k);
#endif
