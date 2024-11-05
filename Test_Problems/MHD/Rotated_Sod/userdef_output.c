#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 ***************************************************************** */
{
  int i, j, k;  
  double ***flag;

  flag = GetUserVar("flag");
  DOM_LOOP(k,j,i){
    flag[k][j][i] = d->flag[k][j][i];
  }

}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

#if PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}





