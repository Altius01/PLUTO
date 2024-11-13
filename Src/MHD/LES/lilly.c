/* ------------------------------------------------------------
    Main analysis loop.
    Compute volume-integrated magnetic pressure and stresses
   ------------------------------------------------------------ */

vy2_av = vy_max = vy_min = 0.0;
DOM_LOOP(k, j, i) {
  dV = dx[i] * dy[j] * dz[k];
  vy2_av += d->Vc[VX2][k][j][i] * d->Vc[VX2][k][j][i] * dV;
  vy_max = MAX(d->Vc[VX2][k][j][i], vy_max);
  vy_min = MIN(d->Vc[VX2][k][j][i], vy_min);
}

/* -- divide by total volume -- */

vy2_av /= tot_vol;

/* -- parallel reduce -- */

#ifdef PARALLEL
MPI_Allreduce(&vy2_av, &scrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
vy2_av = scrh;

MPI_Allreduce(&vy_max, &scrh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
vy_max = scrh;

MPI_Allreduce(&vy_min, &scrh, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
vy_min = scrh;
MPI_Barrier(MPI_COMM_WORLD);
#endif
