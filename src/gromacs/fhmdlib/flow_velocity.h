#ifndef FLOW_VELOCITY_H
#define FLOW_VELOCITY_H

#include "data_structures.h"

/*for AFM*/
int read_flow_velocity (FHMD *fh, char const *fname_in_velocity);
void skip_coordinate(FILE *fvel);
void compute_tem(FHMD *fh,const rvec x[], rvec v[], real mass[], int N_atoms, t_commrec *cr, gmx_mtop_t *mtop);

/*for protein diffused in shear flow*/
void compute_vx_instant_mpi (FHMD *fh, rvec x[], rvec v[], real mass[], int N_atoms,t_commrec *cr);
void compute_average_vx_instant (FHMD *fh);


int find_ind(int n, const rvec x[],FHMD *fh);

#endif
