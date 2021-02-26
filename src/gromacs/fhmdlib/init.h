#ifndef FHMD_INIT_H_
#define FHMD_INIT_H_

int fhmd_init(matrix box, int N_atoms, real mass[], rvec x[], double dt_md, t_grpopts *opts ,gmx_mtop_t *mtop, t_commrec *cr,gmx_fhmd_global_stat *g_fhmd_s,FHMD *fh);

#endif /* FHMD_INIT_H_ */
