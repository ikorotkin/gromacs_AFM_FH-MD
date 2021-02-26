#ifndef FHMD_NEW_UPDATE_H_
#define FHMD_NEW_UPDATE_H_

#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/constr.h"

void fhmd_update_coords(FILE        *fplog,
                   gmx_int64_t       step,
                   t_inputrec       *inputrec,  /* input record and box stuff   */
                   t_mdatoms        *md,
                   t_state          *state,
                   rvec             *f,    /* forces on home particles */
                   t_fcdata         *fcd,
                   gmx_ekindata_t   *ekind,
                   matrix            M,
                   gmx_update_t     *upd,
                   int               UpdatePart,
                   t_commrec        *cr, /* these shouldn't be here -- need to think about it */
                   gmx_constr_t      constr,
				   gmx_mtop_t       *mtop,
                   FHMD             *fh);

void fhmd_do_update_md(int start, int nrend,
                       double dt, int nstpcouple,
                       t_grp_tcstat *tcstat,
                       double nh_vxi[],
                       gmx_bool bNEMD, t_grp_acc *gstat, rvec accel[],
                       ivec nFreeze[],
                       real invmass[],
                       unsigned short ptype[], unsigned short cFREEZE[],
                       unsigned short cACC[], unsigned short cTC[],
                       rvec x[], rvec xprime[], rvec v[],
                       rvec f[], matrix M,
                       gmx_bool bNH, gmx_bool bPR, t_commrec *cr, gmx_mtop_t *mtop,FHMD *fh);

void update_FH_MD_tcouple (gmx_int64_t       step,
						   t_inputrec       *inputrec,
						   t_state          *state,
						   gmx_ekindata_t   *ekind,
						   t_extmass        *MassQ,
						   t_mdatoms        *md,
						   FHMD             *fh
						  );

void berendsen_FH_MD_tcoupl (t_inputrec *ir, gmx_ekindata_t *ekind, real dt,gmx_int64_t step,FHMD *fh);

void calc_FHMD_ke_part(t_state *state, t_grpopts *opts, t_mdatoms *md,
                  	   gmx_ekindata_t *ekind, t_nrnb *nrnb, gmx_bool bEkinAveVel,gmx_mtop_t *mtop,FHMD *fh);

#endif /* FHMD_NEW_UPDATE_H_ */
