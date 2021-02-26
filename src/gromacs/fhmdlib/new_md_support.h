#ifndef FHMD_NEW_MD_SUPPORT_H
#define FHMD_NEW_MD_SUPPORT_H

#include "gromacs/fhmdlib/data_structures.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/timing/wallcycle.h"
//#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/constr.h"


struct gmx_constr;
struct gmx_ekindata_t;
struct gmx_enerdata_t;
struct gmx_global_stat;
struct gmx_multisim_t;
struct gmx_signalling_t;
struct t_extmass;
struct t_forcerec;
struct t_grpopts;
struct t_lambda;
struct t_nrnb;
struct t_state;
struct t_trxframe;

/*same as #include "gromacs/mdlib/simulationsignal.h" */

namespace gmx
{

class SimulationSignaller;

}

/* Define a number of flags to better control the information
 * passed to compute_globals in md.c and global_stat.
 */

/* we are computing the kinetic energy from average velocities */
#define CGLO_EKINAVEVEL     (1<<2)
/* we are removing the center of mass momenta */
#define CGLO_STOPCM         (1<<3)
/* bGStat is defined in do_md */
#define CGLO_GSTAT          (1<<4)
/* Sum the energy terms in global computation */
#define CGLO_ENERGY         (1<<6)
/* Sum the kinetic energy terms in global computation */
#define CGLO_TEMPERATURE    (1<<7)
/* Sum the kinetic energy terms in global computation */
#define CGLO_PRESSURE       (1<<8)
/* Sum the constraint term in global computation */
#define CGLO_CONSTRAINT     (1<<9)
/* Reading ekin from the trajectory */
#define CGLO_READEKIN       (1<<10)
/* we need to reset the ekin rescaling factor here */
#define CGLO_SCALEEKIN      (1<<11)
/* After a new DD partitioning, we need to set a flag to schedule
 * global reduction of the total number of bonded interactions that
 * will be computed, to check none are missing. */
#define CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS (1<<12)


/* Compute global variables during integration for FHMD*/

void compute_FHMD_globals(FILE *fplog, gmx_global_stat *gstat, t_commrec *cr, t_inputrec *ir,
		                   t_forcerec *fr, gmx_ekindata_t *ekind,
						   t_state *state, t_mdatoms *mdatoms,
						   t_nrnb *nrnb, t_vcm *vcm, gmx_wallcycle_t wcycle,
						   gmx_enerdata_t *enerd, tensor force_vir, tensor shake_vir, tensor total_vir,
						   tensor pres, rvec mu_tot, gmx_constr *constr,
						   gmx::SimulationSignaller *signalCoordinator,
						   matrix box, int *totalNumberOfBondedInteractions,
						   gmx_bool *bSumEkinhOld, int flags, gmx_mtop_t *mtop,gmx_int64_t step,gmx_fhmd_global_stat *g_fhmd_s,FHMD *fh);

#endif
