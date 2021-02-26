#include "gromacs/fhmdlib/new_md_support.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/fhmdlib/new_update.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/fhmdlib/new_tgroup.h"



void compute_FHMD_globals(FILE *fplog, gmx_global_stat *gstat, t_commrec *cr, t_inputrec *ir,
                     t_forcerec *fr, gmx_ekindata_t *ekind,
                     t_state *state, t_mdatoms *mdatoms,
                     t_nrnb *nrnb, t_vcm *vcm, gmx_wallcycle_t wcycle,
                     gmx_enerdata_t *enerd, tensor force_vir, tensor shake_vir, tensor total_vir,
                     tensor pres, rvec mu_tot, gmx_constr_t constr,
                     gmx::SimulationSignaller *signalCoordinator,
                     matrix box, int *totalNumberOfBondedInteractions,
                     gmx_bool *bSumEkinhOld, int flags,gmx_mtop_t *mtop,gmx_int64_t step,gmx_fhmd_global_stat *g_fhmd_s,FHMD *fh)
{
	tensor   corr_vir, corr_pres;
	gmx_bool bEner, bPres, bTemp;
	gmx_bool bStopCM, bGStat,
	         bReadEkin, bEkinAveVel, bScaleEkin, bConstrain;
	real     prescorr, enercorr, dvdlcorr, dvdl_ekin;

	 /*translate CGLO flags to gmx_booleans*/
	bStopCM       = flags & CGLO_STOPCM;
	bGStat        = flags & CGLO_GSTAT;
	bReadEkin     = (flags & CGLO_READEKIN);
	bScaleEkin    = (flags & CGLO_SCALEEKIN);
	bEner         = flags & CGLO_ENERGY;
	bTemp         = flags & CGLO_TEMPERATURE;
	bPres         = (flags & CGLO_PRESSURE);
	bConstrain    = (flags & CGLO_CONSTRAINT);


	 /*we calculate a full state kinetic energy either with full-step velocity verlet
	       or half step where we need the pressure*/

	bEkinAveVel = (ir->eI == eiVV || (ir->eI == eiVVAK && bPres) || bReadEkin);

	 /*in initalization, it sums the shake virial in vv, and to
	   sums ekinh_old in leapfrog (or if we are calculating ekinh_old) for other reasons*/

	 /*########## Kinetic energy  ##############*/

	if (bTemp)
	{
	/*Non-equilibrium MD: this is parallellized, but only does communication
	* when there really is NEMD.*/


		if (PAR(cr) && (ekind->bNEMD))
	    {
			accumulate_u(cr, &(ir->opts), ekind);
	    }
	    if (!bReadEkin)
	    {
	    	calc_FHMD_ke_part(state, &(ir->opts), mdatoms, ekind, nrnb, bEkinAveVel,mtop,fh);
	    }
	}

	 /* Calculate center of mass velocity if necessary, also parallellized */
	 if (bStopCM)
	 {
		 calc_vcm_grp(0, mdatoms->homenr, mdatoms,
	                 state->x, state->v, vcm);
	 }

	 if (bTemp || bStopCM || bPres || bEner || bConstrain)
	 {
		 if (!bGStat)
	     {
	    		/* We will not sum ekinh_old,
	             * so signal that we still have to do it.
	             */
			 *bSumEkinhOld = TRUE;

	     }
	     else
	     {
	    	 gmx::ArrayRef<real> signalBuffer = signalCoordinator->getCommunicationBuffer();
	         if (PAR(cr))
	         {
	            	wallcycle_start(wcycle, ewcMoveE);
	                global_FHMD_stat(gstat, cr, enerd, force_vir, shake_vir, mu_tot,
	                            ir, ekind, constr, bStopCM ? vcm : NULL,
	                            signalBuffer.size(), signalBuffer.data(),
	                            totalNumberOfBondedInteractions,
	                            *bSumEkinhOld, flags, g_fhmd_s,fh);
	                wallcycle_stop(wcycle, ewcMoveE);
	         }
	         signalCoordinator->finalizeSignals();
	         *bSumEkinhOld = FALSE;
	    }
	 }

	 if (!ekind->bNEMD && debug && bTemp && (vcm->nr > 0))
	 {
		 correct_ekin(debug,
	                  0, mdatoms->homenr,
	                  state->v, vcm->group_p[0],
	                  mdatoms->massT, mdatoms->tmass, ekind->ekin);
	 }

	 /* Do center of mass motion removal */
	 if (bStopCM)
	 {
		 check_cm_grp(fplog, vcm, ir, 1);
	     do_stopcm_grp(0, mdatoms->homenr, mdatoms->cVCM,
	                   state->x, state->v, vcm);
	     inc_nrnb(nrnb, eNR_STOPCM, mdatoms->homenr);
	 }

	 if (bEner)
	 {
		 /* Calculate the amplitude of the cosine velocity profile */
		 ekind->cosacc.vcos = ekind->cosacc.mvcos/mdatoms->tmass;
	 }

	 if (bTemp)
	 {
		 /* Sum the kinetic energies of the groups & calc temp */
	     /* compute full step kinetic energies if vv, or if vv-avek and we are computing the pressure with inputrecNptTrotter */
	     /* three maincase:  VV with AveVel (md-vv), vv with AveEkin (md-vv-avek), leap with AveEkin (md).
	     Leap with AveVel is not supported; it's not clear that it will actually work.
	     bEkinAveVel: If TRUE, we simply multiply ekin by ekinscale to get a full step kinetic energy.
	     If FALSE, we average ekinh_old and ekinh*ekinscale_nhc to get an averaged half step kinetic energy.
	     */
	     enerd->term[F_TEMP] = sum_FHMD_ekin(&(ir->opts), ekind, &dvdl_ekin,
	                                        bEkinAveVel, bScaleEkin, step,fh);
	     enerd->dvdl_lin[efptMASS] = (double) dvdl_ekin;

	     enerd->term[F_EKIN] = trace(ekind->ekin);
	}

	 /* ##########  Long range energy information ###### */

	 if (bEner || bPres || bConstrain)
	 {
		 calc_dispcorr(ir, fr, box, state->lambda[efptVDW],
	                   corr_pres, corr_vir, &prescorr, &enercorr, &dvdlcorr);
	 }

	 if (bEner)
	 {
		 enerd->term[F_DISPCORR]  = enercorr;
	     enerd->term[F_EPOT]     += enercorr;
	     enerd->term[F_DVDL_VDW] += dvdlcorr;
	 }

	 /* ########## Now pressure ############## */
	 if (bPres || bConstrain)
	 {

		 m_add(force_vir, shake_vir, total_vir);

	     /* Calculate pressure and apply LR correction if PPPM is used.
	      * Use the box from last timestep since we already called update().
	      */

	     enerd->term[F_PRES] = calc_pres(fr->ePBC, ir->nwall, box, ekind->ekin, total_vir, pres);

	     /* Calculate long range corrections to pressure and energy */
	     /* this adds to enerd->term[F_PRES] and enerd->term[F_ETOT],
	     and computes enerd->term[F_DISPCORR].  Also modifies the
	     total_vir and pres tesors */

	     m_add(total_vir, corr_vir, total_vir);
	     m_add(pres, corr_pres, pres);
	     enerd->term[F_PDISPCORR] = prescorr;
	     enerd->term[F_PRES]     += prescorr;
	}

}
