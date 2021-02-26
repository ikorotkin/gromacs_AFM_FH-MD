#include "gmxpre.h"

#include <stdio.h>
#include <string.h>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/rbin.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

/*
 * FHMD includes
 */
//#include "gromacs/fhmdlib/data_structures.h"


typedef struct gmx_global_stat
{
    t_bin *rb;
    int   *itc0;
    int   *itc1;
} t_gmx_global_stat;

static int filter_enerdterm(real *afrom, gmx_bool bToBuffer, real *ato,
                            gmx_bool bTemp, gmx_bool bPres, gmx_bool bEner)
{
    int i, to, from;

    from = 0;
    to   = 0;
    for (i = 0; i < F_NRE; i++)
    {
        if (bToBuffer)
        {
            from = i;
        }
        else
        {
            to = i;
        }
        switch (i)
        {
            case F_EKIN:
            case F_TEMP:
            case F_DKDL:
                if (bTemp)
                {
                    ato[to++] = afrom[from++];
                }
                break;
            case F_PRES:
            case F_PDISPCORR:
                if (bPres)
                {
                    ato[to++] = afrom[from++];
                }
                break;
            default:
                if (bEner)
                {
                    ato[to++] = afrom[from++];
                }
                break;
        }
    }

    return to;
}

void global_FHMD_stat(gmx_global_stat_t gs,
                 t_commrec *cr, gmx_enerdata_t *enerd,
                 tensor fvir, tensor svir, rvec mu_tot,
                 t_inputrec *inputrec,
                 gmx_ekindata_t *ekind, gmx_constr_t constr,
                 t_vcm *vcm,
                 int nsig, real *sig,
                 int *totalNumberOfBondedInteractions,
                 gmx_bool bSumEkinhOld, int flags, gmx_fhmd_global_stat *g_fhmd_s, FHMD *fh)
/* instead of current system, gmx_booleans for summing virial, kinetic energy, and other terms */
{
    t_bin     *rb,*rb_bin;
    int       *itc0, *itc1, *itc0_bin, *itc1_bin;
    int        ie    = 0, ifv = 0, isv = 0, irmsd = 0, imu = 0;
    int        idedl = 0, idedlo = 0, idvdll = 0, idvdlnl = 0, iepl = 0, icm = 0, imass = 0, ica = 0, inb = 0;
    int        isig  = -1;
    int        icj   = -1, ici = -1, icx = -1;
    int        inn[egNR];
    real       copyenerd[F_NRE];
    int        nener, j, nener_bin;
    real      *rmsd_data = NULL;
    double     nb;
    gmx_bool   bVV, bTemp, bEner, bPres, bConstrVir, bEkinAveVel, bReadEkin;
    bool       checkNumberOfBondedInteractions = flags & CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS;

    //double     *nrdf;
    AFM_KE     *ke_data;

    bVV           = EI_VV(inputrec->eI);
    bTemp         = flags & CGLO_TEMPERATURE;
    bEner         = flags & CGLO_ENERGY;
    bPres         = (flags & CGLO_PRESSURE);
    bConstrVir    = (flags & CGLO_CONSTRAINT);
    bEkinAveVel   = (inputrec->eI == eiVV || (inputrec->eI == eiVVAK && bPres));
    bReadEkin     = (flags & CGLO_READEKIN);

    rb   = gs->rb;
    itc0 = gs->itc0;
    itc1 = gs->itc1;
    ke_data = fh->ke_data;

    reset_bin(rb);
    /* This routine copies all the data to be summed to one big buffer
     * using the t_bin struct.
     */

    /* First, we neeed to identify which enerd->term should be
       communicated.  Temperature and pressure terms should only be
       communicated and summed when they need to be, to avoid repeating
       the sums and overcounting. */

    nener = filter_enerdterm(enerd->term, TRUE, copyenerd, bTemp, bPres, bEner);

    /* First, the data that needs to be communicated with velocity verlet every time
       This is just the constraint virial.*/
    if (bConstrVir)
    {
        isv = add_binr(rb, DIM*DIM, svir[0]);
        where();
    }

/* We need the force virial and the kinetic energy for the first time through with velocity verlet */
    if (bTemp || !bVV)
    {
        if (ekind)
        {
            for (j = 0; (j < inputrec->opts.ngtc); j++)
            {
                if (bSumEkinhOld)
                {
                    itc0[j] = add_binr(rb, DIM*DIM, ekind->tcstat[j].ekinh_old[0]);
                }
                if (bEkinAveVel && !bReadEkin)
                {
                    itc1[j] = add_binr(rb, DIM*DIM, ekind->tcstat[j].ekinf[0]);
                }
                else if (!bReadEkin)
                {
                    itc1[j] = add_binr(rb, DIM*DIM, ekind->tcstat[j].ekinh[0]);
                }
            }
            /* these probably need to be put into one of these categories */
            where();
            idedl = add_binr(rb, 1, &(ekind->dekindl));
            if (bSumEkinhOld)
            {
                idedlo = add_binr(rb, 1, &(ekind->dekindl_old));
            }
            where();
            ica   = add_binr(rb, 1, &(ekind->cosacc.mvcos));
            where();
        }
    }
    where();

    if (bPres)
    {
        ifv = add_binr(rb, DIM*DIM, fvir[0]);
    }


    if (bEner)
    {
        where();
        ie  = add_binr(rb, nener, copyenerd);
        where();
        if (constr)
        {
            rmsd_data = constr_rmsd_data(constr);
            if (rmsd_data)
            {
                irmsd = add_binr(rb, 2, rmsd_data);
            }
        }
        if (!inputrecNeedMutot(inputrec))
        {
            imu = add_binr(rb, DIM, mu_tot);
            where();
        }

        for (j = 0; (j < egNR); j++)
        {
            inn[j] = add_binr(rb, enerd->grpp.nener, enerd->grpp.ener[j]);
        }
        where();
        if (inputrec->efep != efepNO)
        {
            idvdll  = add_bind(rb, efptNR, enerd->dvdl_lin);
            idvdlnl = add_bind(rb, efptNR, enerd->dvdl_nonlin);
            if (enerd->n_lambda > 0)
            {
                iepl = add_bind(rb, enerd->n_lambda, enerd->enerpart_lambda);
            }
        }
    }

    if (vcm)
    {
        icm   = add_binr(rb, DIM*vcm->nr, vcm->group_p[0]);
        where();
        imass = add_binr(rb, vcm->nr, vcm->group_mass);
        where();
        if (vcm->mode == ecmANGULAR)
        {
            icj   = add_binr(rb, DIM*vcm->nr, vcm->group_j[0]);
            where();
            icx   = add_binr(rb, DIM*vcm->nr, vcm->group_x[0]);
            where();
            ici   = add_binr(rb, DIM*DIM*vcm->nr, vcm->group_i[0][0]);
            where();
        }
    }

    if (checkNumberOfBondedInteractions)
    {
        nb  = cr->dd->nbonded_local;
        inb = add_bind(rb, 1, &nb);
    }
    where();
    if (nsig > 0)
    {
        isig = add_binr(rb, nsig, sig);
    }

    /* Global sum it all */
    if (debug)
    {
        fprintf(debug, "Summing %d energies\n", rb->maxreal);
    }
    sum_bin(rb, cr);
    where();

    /* Extract all the data locally */

    if (bConstrVir)
    {
        extract_binr(rb, isv, DIM*DIM, svir[0]);
    }

    /* We need the force virial and the kinetic energy for the first time through with velocity verlet */
    if (bTemp || !bVV)
    {
        if (ekind)
        {
            for (j = 0; (j < inputrec->opts.ngtc); j++)
            {
                if (bSumEkinhOld)
                {
                    extract_binr(rb, itc0[j], DIM*DIM, ekind->tcstat[j].ekinh_old[0]);
                }
                if (bEkinAveVel && !bReadEkin)
                {
                    extract_binr(rb, itc1[j], DIM*DIM, ekind->tcstat[j].ekinf[0]);
                }
                else if (!bReadEkin)
                {
                    extract_binr(rb, itc1[j], DIM*DIM, ekind->tcstat[j].ekinh[0]);
                }
            }
            extract_binr(rb, idedl, 1, &(ekind->dekindl));
            if (bSumEkinhOld)
            {
                extract_binr(rb, idedlo, 1, &(ekind->dekindl_old));
            }
            extract_binr(rb, ica, 1, &(ekind->cosacc.mvcos));
            where();
        }
    }
    if (bPres)
    {
        extract_binr(rb, ifv, DIM*DIM, fvir[0]);
    }

    if (bEner)
    {
        extract_binr(rb, ie, nener, copyenerd);
        if (rmsd_data)
        {
            extract_binr(rb, irmsd, 2, rmsd_data);
        }
        if (!inputrecNeedMutot(inputrec))
        {
            extract_binr(rb, imu, DIM, mu_tot);
        }

        for (j = 0; (j < egNR); j++)
        {
            extract_binr(rb, inn[j], enerd->grpp.nener, enerd->grpp.ener[j]);
        }
        if (inputrec->efep != efepNO)
        {
            extract_bind(rb, idvdll, efptNR, enerd->dvdl_lin);
            extract_bind(rb, idvdlnl, efptNR, enerd->dvdl_nonlin);
            if (enerd->n_lambda > 0)
            {
                extract_bind(rb, iepl, enerd->n_lambda, enerd->enerpart_lambda);
            }
        }
        where();

        filter_enerdterm(copyenerd, FALSE, enerd->term, bTemp, bPres, bEner);
    }

    if (vcm)
    {
        extract_binr(rb, icm, DIM*vcm->nr, vcm->group_p[0]);
        where();
        extract_binr(rb, imass, vcm->nr, vcm->group_mass);
        where();
        if (vcm->mode == ecmANGULAR)
        {
            extract_binr(rb, icj, DIM*vcm->nr, vcm->group_j[0]);
            where();
            extract_binr(rb, icx, DIM*vcm->nr, vcm->group_x[0]);
            where();
            extract_binr(rb, ici, DIM*DIM*vcm->nr, vcm->group_i[0][0]);
            where();
        }
    }

    if (checkNumberOfBondedInteractions)
    {
        extract_bind(rb, inb, 1, &nb);
        *totalNumberOfBondedInteractions = static_cast<int>(nb+0.5);
    }

    if (nsig > 0)
    {
        extract_binr(rb, isig, nsig, sig);
    }
    where();


    /*for the AFM kinetic enerty term*/

    rb_bin = g_fhmd_s->rb_bin;
    itc0_bin = g_fhmd_s->itc0_bin;
    itc1_bin = g_fhmd_s->itc1_bin;



    reset_bin(rb_bin);

    /*maybe not needed*/
    nener_bin = filter_enerdterm(enerd->term, TRUE, copyenerd, bTemp, bPres, bEner);

     /*We need the force virial and the kinetic energy for the first time through with velocity verlet*/
        if (bTemp || !bVV)
        {
            if (ekind) //we add the KE data in the ke_data as the gromacs add the KE data in ekind
            {
                for (j = 0; (j < fh->nzbin); j++)
                {
                    if (bSumEkinhOld)
                    {
                        itc0_bin[j] = add_binr(rb_bin, DIM*DIM, ke_data[j].ekinh_old[0]);
                    }
                    if (bEkinAveVel && !bReadEkin)
                    {
                        itc1_bin[j] = add_binr(rb_bin, DIM*DIM, ke_data[j].ekinf[0]);
                    }
                    else if (!bReadEkin)
                    {
                        itc1_bin[j] = add_binr(rb_bin, DIM*DIM, ke_data[j].ekinh[0]);
                    }
                }
                 /*these probably need to be put into one of these categories*/
                where();
            }
        }
        where();

        /* Global sum it all */
        if (debug)
        {
        	fprintf(debug, "Summing %d energies\n", rb->maxreal);
        }
        sum_bin(rb_bin, cr);
        where();

        /* Extract all the data locally */

        /* We need the force virial and the kinetic energy for the first time through with velocity verlet */
        if (bTemp || !bVV)
        {
        	if (ekind)
            {
        		for (j = 0; (j < fh->nzbin); j++)
                {
        			if (bSumEkinhOld)
                    {
        				extract_binr(rb_bin, itc0_bin[j], DIM*DIM, ke_data[j].ekinh_old[0]);
                    }
                    if (bEkinAveVel && !bReadEkin)
                    {
                    	extract_binr(rb_bin, itc1_bin[j], DIM*DIM, ke_data[j].ekinf[0]);
                    }
                    else if (!bReadEkin)
                    {
                    	extract_binr(rb_bin, itc1_bin[j], DIM*DIM, ke_data[j].ekinh[0]);
                    }

                }
                where();
            }
       }

       //pass the number degree of freedom in ke_data to nrdf

      for (j = 0; (j < fh->nzbin); j++)
       {
    	   fh->nrdf_mpi[j] = ke_data[j].nrdf;
       }

       gmx_sumd(fh->nzbin, fh->nrdf_mpi, cr);

       //extract the data
       for (j = 0; (j < fh->nzbin); j++)
       {
    	   ke_data[j].nrdf = fh->nrdf_mpi[j];
       }

}
