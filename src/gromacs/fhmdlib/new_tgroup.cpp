#include "gromacs/fhmdlib/new_tgroup.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/fhmdlib/data_structures.h"

real sum_FHMD_ekin(t_grpopts *opts, gmx_ekindata_t *ekind, real *dekindlambda,
              gmx_bool bEkinAveVel, gmx_bool bScaleEkin, gmx_int64_t step,FHMD *fh)
{
    int           i, j, m, ngtc,nzbin ;
    double        T_bin, nrdf_bin ;
    t_grp_tcstat *tcstat;
    AFM_KE       *ke_data;
    real          T, nrdf, nd, *ndf;
    double        nd_bin;

    nzbin = fh->nzbin;
    ke_data = fh->ke_data;
    tcstat = tcstat  = ekind->tcstat;

    T_bin = 0.0;
    nrdf_bin = 0.0;
    ngtc = opts->ngtc;
    ndf  = opts->nrdf;
    T    = 0;
    nrdf = 0;

    clear_mat(fh->ekin);

    /*cal temp for each bin*/

    for(i = 0; i < nzbin; i++)
    {
    	nd_bin = ke_data[i].nrdf;

    	if (nd_bin>0)
    	{
    		if (bEkinAveVel)
    		{
    			if (!bScaleEkin)
    			{
    			 /* in this case, kinetic energy is from the current velocities already */
    			 msmul(tcstat->ekinf, tcstat->ekinscalef_nhc, tcstat->ekinf);
    			 }
    		}

    		else
    		{
    			for (j = 0; j < DIM; j++)
    			{
    				for (m = 0; m < DIM; m++)
    				{
    					ke_data[i].ekinf[j][m] =
    							0.5*(ke_data[i].ekinh[j][m]*tcstat->ekinscaleh_nhc + ke_data[i].ekinh_old[j][m]);
    				}
    			}
    		}

    		m_add(ke_data[i].ekinf,fh->ekin,fh->ekin);

    		//printf("the bin number  = %d	KE = %g		nd_bin = %g \n", i,trace(ke_data[i].ekinf),nd_bin);

    		ke_data[i].Th = calc_temp(trace(ke_data[i].ekinh),nd_bin);
    		ke_data[i].T  = calc_temp(trace(ke_data[i].ekinf),nd_bin);

    	    /* after the scaling factors have been multiplied in, we can remove them */
    		if (bEkinAveVel)
    		{
    			tcstat->ekinscalef_nhc = 1.0;
    		}
    		else
    		{
    			tcstat->ekinscaleh_nhc = 1.0;
    		}

    	}
    	else
    	{
    		ke_data[i].T  = 0;
    		ke_data[i].Th = 0;
    	}

    	//printf("the bin number  = %d	KE = %g		nd_bin = %g \n", i,trace(ke_data[i].ekinf),nd_bin);

    	T_bin    += nd_bin*ke_data[i].T;
    	nrdf_bin += nd_bin;

    	//printf("the bin number is %d and the temp is %g \n", i,ke_data[i].T);
    	//printf("the bin number = %d T_bin = %g\n", i,T_bin);
    }

    if(nrdf_bin > 0)
    {
    	T_bin /=nrdf_bin;
    }

    /*original Gromacs code for temp calc*/

    clear_mat(ekind->ekin);

    for (i = 0; (i < ngtc); i++)
    {

        nd     = ndf[i];
        tcstat = &ekind->tcstat[i];
        /* Sometimes a group does not have degrees of freedom, e.g.
         * when it consists of shells and virtual sites, then we just
         * set the temperatue to 0 and also neglect the kinetic
         * energy, which should be  zero anyway.
         */

        if (nd > 0)
        {
            if (bEkinAveVel)
            {
                if (!bScaleEkin)
                {
                    /* in this case, kinetic energy is from the current velocities already */
                    msmul(tcstat->ekinf, tcstat->ekinscalef_nhc, tcstat->ekinf);
                }
            }
            else
            {
                /* Calculate the full step Ekin as the average of the half steps */
                for (j = 0; (j < DIM); j++)
                {
                    for (m = 0; (m < DIM); m++)
                    {
                        tcstat->ekinf[j][m] =
                            0.5*(tcstat->ekinh[j][m]*tcstat->ekinscaleh_nhc + tcstat->ekinh_old[j][m]);
                    }
                }
            }
            m_add(tcstat->ekinf, ekind->ekin, ekind->ekin);

            tcstat->Th = calc_temp(trace(tcstat->ekinh), nd);
            tcstat->T  = calc_temp(trace(tcstat->ekinf), nd);

            /* after the scaling factors have been multiplied in, we can remove them */
            if (bEkinAveVel)
            {
                tcstat->ekinscalef_nhc = 1.0;
            }
            else
            {
                tcstat->ekinscaleh_nhc = 1.0;
            }
        }
        else
        {
            tcstat->T  = 0;
            tcstat->Th = 0;
        }
        T    += nd*tcstat->T;
        nrdf += nd;
    }
    if (nrdf > 0)
    {
        T /= nrdf;
    }
    if (dekindlambda)
    {
        if (bEkinAveVel)
        {
            *dekindlambda = ekind->dekindl;
        }
        else
        {
            *dekindlambda = 0.5*(ekind->dekindl + ekind->dekindl_old);
        }
    }

    /*for debug only*/
    real ekind_ekin,fh_ekin;
    ekind_ekin = trace(ekind->ekin);
    fh_ekin = trace(fh->ekin);

    if (step !=0)
	{
    	if(fabs(ekind_ekin - fh_ekin) > 1)
    	{
    		//printf(MAKE_RED "FHMD eror:Overall KE calculated by original gromacs and different bins is not equal\n");
    		//printf("the KE by gromacs is %g			the KE by bin is %g\n",ekind_ekin,fh_ekin);
    		//exit(2);
    	}

    	if(fabs(T_bin - T) > 1)
    	{
    		//printf(MAKE_RED "FHMD error:The temperature calculated by original gromacs and different bins is not equal,probably the number of degree of freedom of AFM is not correct\n");
    		//printf("the temperature by gromacs is %g		the temperature by bin is %g\n",T,T_bin);
    	    //exit(3);
    	}
	}

    /*for debug only*/
    return T;
}
