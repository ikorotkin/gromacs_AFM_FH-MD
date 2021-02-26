#include "gromacs/mdtypes/group.h"
#include "gromacs/topology/atoms.h"
#include "data_structures.h"
#include "interpolation.h"
#include "sfunction.h"
#include "macro.h"

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/fhmdlib/new_update.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/topology/mtop_util.h"     /* This is for gmx_mtop_atominfo_global() */
#include "flow_velocity.h"

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
                       gmx_bool bNH, gmx_bool bPR, t_commrec *cr, gmx_mtop_t *mtop,FHMD *fh)
{
    double imass, w_dt;
    int    gf = 0, ga = 0, gt = 0;
    rvec   vrel;
    real   vn, vv, va, vb, vnrel;
    real   lg, vxi = 0, u;
    int    n, d;

    /* FHMD variables */
    FH_arrays   *arr = fh->arr;
    AFM_KE      *ke_data = fh->ke_data;
    int          ind,ibinz;
    double       invro_dt;
    double       S = fh->S;
    double       gamma_u, gamma_x;
    int          nbr[8];
    dvec         xi;
    dvec         f_fh, u_fh, alpha_term, beta_term, grad_ro;
    dvec         u_fh_flow;
    const double g_eps = 1e-10;

    char  *atomname, *resname;
    int res_nr;
    //double v_d[nrend];
    //double vx[fh->nzbin], mvx[fh->nzbin];
    //int count ;
    const double Vm = 0.05;
    const double ax = 2.0*Vm/fh->box[ZZ];
    const double bx = -Vm;

    if (bNH || bPR)
    {
        /* Update with coupling to extended ensembles, used for
         * Nose-Hoover and Parrinello-Rahman coupling
         * Nose-Hoover uses the reversible leap-frog integrator from
         * Holian et al. Phys Rev E 52(3) : 2338, 1995
         */

        /* FHMD Error */
        printf(MAKE_RED "\nFHMD: ERROR: FH-MD coupling doesn't support Nose-Hoover and Parrinello-Rahman\n" RESET_COLOR "\n");
        exit(11);

    }
    else if (cFREEZE != NULL ||
             nFreeze[0][XX] || nFreeze[0][YY] || nFreeze[0][ZZ] ||
             bNEMD)
    {
        /* Update with Berendsen/v-rescale coupling and freeze or NEMD */

        /* FHMD Error */
    	if (fh->thermostat_function == bin)
		{
    		printf(MAKE_RED "\nFHMD: ERROR: bin thermostat is not supported with freezed group in mdp file \n" RESET_COLOR "\n");
        	//exit(12);
		}

    	 /* Update with Berendsen/v-rescale coupling and freeze or NEMD */
    	 for (n = start; n < nrend; n++)
    	 {
    		 w_dt = invmass[n]*dt;
    	     if (cFREEZE)
    	     {
    	    	 gf   = cFREEZE[n];
    	     }
    	     if (cACC)
    	     {
    	    	 ga   = cACC[n];
    	     }
    	     if (cTC)
    	     {
    	    	 gt   = cTC[n];
    	     }
    	     lg   = tcstat[gt].lambda;

    	     for (d = 0; d < DIM; d++)
    	     {
    	    	 vn             = v[n][d];
    	         if ((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d])
    	         {
    	             //w_dt     = invmass[n]*dt;
    	             ind      = fh->ind[n];
    	             invro_dt = arr[ind].inv_ro*dt;

    	             trilinear_find_neighbours(x[n], n, xi, nbr, fh);

    	             if(fh->scheme == Two_Way)
    	             trilinear_interpolation(f_fh,   xi, INTERPOLATE(f_fh));
    	             else
    	             clear_dvec(f_fh);
    	             trilinear_interpolation(u_fh,       xi, INTERPOLATE(u_fh));
    	             trilinear_interpolation(alpha_term, xi, INTERPOLATE(alpha_term));
    	             trilinear_interpolation(beta_term,  xi, INTERPOLATE(beta_term));
    	             trilinear_interpolation(grad_ro,    xi, INTERPOLATE(grad_ro));

    	             #ifdef FHMD_DEBUG_INTERPOL
    	                    if(!(n % 10000) && !(fh->step_MD % 100))
    	                    printf("\nStep %d, atom #%d (%g %g %g): %g %g %g\nneighbour cells: %d %d %d %d %d %d %d %d\nrescaled coordinates: %g %g %g\n",
    	                    fh->step_MD, n, x[n][0], x[n][1], x[n][2], u_fh[0], u_fh[1], u_fh[2],
    	                    nbr[0], nbr[1], nbr[2], nbr[3], nbr[4], nbr[5], nbr[6], nbr[7], xi[0], xi[1], xi[2]);
    	             #endif

    	             if(fh->S_function == moving_sphere)
    	             //S = fhmd_Sxyz_r(x[n], fh->protein_com, fh);     // MD/FH sphere follows protein
    	               S = fhmd_Sxyz_protein(x[n], fh->prot_com,fh);
    	             else if(fh->S_function == fixed_sphere)
    	             S = fhmd_Sxyz_r(x[n], fh->box05, fh);           // Fixed MD/FH sphere
    	             else if(fh->S_function == AFM)
    	             S = fhmd_Sxyz_AFM(x[n], fh->box05, fh);
    	             if(fh->scheme == One_Way)
    	             {
    	            	 vv           = lg*vn + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + alpha_term[d] + S*pow((1 - S),fh->gamma)*beta_term[d])*invro_dt;
    	                 v[n][d]      = vv;
    	                 xprime[n][d] = x[n][d] + ((1 - S)*vv + S*u_fh[d])*dt + S*pow((1 - S),fh->gamma)*grad_ro[d]*invro_dt;

    	             }
    	             else
    	             {
    	            	 printf(MAKE_RED "the freezed group in mdp file is oly supported by one-way coupling\n" RESET_COLOR "\n");
    	            	 exit (1);
    	             }
    	         }
    	         else
    	         {
    	        	 v[n][d]        = 0.0;
    	             xprime[n][d]   = x[n][d];
    	         }
    	     }
    	}

    }

    else
    {
        /* Plain update with Berendsen/v-rescale coupling */
        for (n = start; n < nrend; n++)
        {
            if ((ptype[n] != eptVSite) && (ptype[n] != eptShell))
            {
                w_dt     = invmass[n]*dt;
                ind      = fh->ind[n];
                invro_dt = arr[ind].inv_ro*dt;

                trilinear_find_neighbours(x[n], n, xi, nbr, fh);

                if(fh->scheme == Two_Way)
                    trilinear_interpolation(f_fh,   xi, INTERPOLATE(f_fh));
                else
                    clear_dvec(f_fh);
                trilinear_interpolation(u_fh,       xi, INTERPOLATE(u_fh));
                trilinear_interpolation(alpha_term, xi, INTERPOLATE(alpha_term));
                trilinear_interpolation(beta_term,  xi, INTERPOLATE(beta_term));
                trilinear_interpolation(grad_ro,    xi, INTERPOLATE(grad_ro));
                trilinear_interpolation(u_fh_flow,  xi, INTERPOLATE(u_fh_flow));

#ifdef FHMD_DEBUG_INTERPOL
                if(!(n % 10000) && !(fh->step_MD % 100))
                    printf("\nStep %d, atom #%d (%g %g %g): %g %g %g\nneighbour cells: %d %d %d %d %d %d %d %d\nrescaled coordinates: %g %g %g\n",
                            fh->step_MD, n, x[n][0], x[n][1], x[n][2], u_fh[0], u_fh[1], u_fh[2],
                            nbr[0], nbr[1], nbr[2], nbr[3], nbr[4], nbr[5], nbr[6], nbr[7], xi[0], xi[1], xi[2]);
#endif

                if(fh->S_function == moving_sphere)
                    //S = fhmd_Sxyz_r(x[n], fh->protein_com, fh);     // MD/FH sphere follows protein
                	S = fhmd_Sxyz_protein(x[n], fh->prot_com, fh);
                else if(fh->S_function == fixed_sphere)
                    S = fhmd_Sxyz_r(x[n], fh->box05, fh);           // Fixed MD/FH sphere
                else if(fh->S_function == AFM)
                	S = fhmd_Sxyz_AFM(x[n], fh->box05, fh);

                if (S < 0)
                {
                	printf (MAKE_RED "the S is negative S = %g\n",S);
                	exit (4);
                }

                if(fh->thermostat_function == bin)
                {
                	ibinz = find_ind(n,x,fh);
                	lg = ke_data[ibinz].lambda;

                	for (d = 0; d<DIM; d++)
                	{
                		if(fh->scheme == One_Way)
                		{
                			if(fh->S_function == AFM)
                			{
                				//vn = lg*v[n][d] + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + alpha_term[d] + S*(1 - S)*beta_term[d])*invro_dt;
                				vn           = lg*v[n][d] + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + alpha_term[d] + S*pow((1 - S),fh->gamma)*beta_term[d])*invro_dt;

                				v[n][d]      = vn;
                				//xprime[n][d] = x[n][d] + ((1 - S)*vn + S*u_fh[d])*dt + S*(1 - S)*grad_ro[d]*invro_dt;
                				xprime[n][d] = x[n][d] + ((1 - S)*vn + S*u_fh[d])*dt + S*pow((1 - S),fh->gamma)*grad_ro[d]*invro_dt;
                			}
                			else
                			{
                				printf(MAKE_RED "thermostat bin only supports the AFM option\n" RESET_COLOR "\n");
                				exit (3);
                			}
                		}

                		else
                		{
                			printf(MAKE_RED "bin thermostat is only supported by one-way coupling\n");
                			exit (2);
                		}
                	}
                }

                else if (fh->thermostat_function == gromacs)
                {
                	if (cTC)
                	{
                		gt = cTC[n];
                	}
                	lg = tcstat[gt].lambda;                             // Thermostat

                	/* Local thermostat */
                	if(fh->S_berendsen >= 0)
                	{
                		if(S > fh->S_berendsen)
                			lg = 1;
                	}
                	else
                	{
                		lg = lg*(1 - pow(S, -fh->S_berendsen)) + pow(S, -fh->S_berendsen);
                	}

                	for (d = 0; d < DIM; d++)
                	{
                		/* Pure MD: */
                		/* vn           = lg*v[n][d] + f[n][d]*w_dt; */
                		/* v[n][d]      = vn; */
                		/* xprime[n][d] = x[n][d] + vn*dt; */

                		if(fh->scheme == One_Way)
                		{
                			if(fh->S_function == AFM)
                			{
                					vn           = lg*v[n][d] + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + alpha_term[d] + S*pow((1 - S),fh->gamma)*beta_term[d])*invro_dt;
                					v[n][d]      = vn;
                			        xprime[n][d] = x[n][d] + ((1 - S)*vn + S*(u_fh[d] + u_fh_flow[d]))*dt + S*pow((1 - S),fh->gamma)*grad_ro[d]*invro_dt;
                			}

                			else if(fh->S_function == moving_sphere)
                			{
                				if (fh->flow_type == flow)
                				{


										vn = lg*v[n][d] + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + alpha_term[d] + S*(1 - S)*beta_term[d])*invro_dt;
										v[n][d]      = vn;
										xprime[n][d] = x[n][d] + ((1 - S)*vn + S*(u_fh[d] + u_fh_flow[d]))*dt + S*(1 - S)*grad_ro[d]*invro_dt;
										//xprime[n][d] = x[n][d] + ((1 - S)*vn + S*u_fh[d])*dt + S*(1 - S)*grad_ro[d]*invro_dt;



									   /* real velocity dx/dt */
									   if ((xprime[n][d] - x[n][d]) > fh->box05[d])
									   	   fh->vel[n][d] = (xprime[n][d] - x[n][d] - fh->box[d])/dt;
									   else if ((xprime[n][d] - x[n][d]) < -fh->box05[d])
									   	   fh->vel[n][d] = (xprime[n][d] - x[n][d] + fh->box[d])/dt;
									   else
									   	   fh->vel[n][d] = (xprime[n][d] - x[n][d])/dt;
								}

                				else
                				{
                					vn = lg*v[n][d] + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + alpha_term[d] + S*(1 - S)*beta_term[d])*invro_dt;
									v[n][d]      = vn;
									xprime[n][d] = x[n][d] + ((1 - S)*vn + S*u_fh[d])*dt + S*(1 - S)*grad_ro[d]*invro_dt;

                				}
                			}

                			else
                			{
                				vn = lg*v[n][d] + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + alpha_term[d] + S*(1 - S)*beta_term[d])*invro_dt;
                				v[n][d]      = vn;
                				xprime[n][d] = x[n][d] + ((1 - S)*vn + S*u_fh[d])*dt + S*(1 - S)*grad_ro[d]*invro_dt;
                			}
                		}
                		else if(fh->scheme == Two_Way)
                		{
                			gamma_u = fh->gamma_u*S*S*S*S*dt*(fh->stat.std_u_fh[d]*fh->stat.std_u_fh[d]/(fh->std_u*fh->std_u) - 1);
                			gamma_x = fh->gamma_x*S*S*S*S*dt*(fh->stat.std_rho_fh/fh->std_rho - 1);

                			if(fabs(gamma_u) < g_eps) gamma_u = g_eps;
                			if(fabs(gamma_x) < g_eps) gamma_x = g_eps;

                			vn = lg*v[n][d]*exp(-gamma_u) + ((1 - S)*f[n][d]*invmass[n] + (S*f_fh[d] + alpha_term[d] + S*(1 - S)*beta_term[d])*arr[ind].inv_ro)
                                 *(1 - exp(-gamma_u))/gamma_u*dt;

                			v[n][d] = vn;

                			xprime[n][d] = x[n][d] + (1 - S)*vn*(1 - exp(-gamma_u))/gamma_u*dt +
                                           (S*u_fh[d] + S*(1 - S)*grad_ro[d]*arr[ind].inv_ro)*(1 - exp(-gamma_x))/gamma_x*dt;
                		}
                	}//end of for (d)
                }//end of thermostat_function
            }//end of ptype[n]
            else
            {
                for (d = 0; d < DIM; d++)
                {
                    v[n][d]        = 0.0;
                    xprime[n][d]   = x[n][d];
                }
            }
        }
    }
}


void update_FH_MD_tcouple (gmx_int64_t       step,
						   t_inputrec       *inputrec,
						   t_state          *state,
						   gmx_ekindata_t   *ekind,
						   t_extmass        *MassQ,
						   t_mdatoms        *md,
						   FHMD             *fh
						   )
{
	gmx_bool   bTCouple = FALSE;
	real       dttc;
	int        i, offset;

	/* if using vv with trotter decomposition methods, we do this elsewhere in the code */
	if (inputrec->etc != etcNO &&
	        !(inputrecNvtTrotter(inputrec) || inputrecNptTrotter(inputrec) || inputrecNphTrotter(inputrec)))
	{
		/* We should only couple after a step where energies were determined (for leapfrog versions)
	    or the step energies are determined, for velocity verlet versions */

		if (EI_VV(inputrec->eI))
	    {
			offset = 0;
	    }
	    else
	    {
	    	offset = 1;
	    }
	        bTCouple = (inputrec->nsttcouple == 1 ||
	                    do_per_step(step+inputrec->nsttcouple-offset,
	                                inputrec->nsttcouple));
	}


	if (bTCouple)
	{
		dttc = inputrec->nsttcouple*inputrec->delta_t;
		switch (inputrec->etc)
		{
			case etcBERENDSEN:
				berendsen_FH_MD_tcoupl(inputrec, ekind, dttc,step,fh);
				break;

			default :
			{
				printf(MAKE_RED "FHMD error: FHMD thermostat only supports Berendsen thermostat \n");
				exit (1);
			}
		}

	}
}



void berendsen_FH_MD_tcoupl (t_inputrec *ir, gmx_ekindata_t *ekind, real dt,gmx_int64_t step,FHMD *fh)
{
	t_grpopts *opts;
    int        i,j,nzbin;
    real       T, reft = 0, lll;
	AFM_KE    *ke_data;


    opts = &ir->opts;
    nzbin = fh->nzbin;
	ke_data = fh->ke_data;
	//check the input temp and tau of each T-group, as it is also used as reference temp for each bin
    for (i = 0; (i < opts->ngtc); i++)
    {
    	for (j = i; (j< opts->ngtc); j++)
		{
    		if (opts->ref_t[i] != opts->ref_t[j] || opts->tau_t[i] != opts->tau_t[j])
		    {
    			printf("FHMD error: the input temperature or tau for T-group is not equal,please use same temperature and tau for T-group\n");
		    }
		}
    }

    for (i = 0; i<nzbin; i++)
    {
    	//get the S value at the bin center

	   ke_data[i].S =fhmd_Sxyz_d_AFM(fh->bin_c[i], fh->box05, fh);
	   //printf("the number of bin is %d and the S value is %g \n", i,ke_data[i].S);

	   if (ir->eI == eiVV)
	   {
		   T = ke_data[i].T;
	   }

	   else
	   {
		   T = ke_data[i].Th;
	   }

	   if(T > 0.0)
	   {
		   if(T < 0.1 && step == 0)
		   {
			   printf(MAKE_RED "FHMD error: please set the gen-vel in mdp file to yes\n");
			   exit (2);
		   }

		   //reft             = (1 - ke_data[i].S)*opts->ref_t[0] + ke_data[i].S*T;
		    reft             = (1 - pow(ke_data[i].S,fh->eta))*opts->ref_t[0] + pow(ke_data[i].S,fh->eta)*T;
		    ke_data[i].reft  = reft;
		    lll              = std::sqrt(1.0 + (dt/opts->tau_t[0])*(reft/T-1.0));
		    //ke_data[i].lambda       = std::max<real>(std::min<real>(lll, 1.8), 0.3);
		    ke_data[i].lambda = lll;

	   }

	   else
	   {
		   ke_data[i].lambda = 1.0;
	   }

	   //printf("the bin number is %d and the lambda is %g\n",i,ke_data[i].lambda);
    }
}

static void calc_FHMD_ke_part_normal (rvec v[], rvec x[], t_grpopts *opts, t_mdatoms *md,
        						      gmx_ekindata_t *ekind, t_nrnb *nrnb, gmx_bool bEkinAveVel,gmx_mtop_t *mtop,FHMD *fh)
{

	/*original code of Gromacs for KE calculation */
	 int           nthread, thread;
     int           g;
     t_grp_tcstat *tcstat  = ekind->tcstat;
     t_grp_acc    *grpstat = ekind->grpstat;

	     /*three main: VV with AveVel, vv with AveEkin, leap with AveEkin.  Leap with AveVel is also
	       an option, but not supported now.
	       bEkinAveVel: If TRUE, we sum into ekin, if FALSE, into ekinh.*/


	     /*group velocities are calculated in update_ekindata and
	     * accumulated in acumulate_groups.
	     * Now the partial global and groups ekin.*/

	    for (g = 0; (g < opts->ngtc); g++)
	    {
	        copy_mat(tcstat[g].ekinh, tcstat[g].ekinh_old);
	        if (bEkinAveVel)
	        {
	            clear_mat(tcstat[g].ekinf);
	            tcstat[g].ekinscalef_nhc = 1.0;    //need to clear this -- logic is complicated!
	        }
	        else
	        {
	            clear_mat(tcstat[g].ekinh);
	        }
	    }
	    ekind->dekindl_old = ekind->dekindl;
	    nthread            = gmx_omp_nthreads_get(emntUpdate);

	#pragma omp parallel for num_threads(nthread) schedule(static)
	    for (thread = 0; thread < nthread; thread++)
	    {
	        // This OpenMP only loops over arrays and does not call any functions
	        // or memory allocation. It should not be able to throw, so for now
	        // we do not need a try/catch wrapper.
	        int     start_t, end_t, n;
	        int     ga, gt;
	        rvec    v_corrt;
	        real    hm;
	        int     d, m;
	        matrix *ekin_sum;
	        real   *dekindl_sum;

	        start_t = ((thread+0)*md->homenr)/nthread;
	        end_t   = ((thread+1)*md->homenr)/nthread;

	        ekin_sum    = ekind->ekin_work[thread];
	        dekindl_sum = ekind->dekindl_work[thread];

	        for (gt = 0; gt < opts->ngtc; gt++)
	        {
	            clear_mat(ekin_sum[gt]);
	        }
	        *dekindl_sum = 0.0;

	        ga = 0;
	        gt = 0;
	        for (n = start_t; n < end_t; n++)
	        {
	            if (md->cACC)
	            {
	                ga = md->cACC[n];
	            }
	            if (md->cTC)
	            {
	                gt = md->cTC[n];
	            }
	            hm   = 0.5*md->massT[n];

	            for (d = 0; (d < DIM); d++)
	            {
	                v_corrt[d]  = v[n][d]  - grpstat[ga].u[d];
	            }
	            for (d = 0; (d < DIM); d++)
	            {
	                for (m = 0; (m < DIM); m++)
	                {
	                     /*if we're computing a full step velocity, v_corrt[d] has v(t).  Otherwise, v(t+dt/2)*/
	                    ekin_sum[gt][m][d] += hm*v_corrt[m]*v_corrt[d];
	                }
	            }
	            if (md->nMassPerturbed && md->bPerturbed[n])
	            {
	                *dekindl_sum +=
	                    0.5*(md->massB[n] - md->massA[n])*iprod(v_corrt, v_corrt);
	            }
	        }
	    }

	    ekind->dekindl = 0;
	    for (thread = 0; thread < nthread; thread++)
	    {
	        for (g = 0; g < opts->ngtc; g++)
	        {
	            if (bEkinAveVel)
	            {
	                m_add(tcstat[g].ekinf, ekind->ekin_work[thread][g],
	                      tcstat[g].ekinf);
	            }
	            else
	            {
	                m_add(tcstat[g].ekinh, ekind->ekin_work[thread][g],
	                      tcstat[g].ekinh);
	            }
	        }

	        ekind->dekindl += *ekind->dekindl_work[thread];
	    }

	    inc_nrnb(nrnb, eNR_EKIN, md->homenr);

        /*KE calculation for each bin*/

	    int    i, j, k, n, d, m, ibinz, ga, gt, res_nr,nwater_atoms,ngraphene_atoms;
	    int    start = 0, homenr = md->homenr;
	    real   dekindl,hm;
	    double   count[fh->nzbin],count_water[fh->nzbin],count_graphene[fh->nzbin];
	    rvec   v_corrt;
	    int    nzbin, ngtc;
	    char  *atomname, *resname;

	    AFM_KE       *ke_data = fh->ke_data;
	    tensor         **ke_data_work_ekinh = fh->ke_data_work_ekinh;
	    tensor         **ke_data_work_ekinh_old = fh->ke_data_work_ekinh_old;
	    //t_grp_acc    *grpstat = ekind->grpstat;
	    //t_grp_tcstat *tcstat  = ekind->tcstat;

	    atomname = NULL;
	    resname  = NULL;

	    ngtc = opts->ngtc;
	    nzbin = fh->nzbin;
	    nwater_atoms = 0;
	    ngraphene_atoms = 0;


	    for(i = 0; i< nzbin; i++)
	    {
	    	copy_mat(ke_data[i].ekinh,ke_data[i].ekinh_old);
	        clear_mat(ke_data[i].ekinh);
	        count[i] = 0;
	        count_water[i] = 0;
	        count_graphene[i] = 0;
	    }

	    for(i = 0; i< nzbin; i++)
	    {
	    	for (j =0; j< ngtc; j++)
	        {
	    		copy_mat(ke_data_work_ekinh[i][j], ke_data_work_ekinh_old[i][j]);
	        	clear_mat(ke_data_work_ekinh[i][j]);
	        }
	    }

	    ekind->dekindl_old = ekind->dekindl;
	    dekindl = 0;

	    ga = 0;
	    gt = 0;

	    for (n = start; n < start+homenr; n++)
	    {
	    	ibinz = find_ind(n,x,fh); /*assign each atom to each bin*/
	    	count[ibinz]++;

	    	/*count the number of water molecules in each bin*/
	    	gmx_mtop_atominfo_global(mtop, n, &atomname, &res_nr, &resname);
	    	if(!(strcmp(resname,"SOL")))
	    	{
	    		count_water[ibinz]++;
	    		nwater_atoms++;
	    	}

	    	gmx_mtop_atominfo_global(mtop, n, &atomname, &res_nr, &resname);
	    	if(!(strcmp(resname,"GRA")))
	    	{
	    		count_graphene[ibinz]++;
	    		ngraphene_atoms++;
	    	}

	    	if (md->cACC)
	    	{
	    		ga = md->cACC[n];
	    	}

	    	if (md->cTC)
	    	{
	    		gt = md->cTC[n];
	    	}

	    	hm   = 0.5*md->massT[n];

	    	for (d = 0; (d < DIM); d++)
	    	{
	    		//printf("The velocity is =%g \n",v[n][d]);
	    		v_corrt[d]  = v[n][d]  - grpstat[ga].u[d];
	    	}

	    	for (d = 0; (d < DIM); d++)
	    	{
	    		for (m = 0; (m < DIM); m++)
	    		{

	    			ke_data_work_ekinh[ibinz][gt][m][d] += hm*v_corrt[m]*v_corrt[d];
	    		}
	    	}

	    	if (md->nPerturbed && md->bPerturbed[n])
	    	{
	    		/* The minus sign here might be confusing.
	    		* The kinetic contribution from dH/dl doesn't come from
	    		* d m(l)/2 v^2 / dl, but rather from d p^2/2m(l) / dl,
	    		* where p are the momenta. The difference is only a minus sign.
	    		*/
	    		 dekindl -= 0.5*(md->massB[n] - md->massA[n])*iprod(v_corrt, v_corrt);
	    	}
	    }

	    /*calculate the degree of freedom of each bin*/

	    for(i = 0; i<nzbin; i++)
	    {
	    	//ke_data[i].nrdf = 3*count[i]*(3*atoms->nr - 6)/3*atoms->nr;
	    	//ke_data[i].nrdf = (3*count[i] - count_water[i])*(static_cast<double> (3*md->homenr - nwater_atoms - 3)/(3*md->homenr - nwater_atoms));
	    	ke_data[i].nrdf = (3*count[i] - count_water[i]- 3*count_graphene[i])*\
	    				(static_cast<double> (3*md->homenr - nwater_atoms - 3*ngraphene_atoms - 3)/(3*md->homenr - nwater_atoms - 3*ngraphene_atoms));
	    	//printf("the bin num is = %d  the nrdf = %g \n",i,ke_data[i].nrdf);
	    }


	    //pass the KE in ke_data_work to ke_data and tcstat


	    for(i = 0; i<nzbin; i++)
	    {
	    	for (j = 0; j<ngtc; j++)
	    	{
	    		//m_add(tcstat[j].ekinh,ke_data_work_ekinh[i][j],
	    		//tcstat[j].ekinh);
	    		m_add(ke_data[i].ekinh,ke_data_work_ekinh[i][j],
	    			  ke_data[i].ekinh);
	    	}

	    }

	    	//ekind->dekindl = dekindl;
	    	//inc_nrnb(nrnb, eNR_EKIN, homenr);
}




void calc_FHMD_ke_part(t_state *state, t_grpopts *opts, t_mdatoms *md,
						gmx_ekindata_t *ekind, t_nrnb *nrnb, gmx_bool bEkinAveVel,gmx_mtop_t *mtop,FHMD *fh)
{
	if (ekind->cosacc.cos_accel == 0)
	{
		calc_FHMD_ke_part_normal(state->v, state->x,opts, md, ekind, nrnb, bEkinAveVel,mtop,fh);

		/*pass the KE in ke_data_work to ke_data and tcstat*/
		//pass_ke_data(ekind,opts,fh);
	}

	else
	{
		printf(MAKE_RED "FHMD error: NEMD is not supported\n");
	}
}



