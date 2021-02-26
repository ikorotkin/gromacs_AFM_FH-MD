#include "data_structures.h"
#include "parser.h"
#include "fh_functions.h"
#include "estimate.h"
#include "sfunction.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdlib/rbin.h"
#include "flow_velocity.h"
//#include "gromacs/math/vec.h"
#include "macro.h"
#include "gromacs/utility/smalloc.h"

int fhmd_init(matrix box, int N_atoms, real mass[], rvec x[], double dt_md, t_grpopts *opts ,gmx_mtop_t *mtop, t_commrec *cr,gmx_fhmd_global_stat *g_fhmd_s ,FHMD *fh)
{
    int N_atoms_th = N_atoms;
    int ngtc;

    ngtc = opts->ngtc;

    FILE *fw;

    if(MASTER(cr))
    {
        char const *fname_in  = "coupling.prm";
        char const *fname_out = "coupling_out.prm";


        /* Initial output */

        printf(MAKE_GREEN "\n  Aston University, Department of Mathematics, Dmitry Nerukh Research Group\n");
        printf("  Queen Mary University of London, School of Engineering and Material Science\n\n");
        printf(MAKE_PURPLE "     Hybrid Molecular Dynamics - 2-Way Coupling Parallel VERSION %4.2f\n\n", FHMD_VERSION);

        /* Default values of FHMD parameters */

        fh->scheme      = 2;
        fh->S           = 0;
        fh->R1          = 0.5;
        fh->R2          = 1;
        fh->z1          = 0.1;
        fh->z2          = 0.5;
        fh->nzbin       = 3;
        fh->Smin        = 0;
        fh->Smax        = 0.5;
        fh->alpha       = 50;
        fh->beta        = 20;
        fh->gamma_x     = 0;
        fh->gamma_u     = 1;
        fh->eps_rho     = 0.05;
        fh->eps_mom     = 0.05;
        fh->S_berendsen = 1;
        fh->N[0]        = 9;
        fh->N[1]        = 9;
        fh->N[2]        = 9;
        fh->N_md[0]     = 5; //fh->N[0];
        fh->N_md[1]     = 5; //fh->N[1];
        fh->N_md[2]     = 5; //fh->N[2];
        fh->FH_EOS      = 1;
        fh->FH_step     = 10;
        fh->FH_equil    = 10000;
        fh->FH_dens     = 600;
        fh->FH_temp     = 298.15;
        fh->FH_blend    = 0.005;
        fh->Noutput     = 100;

        /* Read FHMD parameters */

        printf(MAKE_GREEN "FHMD: Reading parameters from %s..." RESET_COLOR " ", fname_in);

        int ok = parse_prm(fname_in, fh);

        if(ok == 1) {
            printf(MAKE_GREEN "...OK\n" RESET_COLOR "\n");
            fw = fopen(fname_out, "w");
        } else if(ok == -1) {
            printf(MAKE_RED "\nFHMD: File %s not found. Generating default parameter file...\n" RESET_COLOR "\n", fname_in);
            fw = fopen(fname_in, "w");
        } else {
            printf(MAKE_RED "\nFHMD: ERROR in %s file\n" RESET_COLOR "\n", fname_in);
            exit(2);
        }

        /* Print FHMD parameters to the screen and output file */

        fprintf(fw, "; Hybrid Molecular Dynamics - 2-Way Coupling Parallel VERSION %4.2f\n\n", FHMD_VERSION);

        fprintf(fw, "scheme = %d              ; 0 - Pure MD, 1 - One-way coupling, 2 - Two-way coupling\n\n", fh->scheme);

        fprintf(fw, "S = %g                   ; Parameter S -1 - fixed sphere, -2 - moving sphere, -3 - AFM\n\n", fh->S);
        fprintf(fw, "R1   = %g              ; MD sphere radius for variable S, [0..1]\n", fh->R1);
        fprintf(fw, "R2   = %g              ; FH sphere radius for variable S, [0..1]\n", fh->R2);
        fprintf(fw, "Smin = %g                ; Minimum S for variable S\n", fh->Smin);
        fprintf(fw, "Smax = %g             ; Maximum S for variable S\n\n", fh->Smax);

        switch(fh->scheme)
        {
        case Pure_MD:
            printf(MAKE_RED "FHMD: Starting Pure MD simulation (without MD/FH coupling)\n" RESET_COLOR "\n");
            break;
        case One_Way:
            printf(MAKE_YELLOW "FHMD: One-way MD/FH coupling\n" RESET_COLOR "\n");
            break;
        case Two_Way:
            break;
        }

        if(fh->S >= 0.0)
        {
            printf(MAKE_PURPLE "FHMD: S = %g\n", fh->S);
            fh->S_function = constant_S;
        }
        else
        {
            printf(MAKE_PURPLE "FHMD: S = S(x,y,z) = [%g, %g] with R1 = %g, R2 = %g\n", fh->Smin, fh->Smax, fh->R1, fh->R2);
            if(fh->S > -1.5)
            {
            	fh->S_function = fixed_sphere;
            }

            //fh->z1 *= box[2][2];
            //fh->z2 *= box[2][2];
            fh->R1 *= box[0][0]*0.5;
            fh->R2 *= box[0][0]*0.5;
            fh->R12 = fh->R1*fh->R1;
            fh->R22 = fh->R2*fh->R2;
            fh->RS  = (fh->Smax - fh->Smin)/(fh->R2 - fh->R1);
            printf(MAKE_GREEN "FHMD: Absolute values of R [nm]: R1 = %f, R2 = %f\n", fh->R1, fh->R2);

            if(fh->S < -1.5 && fh->S > -2.5)
            {
                printf(MAKE_PURPLE "FHMD: The MD/FH sphere will follow the protein\n");
                fh->S_function = moving_sphere;
             }

            else if(fh->S < -2.5)
            {
            	fh->S_function = AFM;

            	fh->z1 *= box[2][2]*0.5;
            	fh->z2 *= box[2][2]*0.5;
            	fh->zs  = (fh->Smax - fh->Smin)/(fh->z2 - fh->z1);
            	fh->z_c_ball *= box[2][2]*0.5;

            	if(fh->N[2]%2 !=0)
            	{
            		printf(MAKE_RED "The grids in Z direction is not even,use a even number" RESET_COLOR "\n");
            		exit(1);

            	}


            	 printf(MAKE_GREEN "FHMD: Absolute values of z [nm]: z1 = %f, z2 = %f ,R1 = %f, R2 = %f\n", fh->z1, fh->z2,fh->R1,fh->R2);
            	 printf(MAKE_GREEN "FHMD: Absolute values of z [nm]: z_c_ball = %f \n", fh->z_c_ball);
            	 printf(MAKE_GREEN "gamma = %g ; the factor to adjust the alpha and beta term\n", fh->gamma);
            	 fprintf(fw, "R1  = %g              ; the  MD sphere radius AFM tip\n",fh->R1);
            	 fprintf(fw, "R2  = %g              ; the  MD sphere radius AFM tip\n",fh->R2);
            	 fprintf(fw,"z1   = %g              ; upper limit of coordinate of upper MD zone in z direction  for variable\n", fh->z1);
            	 fprintf(fw,"z2   = %g              ; upper limit of coordinate of upper buffer zone in z direction for variable S\n", fh->z2);
            	 fprintf(fw,"z_c_ball = %g              ; center of the AFM tip S\n", fh->z_c_ball);
            	 fprintf(fw,"gamma = %g             ; the factor to adjust the alpha and beta term\n", fh->gamma);
            }
        }

        if (fh->flow == 0)
        {
        	fh->flow_type = no_flow;
            printf(MAKE_PURPLE "no flow field is added\n");
        }

		else if (fh->flow == 1)
		{
			fh->flow_type = flow;
			printf(MAKE_PURPLE "flow field will be added\n");
			fprintf(fw,"nzbin   = %d                ; number of bins in z direction\n", fh->nzbin);
		 }

		else
		{
			printf(MAKE_RED "please chose type of flow" RESET_COLOR "\n");
			exit(3);
		}

		fprintf(fw,"flow_type = %d			; for  moving sphere and AFM in one-way(0 - no flow , 1 - with flow) \n",fh->flow_type);



        if (fh->thermostat == 0)
        {
        	fh->thermostat_function = gromacs;
        	printf(MAKE_GREEN "thermostat = %g		;the original gromacs thermostat is used\n",fh->thermostat);
        	fprintf(fw,"thermostat = %g			;the original gromacs thermostat is used\n", fh->thermostat);
        }

        else if (fh->thermostat == 1)
        {
        	fh->thermostat_function = bin;
        	printf(MAKE_GREEN "FHMD: thermostat = %g		;Warning:the thermostat according to bin is used,but results are not tested\n",fh->thermostat);
        	printf(MAKE_GREEN "FHMD: eta = %g    		    ;the factor to adjust the thermostat\n", fh->eta);
        	printf(MAKE_GREEN "FHMD: bins in z direction for thermostat nzbin = %d\n", fh->nzbin);
        	fprintf(fw,"thermostat = %g			;the thermostat according to bin is used\n", fh->thermostat);
        	fprintf(fw,"eta = %g               		;the factor to adjust the thermostat\n", fh->eta);
        	fprintf(fw,"nzbin   = %d                ; number of bins in z direction\n", fh->nzbin);
        }

        else
        {
        	printf(MAKE_RED "please chose type of thermostat" RESET_COLOR "\n");
        	exit(4);
        }

        printf(MAKE_GREEN "FHMD: alpha = %g [nm^2/ps], beta = %g [ps^-1]\n", fh->alpha, fh->beta);
        fprintf(fw, "alpha   = %g           ; Alpha parameter for dx/dt and du/dt equations, nm^2/ps\n", fh->alpha);
        fprintf(fw, "beta    = %g           ; Beta parameter for du/dt equation, ps^-1\n\n", fh->beta);

        if(fh->scheme == Two_Way)
        {
            printf("FHMD: MD dissipator: gamma_x = %g [ps^-1], gamma_u = %g [ps^-1]\n", fh->gamma_x, fh->gamma_u);
            printf("FHMD: FH dissipator: eps_rho = %g [--], eps_mom = %g [--]\n", fh->eps_rho, fh->eps_mom);
            fprintf(fw, "gamma_x = %g             ; Gamma_x parameter (MD density fluctuations dissipator), ps^-1\n", fh->gamma_x);
            fprintf(fw, "gamma_u = %g             ; Gamma_u parameter (MD velocity fluctuations dissipator), ps^-1\n", fh->gamma_u);
            fprintf(fw, "eps_rho = %g             ; Eps_rho parameter (FH density fluctuations dissipator)\n", fh->eps_rho);
            fprintf(fw, "eps_mom = %g             ; Eps_mom parameter (FH momentum fluctuations dissipator)\n\n", fh->eps_mom);
        }

        if(fh->S_berendsen >= 0)
            printf("FHMD: Berendsen thermostat works for S <= %g\n", fh->S_berendsen);
        else
            printf("FHMD: Berendsen thermostat with (1 - S^%g) multiplier\n", -fh->S_berendsen);
        fprintf(fw, "S_berendsen = %g         ; If S_berendsen >= 0, Berendsen thermostat will work for S <= S_berendsen,\n", fh->S_berendsen);
        fprintf(fw, "                        ; otherwise factor (1-S^(-S_berendsen)) will be applied (local thermostat)\n\n");

        for(int d = 0; d < DIM; d++)
        {
            fh->box[d]         = box[d][d];
            fh->box05[d]       = 0.5*fh->box[d];
            //fh->protein_com[d] = fh->box05[d];
            fh->N_shift[d]     = (fh->N[d] - fh->N_md[d])/2;
        }

        fh->box_volume = fh->box[0]*fh->box[1]*fh->box[2];

        printf("FHMD: MD box size:  %g x %g x %g [nm]\n", fh->box[0], fh->box[1], fh->box[2]);
        printf("FHMD: FH grid size: %d x %d x %d\n", fh->N[0], fh->N[1], fh->N[2]);
        fprintf(fw, "Nx = %d                  ; Number of FH cells along X axis\n", fh->N[0]);
        fprintf(fw, "Ny = %d                  ; Number of FH cells along Y axis\n", fh->N[1]);
        fprintf(fw, "Nz = %d                  ; Number of FH cells along Z axis\n\n", fh->N[2]);

        if(fh->scheme == Two_Way)
        {
            printf("FHMD: Small-scale MD/FH grid size: %d x %d x %d\n", fh->N_md[0], fh->N_md[1], fh->N_md[2]);
            fprintf(fw, "NxMD = %d                ; Number of small-scale MD-FH cells along X axis\n", fh->N_md[0]);
            fprintf(fw, "NyMD = %d                ; Number of small-scale MD-FH cells along Y axis\n", fh->N_md[1]);
            fprintf(fw, "NzMD = %d                ; Number of small-scale MD-FH cells along Z axis\n\n", fh->N_md[2]);
        }
        else
        {
            for(int d = 0; d < DIM; d++)
            {
                fh->N_shift[d] = 0;
                fh->N_md[d]    = fh->N[d];
            }
        }

        fh->Ntot    = fh->N[0]*fh->N[1]*fh->N[2];
        fh->Ntot_md = fh->N_md[0]*fh->N_md[1]*fh->N_md[2];

        switch(fh->FH_EOS)
        {
        case 0:
            fh->eos = eos_argon;
            printf("FHMD: Equation of state: Liquid Argon (300K)\n");
            break;
        case 1:
            fh->eos = eos_spce;
            printf("FHMD: Equation of state: Rigid SPC/E water\n");
            break;
        default:
            printf(MAKE_RED "\nFHMD: Unknown equation of state (%d) in %s\n" RESET_COLOR "\n", fh->FH_EOS, fname_in);
            exit(18);
        }

        fprintf(fw, "FH_EOS   = %d            ; EOS: 0 - Liquid Argon, 1 - SPC/E water\n", fh->FH_EOS);

        fh->dt_FH = (double)(fh->FH_step)*dt_md;

        printf("FHMD: FH time step dt_FH = %d * dt_MD = %g [ps]\n", fh->FH_step, fh->dt_FH);
        fprintf(fw, "FH_step  = %d           ; FH time step dt_FH = FH_step * dt_MD\n", fh->FH_step);

        if(fh->scheme == Two_Way) fh->FH_equil = 0;
        printf("FHMD: FH equilibration steps: %d\n", fh->FH_equil);
        fprintf(fw, "FH_equil = %d        ; Number of time steps for the FH model equilibration (for 1-way coupling)\n", fh->FH_equil);

        printf("FHMD: FH Density = %g [amu/nm^3], FH Temperature = %g [K]\n", fh->FH_dens, fh->FH_temp);
        fprintf(fw, "FH_dens  = %g          ; FH mean density\n", fh->FH_dens);
        fprintf(fw, "FH_temp  = %g       ; FH mean temperature\n", fh->FH_temp);
        fprintf(fw, "FH_blend = %g        ; FH Blending: -1 - dynamic, or define static blending parameter (0..1)\n\n", fh->FH_blend);

        printf("FHMD: MD/FH arrays will be written every %d MD time steps\n", fh->Noutput);
        fprintf(fw, "Noutput  = %d           ; Write arrays to files every Noutput MD time steps (0 - do not write)\n", fh->Noutput);

        printf(RESET_COLOR "\n");

        fflush(stdout);
    } // if(MASTER(cr))

    fhmd_reset_statistics(fh);

    fh->total_density = 0;
    for(int i = 0; i < N_atoms_th; i++)
        fh->total_density += mass[i];

    /* Broadcast parameters to all threads */



    if(PAR(cr))
    {
    	//printf("The MPI code stopted after PAR(cr) before gmx_sumd\n");

        gmx_sumd(1, &fh->total_density, cr);
        gmx_sumi(1, &N_atoms, cr);
        gmx_bcast(sizeof(FHMD), fh, cr);
    }

    fh->total_density /= fh->box_volume;


	if (fh->S_function == moving_sphere)
	{

		fh->prot_com = (dvec*)calloc(fh->prot_num,sizeof(dvec));
		fh->prot_mass = (double*)calloc(fh->prot_num,sizeof(double));
		fh->prot_N = (int*)calloc(fh->prot_num,sizeof(int));
		fh->prot_name = (char **)calloc(fh->prot_num,sizeof(fh->tem));

		for (int i = 0; i < fh->prot_num; i++)
		{
			fh->prot_name[i] = (char*)calloc(1,sizeof(fh->tem));
		}
		if(fh->prot_com == NULL || fh->prot_mass == NULL || fh->prot_N == NULL || fh->prot_name == NULL)
		{
			printf(MAKE_RED "FHMD error:Allocation of memory for protein failed" RESET_COLOR "\n");
			exit(4);
		}

		char delim[] = " ";
		char *ptr = strtok(fh->tem, delim);
		char    name[fh->prot_num][800];

		for (int i =0; i < fh->prot_num; i++)
		{
			strcpy(fh->prot_name[i],ptr);
			ptr = strtok(NULL, delim);
			fprintf(fw,"protein_name = %s \n", fh->prot_name[i]);

		}


		fhmd_find_protein(mtop, N_atoms_th, mass, cr, fh);
		fhmd_find_protein_com(mtop, N_atoms_th, x, mass, cr, fh);
	}

	//fhmd_find_protein_com(top_global, mdatoms->homenr, state->x, mdatoms->massT, cr, &fhmd);

    if (fh->eos != eos_spce)
    {
    	printf(MAKE_RED "Water must be used for AFM and moving sphere" RESET_COLOR "\n");
    	exit (5);
    }

    if(MASTER(cr))
    {
        printf(MAKE_GREEN "FHMD: Total number of atoms in the box: %d\n", N_atoms);
        printf("FHMD: Total density of the box: %g [amu/nm^3]\n", fh->total_density);

        //if(fh->protein_N > 0)
            //printf(MAKE_PURPLE "FHMD: Found protein: %d atoms, mass = %g [amu], COM = (%g, %g, %g) [nm]\n",
                    //fh->protein_N, fh->protein_mass, fh->protein_com[0], fh->protein_com[1], fh->protein_com[2]);

        if (fh->S_function == moving_sphere)
        {
			for(int i = 0; i < fh->prot_num; i++)
			{
				if(fh->prot_N[i] >= 0)
				{
					printf(MAKE_PURPLE "FHMD: Found protein: %d atoms, mass = %g [amu], COM = (%g, %g, %g) [nm]\n",
							fh->prot_N[i], fh->prot_mass[i],fh->prot_com[i][0], fh->prot_com[i][1], fh->prot_com[i][2]);
				}

			}
        }
        printf(RESET_COLOR "\n");
        fflush(stdout);

        fprintf(fw, "\n; You may consider to use FH_dens = %g since this is the total MD density of the box.\n", fh->total_density);
        fprintf(fw, "; NB: Please use spaces before and after '=' in this file, e.g. 'S = 0' (not 'S=0').\n");
        fclose(fw);
    }


    if(fh->scheme == Pure_MD) return 0;     // Start Pure MD simulation

    /* Allocate memory */

    fh->arr  = (FH_arrays*)calloc(fh->Ntot, sizeof(FH_arrays));
    fh->ind  = (int*)calloc(N_atoms, sizeof(int));
    fh->indv = (ivec*)calloc(N_atoms, sizeof(ivec));



    fh->mpi_linear = (double*)malloc(8*fh->Ntot*sizeof(double));   // 8 components: ro_md, uro_md[3], ro_md_s, uro_md_s[3]

    if(fh->arr == NULL || fh->ind == NULL || fh->indv == NULL || fh->mpi_linear == NULL)
    {
        if(MASTER(cr)) printf(MAKE_RED "\nFHMD: ERROR: Out of memory (array allocator)\n" RESET_COLOR "\n");
        fflush(stdout);
        exit(3);
    }

    if (fh->flow_type == flow)
    {
    	fh->vx_m = (double*)calloc(fh->nzbin,sizeof(double));
		fh->sum_vx_m = (double*)calloc(fh->nzbin,sizeof(double));
		fh->vx = (double*)calloc(fh->nzbin,sizeof(double));
		fh->mass_bin = (double*)calloc(fh->nzbin,sizeof(double));

		fh->vel_mass = (double*)calloc(fh->nzbin,sizeof(double));
		fh->sum_vel = (double*)calloc(fh->nzbin,sizeof(double));
		fh->vel_m = (double*)calloc(fh->nzbin,sizeof(double));

		fh->vel  = (rvec*)calloc(N_atoms, sizeof(dvec));

		for(int i = 0; i< fh->nzbin; i++)
		{
			fh->vx_m[i] = 0.0;
			fh->sum_vx_m[i] = 0.0;
			fh->sum_vel[i] = 0.0;
		}

		if(fh->vx_m == NULL || fh->sum_vx_m == NULL || fh->vx == NULL || fh->mass_bin == NULL || fh->vel == NULL)
		{
			printf(MAKE_RED "FHMD error:Allocation of memory for flow filed failed" RESET_COLOR "\n");
			exit(5);
		}

		if(fh->vel_mass == NULL || fh->sum_vel == NULL || fh->vel_m == NULL)
		{
			printf(MAKE_RED "FHMD error:Allocation of memory for flow field failed" RESET_COLOR "\n");
			exit(6);
		}
    }



    fh->grid.c    = (dvec*)malloc(fh->Ntot*sizeof(dvec));
    fh->grid.n    = (dvec*)malloc(fh->Ntot*sizeof(dvec));
    fh->grid.h    = (dvec*)malloc(fh->Ntot*sizeof(dvec));
    fh->grid.vol  = (double*)malloc(fh->Ntot*sizeof(double));
    fh->grid.ivol = (double*)malloc(fh->Ntot*sizeof(double));
    fh->grid.md   = (FHMD_CELL*)malloc(fh->Ntot*sizeof(FHMD_CELL));

    if(fh->grid.c == NULL || fh->grid.n == NULL || fh->grid.h == NULL || fh->grid.vol == NULL || fh->grid.ivol == NULL)
    {
        if(MASTER(cr)) printf(MAKE_RED "\nFHMD: ERROR: Out of memory (FH grid allocator)\n" RESET_COLOR "\n");
        fflush(stdout);
        exit(3);
    }

    fh->stat.avg_rho_md_cell = (double*)calloc(fh->Ntot, sizeof(double));
    fh->stat.avg_rho_fh_cell = (double*)calloc(fh->Ntot, sizeof(double));

    if(fh->stat.avg_rho_md_cell == NULL || fh->stat.avg_rho_fh_cell == NULL)
    {
        if(MASTER(cr)) printf(MAKE_RED "\nFHMD: ERROR: Out of memory (Statistics allocator)\n" RESET_COLOR "\n");
        fflush(stdout);
        exit(3);
    }

    /* Allocate memory for AFM simulation*/
    //if(fh->S_function == AFM)
    if(fh->thermostat_function == bin)
    {
       	fh->ke_data = (AFM_KE*)calloc(fh->nzbin,sizeof(AFM_KE));
       	fh->nrdf_mpi = (double*)calloc(fh->nzbin,sizeof(double));
       	fh->bin_c = (dvec*)calloc(fh->nzbin,sizeof(dvec));


       	for(int i = 0; i< fh->nzbin; i++)
        {
           	fh->ke_data[i].nrdf = 0.0;
           	fh->ke_data[i].Th = 0.0;
           	fh->ke_data[i].T = 0.0;
           	fh->ke_data[i].lambda = 0.0;
        	fh->ke_data[i].S = 0.0;
        	fh->ke_data[i].reft = 0.0;
           	clear_mat(fh->ke_data[i].ekinh);
           	clear_mat(fh->ke_data[i].ekinh_old);
           	clear_mat(fh->ke_data[i].ekinf);
           	clear_dvec(fh->bin_c[i]);
        }

       	fh->ke_data_work_ekinh = (tensor**)calloc(fh->nzbin,ngtc*sizeof(tensor));
       	fh->ke_data_work_ekinh_old = (tensor**)calloc(fh->nzbin,ngtc*sizeof(tensor));

       	for(int i = 0; i<fh->nzbin; i++)
       	{
       		fh->ke_data_work_ekinh[i] = (tensor *)calloc(ngtc,sizeof(tensor));
       		fh->ke_data_work_ekinh_old[i] = (tensor*)calloc(ngtc,sizeof(tensor));
       	}

       	for(int i = 0; i<fh->nzbin; i++)
       	{
       		for (int j = 0; j < ngtc; j++)
       		{
       			clear_mat(fh->ke_data_work_ekinh[i][j]);
       			clear_mat(fh->ke_data_work_ekinh_old[i][j]);
       		}
       	}

       	//g_fhmd_s = (gmx_fhmd_global_stat*)calloc(1,sizeof(gmx_fhmd_global_stat));
       	g_fhmd_s->rb_bin = (t_bin*)calloc(1,sizeof(t_bin));
       	g_fhmd_s->itc0_bin = (int*)calloc(fh->nzbin,sizeof(int));
       	g_fhmd_s->itc1_bin = (int*)calloc(fh->nzbin,sizeof(int));

       	if(fh->ke_data == NULL || fh->ke_data_work_ekinh == NULL || fh->ke_data_work_ekinh_old == NULL || fh->nrdf_mpi == NULL || fh->bin_c == NULL)
       	{
       		if(MASTER(cr)) printf(MAKE_RED "\nFHMD: ERROR: Out of memory (AFM kinetic energy data allocator)\n" RESET_COLOR "\n");
       		fflush(stdout);
       		exit(3);
       	}

       	if(g_fhmd_s->rb_bin == NULL)
       	{
       		if(MASTER(cr)) printf(MAKE_RED "\nFHMD: ERROR: Out of memory (AFM global state data rb_bin allocator)\n" RESET_COLOR "\n");
       		fflush(stdout);
       		exit(3);
       	}

       	if(g_fhmd_s->itc0_bin == NULL)
       	{
       		if(MASTER(cr)) printf(MAKE_RED "\nFHMD: ERROR: Out of memory (AFM global state data itc0_bin allocator)\n" RESET_COLOR "\n");
       		fflush(stdout);
       		exit(3);
       	}

       	if(g_fhmd_s->itc1_bin == NULL)
       	{
       		if(MASTER(cr)) printf(MAKE_RED "\nFHMD: ERROR: Out of memory (AFM global state data itc1_bin allocator)\n" RESET_COLOR "\n");
       		fflush(stdout);
       		exit(3);
       	}

   }

       /*end of allocating memory for AFM*/

    /* Create FH grid and initialize FH solver */

    define_FH_grid(cr, fh);

    if (fh->thermostat_function == bin)
    {
    	find_bin_center(fh);
    }


	 if (fh->S_function == AFM)
	 {
		char const *fname_in_velocity = "vel2fhmd.txt";
		if (MASTER(cr))
		{
			printf(MAKE_GREEN "Reading the flow velocity for AFM from %s..." RESET_COLOR "\n",fname_in_velocity);
		}

		int is_flow = read_flow_velocity(fh,fname_in_velocity);

		if (is_flow != fh->flow_type)
		{
			if(MASTER(cr)) printf(MAKE_RED "Inconsistency between flow_type in input and flow velocity field" RESET_COLOR "\n");
			exit(5);
		}
	}

	if (fh->S_function == moving_sphere)
	{
		char const *fname_vel = "vel_poiseuille.txt";
		if(MASTER(cr))
		{
			printf(MAKE_GREEN "Reading the flow velocity for Couette flow from %s..." RESET_COLOR "\n",fname_vel);
		}

		int is_flow = read_flow_velocity(fh,fname_vel);

		if (is_flow != fh->flow_type)
		{
			if(MASTER(cr)) printf(MAKE_RED "Inconsistency between flow_type in input and flow velocity field" RESET_COLOR "\n");
			exit(6);
		}

	}




    if(fh->scheme == One_Way)
        FH_init(fh, cr);

    for(int i = 0; i < fh->Ntot; i++)
    {
        fh->arr[i].S = 1;
        for(int d = 0; d < DIM; d++)
        {
            fh->arr[i].Sf[d] = 1;
        }
    }

    if(MASTER(cr))
    {
        if(fh->scheme == One_Way)
        {
            FH_equilibrate(fh);
            printf(MAKE_GREEN "FHMD: Initialization finished. Starting MD/FH solver...\n" RESET_COLOR "\n");
        }
        fflush(stdout);
    }

    return 1;   // Success
}
