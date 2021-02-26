#ifndef FHMD_DATA_STRUCTURES_H_
#define FHMD_DATA_STRUCTURES_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "params.h"
#include "gromacs/topology/topology.h"      /* GROMACS gmx_mtop_t  structure definition */
#include "gromacs/domdec/domdec_struct.h"   /* GROMACS gmx_domdec_t structure definition */
#include "gromacs/mdtypes/commrec.h"        /* GROMACS MPI definitions, 't_commrec', MASTER(), PAR(), etc. */
#include "gromacs/gmxlib/network.h"         /* GROMACS MPI functions, gmx_bcast(), gmx_sumf(), etc. */
#include "gromacs/math/vectypes.h"          /* GROMACS vector types: rvec, dvec, ivec, etc. */
#include "gromacs/math/vec.h"               /* GROMACS vector operations: copy_ivec, dvec_add, etc. */
#include "gromacs/mdlib/rbin.h"         /* for gmx_fhmd_global_stat struct*/

typedef struct FH_arrays                    /* FH/MD arrays */
{
    double      ro_md, ro_fh;               /* densities */
    double      inv_ro;                     /* inverse density: 1/ro_md */
    dvec        u_md, u_fh;                 /* velocities */
    dvec        uro_md;                     /* momentum */
    dvec        f_fh;                       /* FH force */
    dvec        alpha_term;                 /* alpha term for du/dt equation */
    dvec        beta_term;                  /* beta term for du/dt equation */

    double      delta_ro;                   /* delta of MD and FH densities */
    dvec        grad_ro;                    /* grad of density */
    matrix      alpha_u_grad;               /* preliminary alpha-term [u-index][grad-index] */

    double      S;                          /* S parameter in the FH cell centres */
    dvec        Sf;                         /* S parameter in the FH cell faces */

    double      p, pn;                      /* FH pressure */
    dvec        rof, rofn, pf, pfn;         /* FH flux variables */
    matrix      uf, ufn;                    /* FH flux velocities */
    double      ro_fh_n;                    /* FH density (new time layer) */
    dvec        u_fh_n;                     /* FH velocity (new time layer) */
    matrix      rans;                       /* FH random stress */

    double      ro_md_prime;                /* MD density from layer n */
    dvec        uro_md_prime;               /* MD momentum from layer n */

    double      ro_prime, ron_prime, ronn_prime;    /* density prime */
    double      ro_prime_b;                         /* density prime on the boundary */
    double      ro_star, ron_star;                  /* density star */
    dvec        m_prime, mn_prime, mnn_prime;       /* momentum prime */
    dvec        m_prime_b;                          /* momentum prime on the boundary */
    dvec        m_star, mn_star;                    /* momentum star */

    double      ro_md_s, ros_md, ropr_md;           /* sources of MD density */
    dvec        uro_md_s, uros_md, uropr_md;        /* sources of MD momentum */
    dvec        u_fh_flow;                          /* the flow velocity for AFM*/
} FH_arrays;


typedef struct FH_grid          /* Computational grid */
{
    dvec       *c;              /* FH cell centres coordinates */
    dvec       *n;              /* FH cell nodes coordinates */
    dvec       *h;              /* FH cell steps */
    double     *vol;            /* FH cell volume */
    double     *ivol;           /* 1/cellVolume */
    FHMD_CELL  *md;             /* The cell is in the FH_zone, hybrid_zone or boundary */
} FH_grid;


typedef struct MD_stat          /* Particle statistics */
{
    int         N, n;
    double      invN;

    double     *avg_rho_md_cell, *avg_rho_fh_cell;

    double      davg_rho_md,  davg_rho_fh;
    double      davg_rho2_md, davg_rho2_fh;
    dvec        davg_u_md,    davg_u_fh;
    dvec        davg_u2_md,   davg_u2_fh;

    double      avg_rho_md,   avg_rho_fh;
    double      std_rho_md,   std_rho_fh;
    dvec        avg_u_md,     avg_u_fh;
    dvec        std_u_md,     std_u_fh;
} MD_stat;

/*struct for AFM KE*/
typedef struct AFM_KE
{
	real 		Th;                 /* Temperature at half step per bin*/
	real 		T;                  /* Temperature at full step per bin*/
	tensor 		ekinh;              /* Kinetic energy at half step per bin*/
	tensor      ekinh_old;          /* Kinetic energy at old half step per bin*/
	tensor      ekinf;              /* Kinetic energy at full step per bin*/
	real        lambda;             /* Berendsen coupling lambda per bin*/
	double      nrdf;               /* number of degree freedom of each bin */
	double      S;                  /* parameter in the bin centres */
	double      reft;
	double      vx_m;
} AFM_KE;

/*struct for summation of AFM_KE on each core if parallel*/
typedef struct gmx_fhmd_global_stat
{
	t_bin *rb_bin;
	int   *itc0_bin;
	int   *itc1_bin;
} gmx_fhmd_global_stat;

/*struct for CFD flow*/
/*typedef struct FH_arrays_flow
{
	double u_fh_flow;
	double ro_fh_flow;
} FH_arrays_flow;*/

typedef struct FHMD
{
    FH_arrays  *arr;            /* FH/MD arrays */
    FH_grid     grid;           /* FH grid */
    MD_stat     stat;           /* Particle statistics */
    int        *ind;            /* FH cell number for each atom */
    ivec       *indv;           /* 3-component FH cell number for each atom (vector) */
    double     *mpi_linear;     /* Linear array to summarise MDFH arrays */

    FHMD_S      S_function;     /* S = const or S = S(x,y,z) - fixed or moving */
    int         scheme;         /* 0 - Pure MD, 1 - One-way coupling, 2 - Two-way coupling */
    double      S;              /* Parameter S (-1 - fixed sphere, -2 - moving sphere) */
    double      R1;             /* MD sphere radius for variable S, [0..1] */
    double      R2;             /* FH sphere radius for variable S, [0..1] */

    /*parameters for AFM
     *
     */
    double     z1;              /*upper limit of coordinate of upper MD zone in z direction  for variable S */
    double     z2;              /*upper limit of coordinate of upper buffer zone in z direction for variable S*/
    double     z_c_ball;

    double     zs;              /*Derived variables from z1 and z2 */

    int        nzbin;           /*number of the bins used in bufferzone to thermostat the temperature*/
    AFM_KE    *ke_data;         /*kinetic energy data of AFM */
    tensor    ekin;
    tensor    **ke_data_work_ekinh; 		/*working array for the ke_data, the first dimension is number of bins,
    										/*the second is T-coupling group */
    tensor    **ke_data_work_ekinh_old;     /*working array for the ke_data, the first dimension is number of bins,
    										/*the second is T-coupling group */
    double    *nrdf_mpi;
    dvec      *bin_c;                       /*the z coordinate of center of each bin*/
    double     gamma;                        /*factor to adjust the alpha and beta*/
    double     eta;                          /*factor to adjust the thermostat*/
    double     thermostat;                   /*flag to choose the thermostat, 0 - original gromacs thermostat,
                                              1 - thermostat according to the bin*/
    THERMOSTAT thermostat_function;        /*0 - original gromacs thermostat,1 - thermostat according to the bin*/
    double AFM_tem; 						   /*temperature for water on which rdf is calculated*/


    double      Smin;           /* Minimum S for variable S */
    double      Smax;           /* Maximum S for variable S */
    double      R12, R22, RS;   /* Derived variables from R1, R2, Smin, Smax */
    double      alpha;          /* Alpha parameter for dx/dt and du/dt equations, nm^2/ps */
    double      beta;           /* Beta parameter, nm^2/ps or ps^-1 depending on the scheme */
    double      gamma_x;        /* Gamma_x parameter (density fluctuations dissipator), ps^-1 */
    double      gamma_u;        /* Gamma_u parameter (velocity fluctuations dissipator), ps^-1 */
    double      eps_rho;        /* Eps_rho parameter for ro_prime FH equation (dissipator factor, 0 <= eps_rho <= 1) */
    double      eps_mom;        /* Eps_mom parameter for m_prime FH equation (dissipator factor, 0 <= eps_mom <= 1) */
    double      S_berendsen;    /* If S_berendsen >= 0, Berendsen thermostat works for S <= S_berendsen, otherwise factor (1-S^(-S_berendsen)) is applied */

    ivec        N, N_md;        /* Number of FH cells along each direction */
    ivec        N_shift;        /* N_shift = (N - N_md)/2 */
    dvec        box;            /* Box size */
    dvec        box05;          /* Half of box size */

    /*parameters for protein
     *
     *
     */
    dvec        protein_com;    /* Protein COM coordinates */
    double      protein_mass;   /* Protein mass */
    int         protein_N;      /* Number of atoms in the protein */

    int           prot_num;
    char          tem[800];
    char         **prot_name;

    dvec         *prot_com;
    double       *prot_mass;
    int          *prot_N;
    int           flow; /*flag for with flow or not*/
    FLOW_TYPE flow_type;

    double 		 *vx_m;           /*instant average velocity*/
    double       *vel_m;
    double       *sum_vx_m;

    double		 *vx;
    double       *mass_bin;
    rvec       	 *vel;           /* Real velocity of each particle (dx/dt) */
    double       *vel_mass;         /* momentum of real velocity */
    double       *sum_vel;
    /*
     *
     */

    double      box_volume;     /* Volume of the box, nm^3 */
    double      total_density;  /* Total density of the box, a.m.u./nm^3 */
    int         Ntot, Ntot_md;  /* Total number of FH cells */
    int         step_MD;        /* Current MD time step */
    int         Noutput;        /* Write arrays to files every Noutput MD time steps (0 - do not write) */

    FHMD_EOS    eos;            /* Equation of state */
    int         FH_EOS;         /* EOS: 0 - Liquid Argon, 1 - SPC/E water */
    int         FH_step;        /* dt_FH = FH_step * dt_MD */
    int         FH_equil;       /* Number of time steps for the FH model equilibration */
    double      FH_dens;        /* FH mean density */
    double      FH_temp;        /* FH mean temperature */
    double      FH_blend;       /* FH Blending: -1 - dynamic, or define static blending parameter (0.0 = Central Diff., 1.0 = Classic CABARET) */
    double      dt_FH;          /* FH time step */
    double      std_rho;        /* Analytical STD of density */
    double      std_u;          /* Analytical STD of velocity */
} FHMD;

#endif /* FHMD_DATA_STRUCTURES_H_ */
