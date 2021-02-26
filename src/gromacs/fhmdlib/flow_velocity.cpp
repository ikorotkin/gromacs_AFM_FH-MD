#include "data_structures.h"
#include "flow_velocity.h"
#include "macro.h"
#include <stdio.h>
#include <stdlib.h>
#include "gromacs/topology/mtop_util.h"
#include "gromacs/math/units.h"
//#include "threads.h"
//#define FLOW_VELOCITY_DEBUG


int read_flow_velocity(FHMD *fh, char const *fname_in_velocity)
{
	FILE *fvel;
	ivec ind;
	int counter = 0;
	FH_arrays *arr = fh->arr;
	char line[300];

	fvel = fopen(fname_in_velocity,"r");

	if (fvel == NULL)
	{
		printf(MAKE_GREEN "The file of flow velocity does not exist,set the flow velocity to zero" RESET_COLOR "\n");

		for(int k = 0; k < NZ; k++)
		    {
		        for(int j = 0; j < NY; j++)
		        {
		            for(int i = 0; i < NX; i++)
		            {
		                ASSIGN_IND(ind, i, j, k);

		                //for(int d = 0; d < DIM; d++)
		                //{
		                	clear_dvec(arr[C].u_fh_flow);
		                //}
		            }
		        }
		    }

		return 0;
	}

	else
	{
		while (fscanf(fvel,"%s",line)!=-1)
		{

			//skip the first line
			if(fgets(line,sizeof(line),fvel)==NULL)
			{
				printf(MAKE_RED "FHMD error: error when skipping the header of flow velocity file");
				exit(1);
			}

			for(int k = 0; k < NZ; k++)
			{
				for(int j = 0; j< NY; j++)
				{
					for(int i = 0; i< NX; i++)
					{
						ASSIGN_IND(ind, i, j, k);
						skip_coordinate(fvel);
						counter++;
						for (int d = 0; d < DIM; d++)
						{
							int ok = fscanf(fvel,"%le,",&arr[C].u_fh_flow[d]);
						}
					}
				}
			}
		}

		if (counter != NX*NY*NZ)
		{
			printf(MAKE_RED "The grids of flow velocity are not equal to grids of FH-MD" RESET_COLOR "\n");
			exit(0);
		}
	}

#ifdef FLOW_VELOCITY_DEBUG
	printf("U,\t V,\t W\t\n");
	for(int k = 0; k < NZ; k++)
	{
		for(int j = 0; j < NY; j++)
		{
			for (int i = 0; i < NX; i++)
			{
				ASSIGN_IND(ind, i, j, k);

				for (int d = 0; d < DIM; d++)
				{
					printf("%le,\t",arr[C].u_fh_flow[d]);
				}

				printf("\n");
			}

		}
	}
#endif

	return 1;
}



void skip_coordinate(FILE *fvel)
{
	double coor;

	for(int d = 0; d < DIM; d++)
	{
		int ok = fscanf(fvel,"%lf,",&coor);
	}

}


void compute_vx_instant_mpi (FHMD *fh, rvec x[], rvec v[], real mass[], int N_atoms, t_commrec *cr)
{
	 int nzbin = fh->nzbin;
	 int vbin,N;
	 double vx[nzbin];
	 double mass_bin[nzbin];
	 double vel_mass[nzbin];
	 rvec       	 *vel;

	 vel = fh->vel;

	 N = N_atoms;

	 //if(MASTER(cr))
	 //{

		 for (int i = 0; i < nzbin; i++)
	     {
			 vx[i] = 0;
	         mass_bin[i] = 0;
	         vel_mass[i] = 0;
	         //vx_m[i] = 0 ;
	     }

	 //}

	 for (int n = 0; n < N; n++)
	 {
		 vbin = (int)(x[n][2]/fh->box[2]*(double)(nzbin));
		 //vbin = (int)(x[i][1]/hhmd->box.y*(double)(vx_bins));
		 //vbin = find_ind(n,x,fh);
	     vx[vbin] += mass[n]*v[n][XX];
	     mass_bin[vbin] += mass[n];
	     vel_mass[vbin] += mass[n]*vel[n][XX];
	     //fprintf(fp,"%d        %.4f    %.4f\n",vbin,mass[n],v[n][0]);


	 }

	 for (int i = 0; i < nzbin; i++)
	 {
		 if (PAR(cr))
		 {
			 gmx_sumd(1,&vx[i],cr);
		     gmx_sumd(1,&mass_bin[i],cr);
		     gmx_sumd(1,&vel_mass[i],cr);
		 }

		 fh->vx[i] = vx[i];
		 fh->mass_bin[i] = mass_bin[i];
		 fh->vel_mass[i] = vel_mass[i];
	 }


}

void compute_average_vx_instant (FHMD *fh)
{
	for (int i = 0; i < fh->nzbin; i++)
	{
		fh->vx_m[i] = fh->vx[i]/fh->mass_bin[i];    //instant average velocity
		fh->vel_m[i] = fh->vel_mass[i]/fh->mass_bin[i];
		fh->sum_vx_m[i] = fh->sum_vx_m[i] + fh->vx_m[i];
		fh->sum_vel[i] = fh->sum_vel[i] + fh->vel_m[i];
	}
}


int find_ind(int n, const rvec x[],FHMD *fh)
{
	int ibinz, d = 2;
	dvec xn;
	double prd,invdelta;

	PBC(xn, x[n], fh->box);

		prd = fh->box[d];
		invdelta = fh->nzbin/prd;
		ibinz = static_cast<int> (xn[d] * invdelta);
		return ibinz;

}

void compute_tem (FHMD *fh, const rvec x[], rvec v[], real mass[], int N_atoms, t_commrec *cr, gmx_mtop_t *mtop)
{
	double tem;
	int N,res_nr;
	int ind;
	int counter;
	int ndf;
	//int num_cores;
	double up,lower;
	double ke_total;
	dvec xn;
	dvec ke;
	char    * atomname=NULL, *resname=NULL;

	N = N_atoms;
	up = 1.5;
	lower = 1;
	//num_cores = cr->nnodes;

	clear_dvec(ke);
	counter = 0;
	ndf = 0;
	ke_total = 0;
	tem = 0;

	for (int n = 0; n < N; n++)
	{
		if(PAR(cr) && DOMAINDECOMP(cr))
		{
			ind = cr->dd->gatindex[n];
		} else {
		    ind = n;
		}

		gmx_mtop_atominfo_global(mtop, ind, &atomname, &res_nr, &resname);

		if(!strcmp(resname, "SOL"))
		{
			PBC(xn, x[n], fh->box);
			if (xn[ZZ] >= lower && xn[ZZ] <= up)
			{
				counter ++;
				for (int d = 0; d < DIM; d++)
				{
					ke[d] +=  + 0.5*mass[n]*v[n][d]*v[n][d];
				}
			}
		}

		else
			continue;
	}

	ke_total = ke[XX] + ke [YY] + ke[ZZ];
	ndf = 3*counter - counter - 3;

	 if (PAR(cr))
	 {
		 gmx_sumd(1,&ke_total,cr);
		 gmx_sumi(1,&ndf,cr);
	 }

	 tem = 2*ke_total/(ndf*BOLTZ);
	 fh->AFM_tem = tem;
}




