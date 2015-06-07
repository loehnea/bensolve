/*
This file is part of BENSOLVE - VLP solver

Copyright (C) 2014-2015 Andreas Löhne and Benjamin Weißing

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (see the reference manual). If not,
see <http://www.gnu.org/licenses/>
*/

#include <sys/time.h>	// gettimeofday()
#include <assert.h>		// assert
#include <string.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>

#include "bslv_lists.h"
#include "bslv_poly.h"
#include "bslv_lp.h"
#include "bslv_vlp.h"
#include "bslv_algs.h"


extern struct timeval t_end;

void init_image (poly_args *args, size_t dim, size_t dim_primg_prml, size_t dim_primg_dl, void (*dual_v2primal_h)(double *, int, double *))
{
	poly__set_default_args(args, dim);
	args->eps= EPS_POLY; 
	args->dim_primg_prml=dim_primg_prml; 
	args->dim_primg_dl=dim_primg_dl;
	args->primalV2dualH = NULL;
	args->dualV2primalH = dual_v2primal_h;
	poly__initialise (args);
}

void poly_output(poly_args *args, const opttype *opt, int stdout_at_msg_lev, pre_img_type pre_img, swap_type swap, const char *str1, const char *str2, const char *str3)
{
	polytope *primal, *dual;
	if (swap==NO_SWAP)
	{
		primal = &args->primal;
		dual = &args->dual;
	}
	else
	{
		assert(swap == SWAP);
		primal = &args->dual;
		dual = &args->primal;
	}

	char filename[strlen(opt->filename)+MAX_STR_LNGTH+strlen(str3)+1];

	for (size_t k=0; k<dual->cnt; k++)
		if (IS_ELEM(dual->used, k))
			ST_BT(dual->sltn, k);
	for (size_t k=0; k<primal->cnt; k++)
		if (IS_ELEM(primal->used, k))
			ST_BT(primal->sltn, k);

	permutation prm,prm_dual;
	poly__initialise_permutation (primal,&prm);
	poly__initialise_permutation (dual,&prm_dual);

	if (opt->message_level >= stdout_at_msg_lev)
	{
		printf("%s",str1);
		poly__vrtx2file (primal,&prm,NULL,opt->format==FORMAT_LONG?FORMAT_LONG_STR:FORMAT_SHORT_STR);
		printf("%s",str2);
		poly__vrtx2file (dual,&prm_dual,NULL,opt->format==FORMAT_LONG?FORMAT_LONG_STR:FORMAT_SHORT_STR);
	}

	// upper image
	strcpy(filename,opt->filename);
	strcat(filename, IMG_P_STR);
	strcat(filename, str3);
	poly__vrtx2file (primal,&prm,filename,opt->format==FORMAT_SHORT?FORMAT_SHORT_STR:FORMAT_LONG_STR);

	// lower image
	strcpy(filename,opt->filename);
	strcat(filename, IMG_D_STR);
	strcat(filename, str3);
	poly__vrtx2file (dual,&prm_dual,filename,opt->format==FORMAT_SHORT?FORMAT_SHORT_STR:FORMAT_LONG_STR);

	// adjacency list of upper image (row # = vrtx #)
	strcpy(filename,opt->filename);
	strcat(filename, ADJ_P_STR);
	strcat(filename, str3);
	poly__adj2file (primal, &prm, filename, NULL);

	// adjacency list of lower image (row # = vrtx #)
	strcpy(filename,opt->filename);
	strcat(filename, ADJ_D_STR);
	strcat(filename, str3);
	poly__adj2file (dual, &prm_dual, filename, NULL);

	//facet-vertex incidence list of upper image (row # = facet #)
	strcpy(filename,opt->filename);
	strcat(filename, INC_P_STR);
	strcat(filename, str3);
	poly__inc2file (primal, &prm, &prm_dual, filename,NULL);

	//facet-vertex incidence list of lower image (row # = facet #)
	strcpy(filename,opt->filename);
	strcat(filename, INC_D_STR);
	strcat(filename, str3);
	poly__inc2file (dual, &prm_dual, &prm, filename,NULL);

	if (pre_img == PRE_IMG_ON)
	{
		// pre-image (solution)
		strcpy(filename,opt->filename);
		strcat(filename, PRE_IMG_P_STR);
		strcat(filename, str3);
		if (pre_img==PRE_IMG_ON)
			poly__primg2file (primal, &prm, filename,opt->format==FORMAT_SHORT?FORMAT_SHORT_STR:FORMAT_LONG_STR);
		else
			string_fprint(filename,"Solution (pre-image) was not stored. Use option -s\n");

		// dual pre-image (dual solution)
		strcpy(filename,opt->filename);
		strcat(filename, PRE_IMG_D_STR);
		strcat(filename, str3);	
		if (pre_img==PRE_IMG_ON)
			poly__primg2file (dual, &prm_dual, filename,opt->format==FORMAT_SHORT?FORMAT_SHORT_STR:FORMAT_LONG_STR);
		else
			string_fprint(filename,"Solution (pre-image) was not stored. Use option -s\n");
	}
	poly__kill_permutation (&prm_dual);
	poly__kill_permutation (&prm);
}

void poly_count(poly_args *args, soltype *sol, swap_type swap)
{
	polytope *primal, *dual;
	if (swap == NO_SWAP)
	{
		primal = &args->primal;
		dual = &args->dual;
	}
	else
	{
		assert(swap == SWAP);
		primal = &args->dual;
		dual = &args->primal;
	}
	sol->pp=0;
	sol->pp_dir=0;
	for (size_t k=0; k<primal->cnt; k++)
	{	
		if(IS_ELEM(primal->used,k))
		{	
			if (IS_ELEM(primal->ideal,k))
				sol->pp_dir++;
			else
				sol->pp++;
		}
	}
	sol->dd=0;
	sol->dd_dir=0;
	for (size_t k=0; k<dual->cnt; k++)
	{
		if(IS_ELEM(dual->used,k))
		{	
			if (IS_ELEM(dual->ideal,k))
				sol->dd_dir++;	
			else
				sol->dd++;
		}
	}
}

void poly_chop(poly_args *args)
{
	for(size_t i=0; i<args->primal.cnt; i++)
		if(IS_ELEM(args->primal.used,i))
	{
		for (size_t j=0;j<args->dim;j++)
			if (ABS(args->primal.data[args->dim*i + j])<EPS_OUTPUT_CHOP)
				args->primal.data[args->dim*i + j] = 0.0;
		for (size_t j=0;j<args->primal.dim_primg;j++)
			if (ABS(args->primal.data_primg[args->primal.dim_primg*i + j])<EPS_OUTPUT_CHOP)
				args->primal.data_primg[args->primal.dim_primg*i + j] = 0.0;
	}	
	for(size_t i=0; i<args->dual.cnt; i++)
		if(IS_ELEM(args->dual.used,i))
	{
		for (size_t j=0;j<args->dim;j++)
			if (ABS(args->dual.data[args->dim*i + j])<EPS_OUTPUT_CHOP)
				args->dual.data[args->dim*i + j] = 0.0;
		for (size_t j=0;j<args->dual.dim_primg;j++)
			if (ABS(args->dual.data_primg[args->dual.dim_primg*i + j])<EPS_OUTPUT_CHOP)
				args->dual.data_primg[args->dual.dim_primg*i + j] = 0.0;
	}
}

void poly_minus(poly_args *args, size_t idx_p_1, size_t idx_p_2, size_t idx_d_1, size_t idx_d_2)
{
	for(size_t i=0; i<args->primal.cnt; i++)
		if(IS_ELEM(args->primal.used,i))
			for (size_t j=idx_p_1-1;j<idx_p_2;j++)
				args->primal.data[args->dim*i + j] = -args->primal.data[args->dim*i + j];
	for(size_t i=0; i<args->dual.cnt; i++)
		if(IS_ELEM(args->dual.used,i))
			for (size_t j=idx_d_1-1;j<idx_d_2;j++)
				args->dual.data[args->dim*i + j] = -args->dual.data[args->dim*i + j];
}

// transform output for maximization problems and the case of negative c_q
void poly_trans_primal(const vlptype *vlp, const soltype *sol, const opttype *opt, poly_args *args)
{
	if (sol->c_dir==C_DIR_POS && vlp->optdir==-1)	// max / c_q>0
		poly_minus(args,1,vlp->q,vlp->q,vlp->q);	// y -> -y ; y^*_q -> -y^*_q
	if (sol->c_dir==C_DIR_NEG && vlp->optdir==1)	// min / c_q<0
		poly_minus(args,1,vlp->q,1,0);				// y -> -y
	if (sol->c_dir==C_DIR_NEG && vlp->optdir==-1)	// max / c_q<0
		poly_minus(args,1,0,vlp->q,vlp->q);			// y^*_q -> -y^*_q
}

// transform output for maximization problems and the case of negative c_q
void poly_trans_dual(const vlptype *vlp, const soltype *sol, const opttype *opt, poly_args *args)
{
	if (sol->c_dir==C_DIR_POS && vlp->optdir==-1)	// max / c_q>0
		poly_minus(args,vlp->q,vlp->q,1,vlp->q);	// y -> -y ; y^*_q -> -y^*_q
	if (sol->c_dir==C_DIR_NEG && vlp->optdir==1)	// min / c_q<0
		poly_minus(args,1,0,1,vlp->q);				// y -> -y
	if (sol->c_dir==C_DIR_NEG && vlp->optdir==-1)	// max / c_q<0
		poly_minus(args,vlp->q,vlp->q,1,0);			// y^*_q -> -y^*_q
}

void poly_normalize_dir(poly_args *args)
{
	for(size_t i=0; i<args->primal.cnt; i++)
	{
		if(IS_ELEM(args->primal.used,i)&&IS_ELEM(args->primal.ideal,i))
		{
			double max = 0;
			for (size_t j=0;j<args->dim;j++)
				if (ABS(args->primal.data[args->dim*i + j]) > max)
					max = ABS(args->primal.data[args->dim*i + j]);
			if (max > 1e-9)
				for (size_t j=0;j<args->dim;j++)
					args->primal.data[args->dim*i + j] /= max;
			else
				for (size_t j=0;j<args->dim;j++)
					args->primal.data[args->dim*i + j]=0;
		}
	}

	for(size_t i=0; i<args->dual.cnt; i++)
	{		
		if(IS_ELEM(args->dual.used,i)&&IS_ELEM(args->dual.ideal,i))
		{
			double max = 0;
			for (size_t j=0;j<args->dim;j++)
				if (ABS(args->dual.data[args->dim*i + j]) > max)
					max = ABS(args->dual.data[args->dim*i + j]);
			if (max > 1e-9)
				for (size_t j=0;j<args->dim;j++)
					args->dual.data[args->dim*i + j] /= max;
			else
				for (size_t j=0;j<args->dim;j++)
					args->dual.data[args->dim*i + j]=0;
		}
	}
}

struct 
{
	size_t dim;
	double *ip;
} fnc_prmtr;

static void lowerV2upperH (double *vrtx_lower, int is_dir, double *hp_upper)
{
	if (is_dir)
	{
		for (size_t j=0; j<fnc_prmtr.dim; j++)
			hp_upper[j]=0;
		hp_upper[fnc_prmtr.dim]=-1.0;
	}
	else
	{
		hp_upper[fnc_prmtr.dim-1]=1.0;
		for (size_t j=0; j<fnc_prmtr.dim-1; j++)
		{
			hp_upper[j]=vrtx_lower[j];
			hp_upper[fnc_prmtr.dim-1]-=fnc_prmtr.ip[j]*hp_upper[j];
		}
		hp_upper[fnc_prmtr.dim]=vrtx_lower[fnc_prmtr.dim-1];
	}
}

static void upperV2lowerH (double *vrtx_upper, int is_dir, double *hp_lower)
{
	hp_lower[fnc_prmtr.dim-1]=is_dir?0:-1.0;
	for (size_t j=0; j<fnc_prmtr.dim-1; j++)
		hp_lower[j]=vrtx_upper[j]-vrtx_upper[fnc_prmtr.dim-1]*fnc_prmtr.ip[j];
	hp_lower[fnc_prmtr.dim]=-vrtx_upper[fnc_prmtr.dim-1];
}

static void trnsfrm_plot (double *dual_vrtx, int is_dir, double *hp)
{
	hp[0]=-dual_vrtx[0];
	hp[1]=-dual_vrtx[1];
	hp[2]=-1+dual_vrtx[0]+dual_vrtx[1];
	hp[3]=-dual_vrtx[2];
}

static void trnsfrm_plot_dual (double *vrtx, int is_dir, double *hp)
{
	hp[0]=vrtx[0];
	hp[1]=vrtx[1];
	hp[2]=1-vrtx[0]-vrtx[1];
	hp[3]=vrtx[2];
}

int cone_vertenum(double **prim, lp_idx *n_prim, double **dual, lp_idx *n_dual, double *prim_in, const size_t n_prim_in, const size_t dim, const opttype *opt, cone_out_type output, swap_type swap)
{
	poly_args cone_poly;

	poly__set_default_args(&cone_poly, dim);
	cone_poly.eps= EPS_POLY;
	poly__initialise(&cone_poly);
	UNST_BT (cone_poly.dual.ideal,0);
	cone_poly.dual.data[dim-1] = 0;
		
	for (size_t k=0; k<n_prim_in; k++)
	{
		for (size_t j=0; j<dim; j++)
		{
			cone_poly.val[j]=prim_in[j*n_prim_in+k];
		}
		cone_poly.ideal=1;
		poly__add_vrtx(&cone_poly);
	}
	if (poly__intl_apprx (&cone_poly)==EXIT_FAILURE) // vertex enumeration
		return EXIT_FAILURE;

	// determine number of non-redundant generators of primal cone
	*n_prim=0;
	for (size_t k=0; k<cone_poly.dual.cnt; k++)
	{
		if(IS_ELEM(cone_poly.dual.used,k)&&IS_ELEM(cone_poly.dual.ideal,k))
			++*n_prim;
	}

	// store non-redundant generators of primal cone
	if (*prim != NULL) free(*prim);
	*prim = malloc(*n_prim*dim*sizeof(double));
	for (size_t k=0,i=0; k<cone_poly.dual.cnt; k++)
	{	
		if(IS_ELEM(cone_poly.dual.used,k)&&IS_ELEM(cone_poly.dual.ideal,k))
		{
			for (size_t j=0; j<dim; j++)
			{
				(*prim)[j*(*n_prim)+i]=cone_poly.dual.data[k*dim+j];
			}
			i++;
		}
	}

	// determine number of generators of the dual cone
	*n_dual=0;
	for (size_t k=0; k<cone_poly.primal.cnt; k++)
		if(IS_ELEM(cone_poly.primal.used,k)&&IS_ELEM(cone_poly.primal.ideal,k))
			++*n_dual;

	// store generators of polar cone
	if (*dual != NULL) free(*dual);
	*dual = malloc(*n_dual*dim*sizeof(double));
	size_t i = 0;
	while (!poly__get_vrtx(&cone_poly))
	{
		ST_BT(cone_poly.primal.sltn,cone_poly.idx);
		if (!IS_ELEM(cone_poly.primal.ideal,cone_poly.idx)) // rule out vertices
			continue;
		for(size_t j=0; j<dim; j++)
			(*dual)[j*(*n_dual)+i]=cone_poly.val[j];
		i++;
	}

	poly_chop(&cone_poly);
	poly_normalize_dir(&cone_poly);
	poly__update_adjacence (&cone_poly.dual);

	// write cone to file
	if (output==CONE_OUT_ON)
		poly_output(&cone_poly,opt,2,PRE_IMG_OFF,swap,CONE_P_STR,CONE_D_STR,CONE_ENDING_STR);
	if (POLY_TEST)
		poly__polyck (&cone_poly);
	poly__kill(&cone_poly);
	return EXIT_SUCCESS;
}

void poly_plot_primal(const vlptype *vlp, const soltype *sol, const opttype *opt, poly_args *args)
{
	assert(vlp->q==3);

	// compute normal direction "eta" (with sum-norm = 1) of cutting hyperplane
	// from sum of generators (2-norm = 1) of dual of recession cone of upper image
	double eta[] = {0.0, 0.0, 0.0};
	{
		double prod, sum=0.0;
		for (size_t k=0; k<sol->r; k++)
		{
			prod=0.0;
			for (size_t j=0; j<3; j++)
				prod+=sol->R[j*sol->r+k]*sol->R[j*sol->r+k];
			for(size_t j=0; j<3; j++)
				eta[j]+=sol->R[j*sol->r+k]/sqrt(prod);
		}
		for(size_t j=0; j<3; j++)
			sum+=eta[j];
		for(size_t j=0; j<3; j++)
			eta[j]/=sum;
	}

	// min and max rhs of cutting hyperplane taken over all vertices
	double alph_max=-DBL_MAX;
	double alph_min= DBL_MAX;
	polytope *poly= &args->primal;
	for (size_t j=0; j<poly->cnt; j++)
	if (IS_ELEM(poly->used,j)&&!IS_ELEM(poly->ideal,j)){
		double tmp=0;
		for (size_t k=0; k<3; k++)
			tmp+=eta[k]*poly->data[3*j+k];
		if (tmp>alph_max)
			alph_max=tmp;
		else if (tmp<alph_min)
			alph_min=tmp;
	}

	// set transformation map
	args->dualV2primalH=trnsfrm_plot;
	// search for a direction, which is used as outer (ideal) vertex for the cut
	for (args->idx=0;args->idx<poly->cnt;args->idx++)
		if (IS_ELEM(poly->used,args->idx)&&IS_ELEM(poly->ideal,args->idx))
			break;
	

	args->val[0]=eta[0];
	args->val[1]=eta[1];
	args->val[2]=alph_max+PRIMAL_PLOT_CUT_SHIFT*(alph_max-alph_min<10e-8?1.0:alph_max-alph_min);
	args->ideal=0;
	poly__add_vrtx (args);

	// compute width of resulting polytope for scaling
	double wdth_max[3]={-DBL_MAX,-DBL_MAX,-DBL_MAX};
	double wdth_min[3]={ DBL_MAX, DBL_MAX, DBL_MAX};
	for (size_t j=0;j<poly->cnt;j++)
	{	
		if (IS_ELEM(poly->used,j)&&!IS_ELEM(poly->ideal,j)){
			for (size_t k=0; k<3; k++)
			{	
				if (wdth_max[k]<poly->data[j*vlp->q+k])
					wdth_max[k]=poly->data[j*vlp->q+k];
				else if (wdth_min[k]>poly->data[j*vlp->q+k])
					wdth_min[k]=poly->data[j*vlp->q+k];
			}
		}
	}	
	// transform output 
	poly_trans_primal(vlp,sol,opt,args);
	// write polytope to off-file, write transformation to inst-file
	{
		char filename[strlen(opt->filename)+6+1];
		strcpy(filename,opt->filename);
		strcat(filename, "_p.off");
		poly__plot(poly,filename);
		
		char filename2[strlen(opt->filename)+7+1];
		strcpy(filename2,opt->filename);
		char *ptr=strrchr(filename,'/');
		if (ptr)
			ptr++;
		else
			ptr=filename;
		strcat(filename2, "_p.inst");
		FILE *strm=fopen (filename2,"w");
		fprintf(strm, "INST\ngeom < %s",ptr);
		fprintf(strm, "\ntransform {%f 0 0 0 0 %f 0 0 0 0 %f 0 0 0 0 1.0}\n", 1.0/(wdth_max[0]-wdth_min[0]),1.0/(wdth_max[1]-wdth_min[1]),1.0/(wdth_max[2]-wdth_min[2]));
		fclose(strm);
	}
}

void poly_plot_dual(const vlptype *vlp, const soltype *sol, const opttype *opt, poly_args *args)
{
	const double golden_ratio=(1+sqrt(5))/2;
	double wdth_max[3]={-DBL_MAX,-DBL_MAX,-DBL_MAX};
	double wdth_min[3]={ DBL_MAX, DBL_MAX, DBL_MAX};
	polytope *poly = &args->primal;
	for (size_t j=0;j<poly->cnt;j++)
	{
		if (IS_ELEM(poly->used,j)&&!IS_ELEM(poly->ideal,j)){
			for (size_t k=0; k<3; k++)
			{	
				if (wdth_max[k]<poly->data[j*vlp->q+k])
					wdth_max[k]=poly->data[j*vlp->q+k];
				if (wdth_min[k]>poly->data[j*vlp->q+k])
					wdth_min[k]=poly->data[j*vlp->q+k];
			}
		}
	}

	double hght=wdth_max[2]-wdth_min[2]<1e-8?1.0:(1+DUAL_PLOT_CUT_SHIFT)*(wdth_max[2]-wdth_min[2]);
	
	// cut polyhedron and obtain a polytope
	{
		args->dualV2primalH=&trnsfrm_plot_dual;
		// search for a direction (there should be exactly one),
		// which is used as outer (ideal) vertex for the cut
		for (args->idx=0;args->idx<args->primal.cnt;args->idx++)
			if (IS_ELEM(poly->used,args->idx)&&IS_ELEM(poly->ideal,args->idx))
				break;
		// determine rhs of cutting hyperplane
		args->val[0]=0;
		args->val[1]=0;
		args->val[2]=wdth_max[2]-hght;
		args->ideal=0;
		poly__add_vrtx (args);
	}

	// transform output for maximization problems and the case of negative c_q
	poly_trans_dual(vlp,sol,opt,args);
	// write polytope to off-file, write transformation to inst-file 
	{
		char filename[strlen(opt->filename)+6+1];
		strcpy(filename,opt->filename);
		strcat(filename, "_d.off");
		poly__plot (&args->primal,filename);
		
		char filename2[strlen(opt->filename)+7+1];
		strcpy(filename2,opt->filename);
		char *ptr=strrchr(filename,'/');
		if (ptr)
			ptr++;
		else
			ptr=filename;
		strcat(filename2, "_d.inst");
		FILE *strm=fopen (filename2,"w");
		fprintf(strm, "INST\ngeom < %s",ptr);
		fprintf(strm, "\ntransform {%f 0 0 0 0 %f 0 0 0 0 %f 0 0 0 0 1.0}\n", 1.0/(wdth_max[0]-wdth_min[0]), 1.0/(wdth_max[1]-wdth_min[1]), 1.0/hght/golden_ratio);
		fclose(strm);
	}
}


/* INIT_P2
 *
 *	initialize the lp for the problem
 *		(P_2(v))  min z 	s.t.
 *		rls <= A x <= lhs, lb <= x <= ub (vlp constraints)
 *		- P x +  y  == 0
 *		Z'y - Z'c z <= Z'v
 *		eta'y       <= 1	(if type == HOMOGENEOUS), infty (if type == INHOMOGENEOUS)
 *		where
 *		type = HOMOGENEOUS / INHOMOGENEOUS (vlp costraints of hom./inhom. problem)
 *		assuming that matrix A is set explicitly
 */
static void init_P2(soltype *const sol, const vlptype *vlp, lp_hom_type type)
{
	double *ZR;
	int p;
	/* set ZR=Z for homogeneous and ZR=R for inhomogeneous problem
	 * Z ... generating vectors of dual cone of ordering cone
	 * R ... generating vectors of dual cone of recession cone of upper image
	 */
	if (type == HOMOGENEOUS)
	{
		p = sol->p;
		ZR = sol->Z;
	}
	else
	{
		assert(type == INHOMOGENEOUS);
		p = sol->r;
		ZR = sol->R;
	}
	lp_update_extra_coeffs(p + 1, 1); // update lp dimension
	{	// set lp extra block of constraint matrix: A2 = (0 Z^T -Z^T c; 0 eta^T 0), where -Z^T c = (-1,...,-1)^T
		double *A2 = (double *) calloc((p + 1) * (vlp->n + vlp->q + 1), sizeof(double));
		for (lp_idx i = 0; i < p; i++)
		{
			for(lp_idx j = 0; j < vlp->q; j++)
				A2[i * (vlp->n + vlp->q + 1) + vlp->n + j] = ZR[j * p + i];
			A2[i * (vlp->n + vlp->q + 1) + vlp->n + vlp->q] = -1.0;
		}
		for(lp_idx j = 0; j < vlp->q; j++)
			A2[p * (vlp->n + vlp->q + 1) + vlp->n + j] = sol->eta[j];
		{
			list1d *rowlist = list1d_alloc(vlp->n + vlp->q + 1);
			list1d_init_idx(rowlist, 1);
			for(lp_idx i = 0; i < p + 1; i++)
			{
				vector_to_list1d(rowlist, A2 + i * (vlp->n + vlp->q + 1), vlp->n + vlp->q + 1);
				lp_set_mat_row(0, rowlist, vlp->m + vlp->q + 1 + i);
			}
			list1d_free(rowlist);
		}
		free(A2);
	}
	{ // set lp objective coefficients
		lp_clear_obj_coeffs(0);
		list1d *obj = list1d_calloc(1);
		obj->idx[0] = vlp->n + vlp->q + 1; obj->data[0] = 1;
		lp_set_obj_coeffs(0, obj);
		list1d_free(obj);
	}
	// set row bounds of lp based on vlp constriants of (in)homogeneous problem
	if (type == HOMOGENEOUS)
		lp_set_rows_hom(0, vlp->rows);
	else
	{
		assert(type == INHOMOGENEOUS);
		lp_set_rows(0, vlp->rows);
	}
	{	// set remaining row bounds part 1
		boundlist *rows = boundlist_calloc( vlp->q, 's');
		boundlist_init_idx(rows, vlp->m + 1);
		lp_set_rows(0, rows);
		boundlist_free(rows);
	}
	{	// set remaining row bounds part 2
		boundlist *rows = boundlist_calloc( p + 1, 'u');
		boundlist_init_idx(rows, vlp->m + vlp->q + 1);
		if (type == HOMOGENEOUS)
			rows->ub[p] = 1;
		else
		{
			assert(type == INHOMOGENEOUS);
			rows->type[p] = 'f';
		}
		lp_set_rows(0, rows);
		boundlist_free(rows);
	}
	// set column bounds of lp based on vlp constriants of (in)homogeneous problem
	if (type == HOMOGENEOUS)
		lp_set_cols_hom(0, vlp->cols);
	else
	{
		assert(type == INHOMOGENEOUS);
		lp_set_cols(0, vlp->cols);
	}
	{	// set y_1,...,y_q, z to be a free variables
		boundlist *cols = boundlist_calloc(vlp->q + 1, 'f');
		boundlist_init_idx(cols, vlp->n + 1);
		lp_set_cols(0, cols);
		boundlist_free(cols);
	}
}

/* PHASE0
 *
 *	Phase 0 of Benson Algorithm computes:
 *	eta ∈ int(D^* + K) with eta' c == 1, see [1] Section 5.5 (Algorithm 3)
 *	here we use a slightly (w.r.t. LP duality) modified variant
 *	here k == q-1, where k is the dimension used in [1]
 */
void phase0(soltype *const sol, const vlptype *vlp, const opttype *opt)
{
	init_P2(sol, vlp, HOMOGENEOUS);
	lp_set_options(&(opt->lp), PHASE0);

	double tmp1, tmp2;
	double z[vlp->q - 1];					// to store the vector z (see [1])
	double V[(vlp->q - 1)*(vlp->q - 1)];	// to store the vectors V(0), ... ,V(q-2) in R^{q-1} (see [1])
	double C[(vlp->q - 1) * (vlp->q - 1)];	// to store the vectors C(0), ... ,C(q-2) in R^{q-1} (see [1])
	double ww_red[vlp->q - 1];
	boundlist* rows = boundlist_calloc(sol->p, 'u');	// to reset the upper row bounds Z'y in P_2(y)
	boundlist_init_idx(rows, vlp->m + vlp->q + 1);		// upper bounds for indices m+1 to m+p
	if (opt->message_level >= 3) printf("solve lp\n");
	lp_status_type lp_status = lp_solve(0); 
	if( lp_status == LP_UNBOUNDED)
	{		
		sol->status=VLP_UNBOUNDED;
		boundlist_free(rows);
		return;
	}
	else	
		assert(lp_status == LP_OPTIMAL);
	lp_dual_solution_rows(0, z, vlp->m + 1, vlp->q - 1,1);

	for (size_t i = 0; i < vlp->q - 1; i++)
	{
		// store in column C(i) a non-zero vector being orthogonal to C(0) ... C(i-1)
		orthogonal_vector(C, vlp->q - 1, i);
		// solve P_2( [C(i);0] )
		// rows->ub = Z' * C(i)
		for (size_t j=0; j<sol->p; j++)
		{	
			rows->ub[j]=0;
			for (size_t k=0; k<vlp->q-1; k++)
				rows->ub[j]+=sol->Z[k*sol->p+j]*C[k*(vlp->q - 1)+i];
		}
		lp_set_rows(0, rows);
		if (opt->message_level >= 3) printf("solve lp\n");
		assert(lp_solve(0) == LP_OPTIMAL);
			
		// V(i) = w_red - z
		lp_dual_solution_rows(0, ww_red, vlp->m + 1, vlp->q - 1,1);
		// ww_red = ww_red - z
		for (size_t k=0; k<vlp->q-1; k++)
			ww_red[k]-=z[k];
		
		// V(i) = ww_red
		for (size_t k=0;k<vlp->q-1;k++)
			V[k*(vlp->q - 1)+i]=ww_red[k];
		
		/* if |<C(i), V(i)>| < eps
		 *		solve P_2( -[C(i);0] )
		 *		V(i) = w_red - z
		 */
		tmp1=0;
		for (size_t k=0;k<vlp->q-1;k++)
			tmp1+=C[k*(vlp->q - 1)+i]*V[k*(vlp->q - 1)+i];
		
		if (ABS(tmp1) < opt->eps_phase0)
		{
			// rows->ub = - Z' * C(i)
			for (size_t j=0;j<sol->p;j++)
			{
				rows->ub[j]=0;
				for (size_t k=0;k<vlp->q-1;k++)
					rows->ub[j]-=sol->Z[k*sol->p+j]*C[k*(vlp->q - 1)+i];
			}
			lp_set_rows(0, rows);
			if (opt->message_level >= 3) printf("solve lp\n");
			assert(lp_solve(0) == LP_OPTIMAL);
			lp_dual_solution_rows(0, ww_red, vlp->m + 1, vlp->q - 1, 1);
			for (size_t k=0;k<vlp->q-1;k++)
				ww_red[k]-=z[k];
			
			for (size_t k=0;k<vlp->q-1;k++)
				V[k*(vlp->q - 1)+i]=ww_red[k];
		}

		//if |<C(i), V(i)>| < eps -> upper image has no vertex
		tmp1=0;
		for (size_t k=0;k<vlp->q-1;k++)
			tmp1+=C[k*(vlp->q - 1)+i]*V[k*(vlp->q - 1)+i];
		
		if (ABS(tmp1) < opt->eps_phase0)
		{
			sol->status = VLP_NOVERTEX;
			boundlist_free(rows);
			return;
		}
		// C(i) = V(i) - sum_{j=0}^{i-1}  <C(j),V(i)> / <C(j),C(j)> * C(j) (use storage of eta)
		for (size_t k=0;k<vlp->q;k++)
			sol->eta[k]=0.0;

		for (size_t j = 0; j < i; j++) 
		{
			tmp1=0;
			for (size_t k=0;k<vlp->q-1;k++)
				tmp1+=C[k*(vlp->q - 1)+j]*V[k*(vlp->q - 1)+i];
			tmp2=0;
			for (size_t k=0;k<vlp->q-1;k++)
				tmp2+=C[k*(vlp->q - 1)+j]*C[k*(vlp->q - 1)+j];
			for (size_t k=0;k<vlp->q-1;k++)
				sol->eta[k]-=tmp1/tmp2*C[k*(vlp->q - 1)+j];
		}		 
		for (size_t k=0;k<vlp->q-1;k++)
			sol->eta[k]+=V[k*(vlp->q - 1)+i];
		for (size_t k=0;k<vlp->q-1;k++)
			C[k*(vlp->q-1)+i]=sol->eta[k];
	}

	/* compute mean of 0, V(0), ... , V(q-2)
	 * add z
	 * determine eta(q-1) such that c' eta == 1
	 */
	for (size_t k=0; k<vlp->q-1; k++)
		sol->eta[k]=V[k*(vlp->q-1)+0];
	for (size_t i=1; i<vlp->q-1; i++)
		for (size_t k=0; k<vlp->q-1; k++)
			sol->eta[k]+=V[k*(vlp->q-1)+i];
	for (size_t k=0; k<vlp->q-1; k++)
		sol->eta[k]/=vlp->q;
	for (size_t k=0; k<vlp->q-1; k++)
		sol->eta[k]+=z[k];
	sol->eta[vlp->q-1]=1.0;
	for (size_t k=0; k<vlp->q-1; k++)
		sol->eta[vlp->q-1]-=sol->c[k]*sol->eta[k];
	boundlist_free(rows);
} // end of phase0

/*
 *
 *
 * PHASE 1 -- PRIMAL
 * 
 *
 *
 */

void phase1_primal(soltype *const sol, const vlptype *vlp, const opttype *opt)
{
	poly_args upper_image;

	fnc_prmtr.dim=sol->q;
	fnc_prmtr.ip=sol->c;

	init_image (&upper_image, vlp->q, 0, 0, &lowerV2upperH);

	init_P2(sol, vlp, HOMOGENEOUS); // assume that appropriate eta is stored in sol->eta
	lp_set_options(&(opt->lp), PHASE1_PRIMAL);

	/* PHASE 1 -- PRIMAL -- PART 1
	 *
	 * determine first outer approximation of upper image
	 *
	 */
	{
		upper_image.ideal=0; // mark y^* (dual vertex to be added) as a point

		boundlist *rows = boundlist_calloc(sol->p, 'f');
		boundlist_init_idx(rows, vlp->m + vlp->q + 1);
		for (lp_idx j = 0; j < sol->p; j++)
		{
			rows->type[j] = 'u';
			lp_set_rows(0, rows);
			rows->type[j] = 'f';
			for (size_t k=0; k<vlp->q-1; k++)
				upper_image.val[k]=sol->Z[k*sol->p+j];
			if (opt->message_level >= 3) printf("initialization - solve lp\n");
			assert(lp_solve(0) == LP_OPTIMAL);
			upper_image.val[upper_image.dim-1]=lp_obj_val(0);
			poly__add_vrtx (&upper_image);
		}
		boundlist_free(rows);

		poly__intl_apprx (&upper_image);
	}
	
	/* PHASE 1 -- PRIMAL -- PART 2
	 *
	 * main loop
	 *
	 */
	{
		double alpha;
		double ww[vlp->q];
		
		upper_image.ideal=0; // mark y^* (dual vertex to be added) as a point

		boundlist* rows = boundlist_calloc(sol->p, 'u'); // to update the upper row bounds ub = Z' y in (P_2eta(y))
		boundlist_init_idx(rows, vlp->m + vlp->q + 1); // upper bounds for indices m+1 to m+p
		for(;;)
		{
			if (poly__get_vrtx(&upper_image))
				break;
			if (upper_image.ideal){
				ST_BT(upper_image.primal.sltn,upper_image.idx);
				continue;
			}
			for (size_t j=0;j<sol->p;j++)
			{
				rows->ub[j]=0;
				for (size_t k=0;k<vlp->q;k++)
					rows->ub[j] += sol->Z[k*sol->p+j] * upper_image.val[k];
			}
			if (opt->message_level >= 3) printf("process primal vertex - solve lp\n");
			lp_set_rows(0, rows);
			assert(lp_solve(0) == LP_OPTIMAL);
			lp_dual_solution_rows(0, ww, vlp->m + 1, vlp->q, 1);
			lp_dual_solution_rows(0, &alpha, vlp->m + vlp->q + sol->p + 1, 1, 1);
			// Z.v = w + alpha.eta
			for(size_t k=0; k<upper_image.dim-1; k++)
				upper_image.val[k]=ww[k]+alpha*sol->eta[k];
			upper_image.val[upper_image.dim-1]=alpha;
			if (lp_obj_val(0) > opt->eps_benson_phase1)
			{
				if (opt->message_level >= 3) printf("add dual vertex\n");
				poly__add_vrtx (&upper_image);
			}
			else
			{
				ST_BT (upper_image.primal.sltn,upper_image.idx);
			}
		}
		boundlist_free(rows);
	}
	/* PHASE 1 -- PRIMAL -- PART 3
	 *
	 *	compute matrix R and store it in sol->R:
	 *		select vertices y^* of lower image of dual vlp with last component zero
	 *		update last component such that c' y^* == 1
	 *		store result in matrix R
	 *
	 *	compute matrix H and store it in sol->H
	 *		compute dual cone of cone(R) and store directions as columns of H
	 */
	{
		double arr[vlp->q * upper_image.dual.cnt];
		size_t pp=0;
		for (size_t l=0; l<upper_image.dual.cnt; l++)
		{
			if (IS_ELEM(upper_image.dual.used,l)&&!IS_ELEM(upper_image.dual.ideal,l)&&ABS((VRTX_VAL(&upper_image.dual,l))[upper_image.dim-1])<opt->eps_phase1)
			{ 
				arr[pp*vlp->q+vlp->q-1]=1.0;
				for (size_t j=0; j<vlp->q-1; j++)
				{
					arr[pp*vlp->q+j]=upper_image.dual.data[l*vlp->q+j];
					arr[pp*vlp->q+vlp->q-1]-=sol->c[j]*arr[pp*vlp->q+j];
				}
				pp++;
			}
		}
		double arr_trans[vlp->q * pp];
		for (size_t l = 0; l < pp; l++)
			for (size_t k = 0; k < vlp->q; k++)
				arr_trans[k*pp+l]=arr[l*vlp->q+k];
		cone_vertenum(&sol->R, &sol->r,&sol->H, &sol->h,arr_trans,pp,vlp->q,opt,CONE_OUT_OFF,NO_SWAP);
	}
	if (POLY_TEST)
		poly__polyck (&upper_image);
	poly__kill (&upper_image);
} // end of phase1_primal
	
/*
 *
 *
 * PHASE 2 -- PRIMAL
 * 
 *
 *
 */
void phase2_init(soltype *sol, const vlptype *vlp)	// initialization for bounded case
{
	sol->R = calloc(sol->q * sol->p, sizeof(double));
	for (size_t i = 0; i < sol->p; i++)
		for (size_t j = 0; j < sol->q; j++)
			sol->R[j * sol->p + i] = sol->Z[j * sol->p + i];
	sol->r = sol->p;	
	
	sol->H = calloc(sol->q * sol->o, sizeof(double));
	for (size_t i = 0; i < sol->o; i++)
		for (size_t j = 0; j < sol->q; j++)
			sol->H[j * sol->o + i] = sol->Y[j * sol->o + i];
	sol->h = sol->o;
}

void phase2_primal(soltype *const sol, const vlptype *vlp, const opttype *opt)
{
	poly_args upper_image;

	fnc_prmtr.dim=vlp->q;
	fnc_prmtr.ip=sol->c;

	if (opt->solution == PRE_IMG_ON)
		init_image (&upper_image, vlp->q, vlp->n, vlp->m+vlp->q, &lowerV2upperH);
	else
	{
		assert(opt->solution==PRE_IMG_OFF);
		init_image (&upper_image, vlp->q, 0, 0, &lowerV2upperH);
	}

	init_P2(sol, vlp, INHOMOGENEOUS);
	lp_set_options(&(opt->lp), PHASE2_PRIMAL);

	{
		upper_image.ideal=0;
	
		boundlist *rows = boundlist_calloc(sol->r + 1, 'f');
		boundlist_init_idx(rows, vlp->m + vlp->q + 1);
		for (lp_idx j = 0; j < sol->r; j++)
		{
			rows->type[j] = 'u';
			lp_set_rows(0, rows);
			rows->type[j] = 'f';
			
			for(size_t k=0; k<upper_image.dim; k++)
				upper_image.val[k]=sol->R[k*sol->r+j];
			
			if (opt->message_level >= 3) printf("initialization - solve lp\n");
			lp_status_type lp_status = lp_solve(0);
			if (lp_status != LP_OPTIMAL)
			{
				if (lp_status == LP_INFEASIBLE)
					sol->status = VLP_INFEASIBLE;
				else
					sol->status = VLP_UNBOUNDED;
				break;
			}
			if (opt->solution == PRE_IMG_ON)
			{
				lp_dual_solution_rows (0,upper_image.val_primg_dl,1,vlp->m,(vlp->optdir==1?1:-1)); // u
				for (size_t k=0;k<vlp->q;k++) 
					upper_image.val_primg_dl[vlp->m+k]=(sol->c_dir==C_DIR_POS?1:-1)*upper_image.val[k]; // w
			}
			upper_image.val[vlp->q-1]=lp_obj_val(0);
			poly__add_vrtx (&upper_image);
		}
		boundlist_free(rows);
		poly__intl_apprx (&upper_image);
		if (sol->status == VLP_INFEASIBLE || sol->status == VLP_UNBOUNDED)
		{
			if (POLY_TEST)
				poly__polyck (&upper_image);
			poly__kill (&upper_image);
			return;
		}
	}
	
	/* PHASE 2 -- PRIMAL -- PART 2
	 *
	 *	main loop
	 *
	 */
	{
		double ww[vlp->q];
		boundlist* rows = boundlist_calloc(sol->r, 'u');	// list to update the upper row bounds Z'y in P_2(y)
		boundlist_init_idx(rows, vlp->m + vlp->q + 1);		// upper bounds for indices m+q+1 to m+q+r
		double yy[sol->q];									// vector to store yy = P x
		for(;;)
		{
			if (poly__get_vrtx (&upper_image))
			{
				break;
			}
			if (upper_image.ideal)  // rule out directions
			{
				ST_BT(upper_image.primal.sltn,upper_image.idx);
				continue;
			}
			for (size_t j=0;j<sol->r;j++)
			{
				rows->ub[j]=0;
				for (size_t k=0;k<vlp->q;k++)
					rows->ub[j]+=sol->R[k*sol->r+j]*upper_image.val[k];
			}
			if (opt->message_level >= 3) printf("process primal vertex - solve lp\n");
			lp_set_rows(0, rows);
			assert(lp_solve(0) == LP_OPTIMAL);
			lp_dual_solution_rows(0, ww, vlp->m + 1, vlp->q, 1);
			for (size_t k=0; k<vlp->q-1; k++)
			{
				upper_image.val[k]=ww[k];
			}	
			lp_primal_solution_cols(0, yy, vlp->n + 1, vlp->q, 1);
			double z = lp_obj_val(0);
			upper_image.val[upper_image.dim-1]=0;
			upper_image.ideal=0;
			for (size_t i = 0; i < vlp->q; i++) // y_star_q = w'* P * x
			{
				upper_image.val[upper_image.dim-1]+=yy[i]*ww[i];
			}
			if (z > opt->eps_benson_phase2)
			{
				if (opt->message_level >= 3) printf("add dual vertex\n");
				if (opt->solution == PRE_IMG_ON)
				{
					lp_dual_solution_rows (0, upper_image.val_primg_dl,1,vlp->m, (vlp->optdir==1?1:-1)); // u
					for (size_t k=0;k<vlp->q;k++)
						upper_image.val_primg_dl[vlp->m+k]=(sol->c_dir==C_DIR_POS?1:-1)*ww[k]; // w
				}
				poly__add_vrtx (&upper_image);
			}
			else
			{
				ST_BT(upper_image.primal.sltn,upper_image.idx);
				if (opt->solution == PRE_IMG_ON)
					lp_primal_solution_cols (0, upper_image.primal.data_primg+upper_image.primal.dim_primg*upper_image.idx,1,vlp->n, 1);
			}
		}
		boundlist_free(rows);
	}

	// write pre-image vectors for directions of upper image
	if (opt->solution == PRE_IMG_ON)
	{
		init_P2(sol, vlp, HOMOGENEOUS);
		// disable last constraint
		{
			boundlist* rows = boundlist_calloc(1, 'f');
			boundlist_init_idx(rows, vlp->m + vlp->q + sol->p + 1);
			lp_set_rows(0, rows);
			boundlist_free(rows);
		}
		boundlist* rows = boundlist_calloc(sol->p, 'u');	// list to update the upper row bounds Z'y in P_2(y)
		boundlist_init_idx(rows, vlp->m + vlp->q + 1);		// upper bounds for indices m+q+1 to m+q+r

		for(size_t i=0; i<upper_image.primal.cnt; i++)
		{
			if(IS_ELEM(upper_image.primal.used,i)&&IS_ELEM(upper_image.primal.ideal,i))
			{
				for (size_t j=0;j<sol->p;j++)
				{
					rows->ub[j]=0;
					for (size_t k=0; k<vlp->q; k++)
						rows->ub[j]+=sol->Z[k*sol->p+j]*upper_image.primal.data[i*vlp->q+k];
				}
				lp_set_rows(0, rows);
				assert(lp_solve(0) == LP_OPTIMAL);
				lp_primal_solution_cols (0, upper_image.primal.data_primg+upper_image.primal.dim_primg*i,1,vlp->n, 1);
			}
		}
		boundlist_free(rows);
	}

	// write pre-image vector for direction of lower image
	if (opt->solution == PRE_IMG_ON)
	{
		for(size_t i=0; i<upper_image.dual.cnt; i++)
			if(IS_ELEM(upper_image.dual.used,i)&&IS_ELEM(upper_image.dual.ideal,i))
				for (size_t j=0;j<upper_image.dual.dim_primg;j++)
					upper_image.dual.data_primg[upper_image.dual.dim_primg*i + j] = 0;
	}
	
	// save dual description of image for dual plot
	poly_args lower_image;
	if (opt->plot)
	{
		upper_image.primalV2dualH = &upperV2lowerH;
		init_image (&lower_image,upper_image.dim,0,0,upper_image.primalV2dualH);
		poly__swap (&upper_image,&lower_image);
		poly_chop(&lower_image);
		poly_normalize_dir(&lower_image);
	}
	
	// transform output for maximization problems and the case of negative c_q
	poly_trans_primal(vlp, sol, opt, &upper_image);
	
	// end of computations - stop timer
	gettimeofday(&t_end, NULL);
	
	poly_chop(&upper_image);
	poly_normalize_dir(&upper_image);
	poly__update_adjacence (&upper_image.dual);
	poly_output(&upper_image,opt,1,opt->solution,NO_SWAP,vlp->optdir==1?MIN_P_STR:MAX_P_STR, vlp->optdir==1?MIN_D_STR:MAX_D_STR,SOL_ENDING_STR);
	poly_count(&upper_image,sol,NO_SWAP);
	
	// write upper image to OFF file
	if (opt->plot)
	{
		// transform back for plotting
		poly_trans_primal(vlp, sol, opt, &upper_image);
		poly_plot_primal(vlp, sol, opt, &upper_image);
		poly_plot_dual(vlp, sol, opt, &lower_image);
		poly__kill (&lower_image);
	}

	if (POLY_TEST)
		poly__polyck (&upper_image);
	poly__kill (&upper_image);
}

/*
 *
 *
 *
 * D U A L
 *
 *
 *
 */

/* INIT_P1
 *
 *	initialize the lp for the problem
 *		(P_1(w))  min w^T y 	s.t.
 *		rls <= A x <= lhs, lb <= x <= ub (vlp constraints)
 *		P x - y == 0
 *		eta'P x + 0 y <= 1 (if type == HOMOGENEOUS) , infty (if type == INHOMOGENEOUS)
 *	where
 *		w = 0
 *		type = HOMOGENEOUS / INHOMOGENEOUS (vlp costraints of hom./inhom. problem)
 *		assuming that matrix A is set explicitly
 */

static void init_P1(soltype *const sol, const vlptype *vlp, lp_hom_type type)
{
	lp_update_extra_coeffs(1, 0); // update dimension
	{	// set lp constraint matrix
		list1d *rowlist = list1d_calloc(vlp->n + vlp->q);
		list1d_init_idx(rowlist, 1);
		for (lp_idx i = 0; i < vlp->q; i++)
			rowlist->data[vlp->n + i] = sol->eta[i];
		lp_set_mat_row(0, rowlist, vlp->m + vlp->q + 1);
		list1d_free(rowlist);
	}
	{	// set lp objective coefficients
		lp_clear_obj_coeffs(0);
	}
	// set row bounds of lp based on vlp constriants of (in)homogeneous problem
	if (type == HOMOGENEOUS)
		lp_set_rows_hom(0, vlp->rows);
	else
	{
		assert(type == INHOMOGENEOUS);
		lp_set_rows(0, vlp->rows);
	}
	{	// set remaining row bounds
		boundlist *rows = boundlist_calloc( vlp->q + 1, 's');
		boundlist_init_idx(rows, vlp->m + 1);
		if (type == HOMOGENEOUS)
		{
			rows->type[vlp->q] = 'u';
			rows->ub[vlp->q]=1.0;
		}
		else
		{
			assert(type == INHOMOGENEOUS);
			rows->type[vlp->q] = 'f';
		}
		lp_set_rows(0, rows);
		boundlist_free(rows);
	}
	// set column bounds of lp based on vlp constriants of (in)homogeneous problem
	if (type == HOMOGENEOUS)
		lp_set_cols_hom(0, vlp->cols);
	else
	{
		assert(type == INHOMOGENEOUS);
		lp_set_cols(0, vlp->cols);
	}	
	{	// set y_1,...,y_q to be free variables
		boundlist *cols = boundlist_calloc(vlp->q, 'f');
		boundlist_init_idx(cols, vlp->n + 1);
		lp_set_cols(0, cols);
		boundlist_free(cols);
	}
}

/*
 *
 *
 * PHASE 1 -- DUAL
 * 
 *
 *
 */
void phase1_dual(soltype *const sol, const vlptype *vlp, const opttype *opt)
{
	poly_args lower_image;
	
	fnc_prmtr.dim=vlp->q;
	fnc_prmtr.ip=sol->c;
		
	init_image (&lower_image, vlp->q, 0, 0, &upperV2lowerH);
	
	init_P1(sol, vlp, HOMOGENEOUS);	// assume that appropriate eta is stored in sol->eta
	lp_set_options(&(opt->lp), PHASE1_DUAL);

	/* PHASE 1 -- DUAL -- PART 1
	 *
	 * determine first outer approximation of upper image
	 *
	 */
	{
		double w[vlp->q];
		for(lp_idx i = 0; i < vlp->q; i++) // w = mean of columns of Z
		{
			w[i] = 0;
			for (lp_idx j = 0; j < sol->p; j++)
				w[i] += sol->Z[i * sol->p + j];
			w[i] /= sol->p;
		}
		lower_image.ideal=0;
		list1d *obj = list1d_calloc(vlp->q);
		list1d_init_idx(obj, vlp->n + 1);
		for(lp_idx i = 0; i < vlp->q; i++)
			obj->data[i] = w[i];
		lp_set_obj_coeffs(0, obj);
		if (opt->message_level >= 3) printf("initialization - solve lp\n");
		assert(lp_solve(0) == LP_OPTIMAL);
		lp_primal_solution_cols(0,lower_image.val,vlp->n+1,vlp->q, 1);
		poly__add_vrtx(&lower_image);
		lower_image.ideal=1;
		for(lp_idx j = 0; j < sol->o; j++)
		{
			for (size_t i=0; i<vlp->q; i++)
				lower_image.val[i]=sol->Y[i * sol->o + j];
			poly__add_vrtx(&lower_image);
		}
		list1d_free(obj);
		poly__intl_apprx (&lower_image);
	}

	/* PHASE 1 -- DUAL -- PART 2
	 *
	 * main loop
	 * 
	 */
	{
		list1d *obj = list1d_calloc(vlp->q);
		list1d_init_idx(obj, vlp->n + 1);

		for(;;)
		{
			if (poly__get_vrtx(&lower_image))
				break;
			if (lower_image.ideal){
				ST_BT(lower_image.primal.sltn,lower_image.idx);
				continue;
			}
			// compute w(y^*)
			obj->data[vlp->q - 1] = 1;
			for (size_t i = 0; i < vlp->q - 1; i++)
			{
				obj->data[i] = lower_image.val[i];
				obj->data[vlp->q - 1] -= lower_image.val[i] * sol->c[i];
			}
			if (opt->message_level >= 3) printf("process dual vertex - solve lp\n");
			lp_set_obj_coeffs(0, obj);
			assert(lp_solve(0) == LP_OPTIMAL);
			double opt_val=lower_image.val[lower_image.dim-1]; // store value before lower_image.val is used for input
			lp_primal_solution_cols(0,lower_image.val,vlp->n+1,vlp->q, 1);
			if (opt_val-lp_obj_val(0)>opt->eps_benson_phase1)
			{
				if (opt->message_level >= 3)
					printf("add primal vertex\n");
				poly__add_vrtx(&lower_image);
			}
			else
				ST_BT (lower_image.primal.sltn,lower_image.idx);
		}
		list1d_free(obj);
	}	

	/* PHASE 1 -- DUAL -- PART 3
	 *
	 *	compute matrix R and store it in sol->R:
	 *		select vertices y^* of lower image of dual vlp with last component zero
	 *		update last component such that c' y^* == 1
	 *		store result in matrix R
	 *
	 *	compute matrix H and store it in sol->H
	 *		compute dual cone of cone(R) and store directions as columns of H
	 */
	{
		double arr[vlp->q * lower_image.primal.cnt];
		size_t pp=0;
		for (size_t l=0; l<lower_image.primal.cnt; l++)
		{
			if (IS_ELEM(lower_image.primal.used,l)&&!IS_ELEM(lower_image.primal.ideal,l)&&ABS((VRTX_VAL(&lower_image.primal,l))[vlp->q-1])<opt->eps_phase1)
			{ 
				arr[pp*vlp->q+vlp->q-1]=1.0;
				for (size_t j=0; j<vlp->q-1; j++)
				{
					arr[pp*vlp->q+j]=lower_image.primal.data[l*vlp->q+j];
					arr[pp*vlp->q+vlp->q-1]-=sol->c[j]*arr[pp*vlp->q+j];
				}
				pp++;
			}
		}
		double arr_trans[vlp->q * pp];
		for (size_t l = 0; l < pp; l++)
			for (size_t k = 0; k < vlp->q; k++)
				arr_trans[k*pp+l]=arr[l*vlp->q+k];
		cone_vertenum(&sol->R,&sol->r,&sol->H,&sol->h,arr_trans,pp,vlp->q,opt,CONE_OUT_OFF,NO_SWAP);
	}	
	if (POLY_TEST)
		poly__polyck (&lower_image);
	poly__kill (&lower_image);
} // end of phase1_dual

/*
 *
 *
 * PHASE 2 -- DUAL
 *
 *
 */

void phase2_dual(soltype *const sol, const vlptype *vlp, const opttype *opt)
{
	poly_args lower_image;

	fnc_prmtr.dim=vlp->q;
	fnc_prmtr.ip=sol->c;

	if (opt->solution == PRE_IMG_ON)
		init_image (&lower_image, vlp->q, vlp->m+vlp->q, vlp->n, &upperV2lowerH);
	else
		init_image (&lower_image, vlp->q, 0, 0, &upperV2lowerH);

	init_P1(sol, vlp, INHOMOGENEOUS); // assume that appropriate eta is stored in sol->eta
	lp_set_options(&(opt->lp), PHASE2_DUAL);

	/* PHASE 2 -- DUAL -- PART 1
	 *
	 * determine first outer approximation of upper image
	 *
	 */
	{
		double w[vlp->q];
		for(lp_idx i = 0; i < vlp->q; i++) // w = mean of columns of R
		{
			w[i] = 0;
			for (lp_idx j = 0; j < sol->r; j++)
				w[i] += sol->R[i * sol->r + j];
			w[i] /= sol->r;
		}
		lower_image.ideal=0;
		list1d *obj = list1d_calloc(vlp->q);
		list1d_init_idx(obj, vlp->n + 1);
		for(lp_idx i = 0; i < vlp->q; i++)
			obj->data[i] = w[i];
		lp_set_obj_coeffs(0, obj);
		if (opt->message_level >= 3) printf("initialization - solve lp\n");
		lp_status_type lp_status = lp_solve(0);
		if (lp_status != LP_OPTIMAL)
		{
			list1d_free(obj);
			if (POLY_TEST)
				poly__polyck (&lower_image);
			poly__kill (&lower_image);
			if (lp_status == LP_INFEASIBLE)
				sol->status = VLP_INFEASIBLE;
			else
				sol->status = VLP_UNBOUNDED;
			return;
		}
		lp_primal_solution_cols(0,lower_image.val,vlp->n+1,vlp->q, 1);
		if (opt->solution == PRE_IMG_ON)
			lp_primal_solution_cols (0, lower_image.val_primg_dl,1,vlp->n, 1);
		poly__add_vrtx(&lower_image);
		/*
		 * add columns of H (generating vectors recession cone of upper image computed in phase 1) as directions of the upper image
		 */
		lower_image.ideal=1;
		for(lp_idx j = 0; j < sol->h; j++)
		{
			for (size_t i=0; i<vlp->q; i++)
				lower_image.val[i] = sol->H[i * sol->h + j]; 
			poly__add_vrtx(&lower_image);
		}
		list1d_free(obj);
		poly__intl_apprx (&lower_image);
	}

	/* PHASE 2 -- DUAL -- PART 2
	 *
	 * main loop
	 *
	 */
	{
		lower_image.ideal=0;
		list1d *obj = list1d_calloc(vlp->q);
		list1d_init_idx(obj, vlp->n + 1);
		int lp_status;
		for(;;)
		{
			if (poly__get_vrtx(&lower_image))
				break;
			if (lower_image.ideal){
				ST_BT(lower_image.primal.sltn,lower_image.idx);
				continue;
			}
			// compute w(y^*)
			obj->data[vlp->q - 1] = 1;
			for (lp_idx i = 0; i < vlp->q - 1; i++)
			{
				obj->data[i] = lower_image.val[i];
				obj->data[vlp->q - 1] -= lower_image.val[i] * sol->c[i];
			}
			if (opt->message_level >= 3) printf("process dual vertex - solve lp\n");
			lp_set_obj_coeffs(0, obj);
			lp_status = lp_solve(0);
			if (lp_status != LP_OPTIMAL)
			{
				assert(lp_status == LP_UNBOUNDED);
				sol->status = VLP_UNBOUNDED; 
				break;
			}
			double opt_val=lower_image.val[lower_image.dim-1]; // store this value before lower_image.val is used for input
			lp_primal_solution_cols(0,lower_image.val,vlp->n+1,vlp->q, 1);
			if (opt_val-lp_obj_val(0)>opt->eps_benson_phase2)
			{
				if (opt->message_level >= 3)
					printf("add primal vertex\n");
				if (opt->solution == PRE_IMG_ON)
					lp_primal_solution_cols (0, lower_image.val_primg_dl,1,vlp->n, 1); //x
				poly__add_vrtx(&lower_image);
			}
			else
			{
				ST_BT (lower_image.primal.sltn,lower_image.idx);
				if (opt->solution == PRE_IMG_ON)
				{
					lp_dual_solution_cols (0, lower_image.primal.data_primg+lower_image.primal.dim_primg*lower_image.idx,1,vlp->m, (vlp->optdir==1?1:-1)); //u
					for (size_t k=0; k<vlp->q; k++)
						lower_image.primal.data_primg[lower_image.primal.dim_primg*lower_image.idx + vlp->m + k]=(sol->c_dir==C_DIR_POS?1:-1)*obj->data[k]; //w
				}
			}
		}
		list1d_free(obj);
	}
	if (sol->status == VLP_UNBOUNDED)
	{
		if (POLY_TEST)
			poly__polyck (&lower_image);
		poly__kill (&lower_image);
		return;
	}
	
	// write pre-image vectors for directions of upper image
	if (opt->solution == PRE_IMG_ON)
	{
		init_P2(sol, vlp, HOMOGENEOUS);
		// disable last constraint
		{
			boundlist* rows = boundlist_calloc(1, 'f');
			boundlist_init_idx(rows, vlp->m + vlp->q + sol->p + 1);
			lp_set_rows(0, rows);
			boundlist_free(rows);
		}
		boundlist* rows = boundlist_calloc(sol->p, 'u');	// list to update the upper row bounds Z'y in P_2(y)
		boundlist_init_idx(rows, vlp->m + vlp->q + 1);		// upper bounds for indices m+q+1 to m+q+r
		
		for(size_t i=0; i<lower_image.dual.cnt; i++)
		{
			if(IS_ELEM(lower_image.dual.used,i)&&IS_ELEM(lower_image.dual.ideal,i))
			{
				for (size_t j=0;j<sol->p;j++)
				{
					rows->ub[j]=0;
					for (size_t k=0; k<vlp->q; k++)
						rows->ub[j]+=sol->Z[k*sol->r+j]*lower_image.dual.data[i*vlp->q+k];
				}
				lp_set_rows(0, rows);
				assert(lp_solve(0) == LP_OPTIMAL);
				lp_primal_solution_cols (0, lower_image.dual.data_primg+lower_image.dual.dim_primg*i,1,vlp->n, 1);
			}
		}
		boundlist_free(rows);
	}
	
	// write pre-image vector for direction of lower image
	if (opt->solution == PRE_IMG_ON)
	{
		for(size_t i=0; i<lower_image.primal.cnt; i++)
			if(IS_ELEM(lower_image.primal.used,i)&&IS_ELEM(lower_image.primal.ideal,i))
				for (size_t j=0;j<lower_image.primal.dim_primg;j++)
					lower_image.primal.data_primg[lower_image.primal.dim_primg*i + j] = 0;
	}
	
	// save dual description of image for dual plot
	poly_args upper_image;
	if (opt->plot)
	{
		lower_image.primalV2dualH = &lowerV2upperH;
		init_image (&upper_image,vlp->q,0,0,lower_image.primalV2dualH);
		poly__swap (&lower_image,&upper_image);
		poly_chop(&upper_image);
		poly_normalize_dir(&upper_image);
	}
	
	// transform output for maximization problems and the case of negative c_q
	poly_trans_dual(vlp,sol,opt,&lower_image);
	poly_chop(&lower_image);
	poly_normalize_dir(&lower_image);
	poly__update_adjacence (&lower_image.dual);
	
	// end of computations - stop timer
	gettimeofday(&t_end, NULL);
	
	poly_output(&lower_image,opt,1,opt->solution,SWAP,vlp->optdir==1?MIN_P_STR:MAX_P_STR, vlp->optdir==1?MIN_D_STR:MAX_D_STR,SOL_ENDING_STR);
	poly_count(&lower_image,sol,SWAP);
	
	if (POLY_TEST)
		poly__polyck (&lower_image);

	// graphics output
	if (opt->plot)
	{
		// transform back for plot
		poly_trans_dual(vlp,sol,opt,&lower_image);
		poly_plot_dual(vlp,sol,opt,&lower_image);
		poly_plot_primal(vlp,sol,opt,&upper_image);
		poly__kill (&upper_image);
	}

	poly__kill (&lower_image);

} // end of phase2_dual
