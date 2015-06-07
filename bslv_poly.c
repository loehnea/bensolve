/*
This file is part of BENSOLVE - VLP solver

Copyright © 2014-2015 Andreas Löhne and Benjamin Weißing

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (see the reference manual). If not,
see <http://www.gnu.org/licenses/>
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "bslv_poly.h"


size_t fnc_dim;

static void cone_polar (double *dual_point, int is_dir, double *primal_hyperplane)
{
	double *out;
	*(primal_hyperplane+fnc_dim)=is_dir?0:-1.0;
	for (out=primal_hyperplane;out<primal_hyperplane+fnc_dim;out++,dual_point++)
		*out=*dual_point;


	return;
}

void poly__set_default_args (poly_args *args, size_t dim)
{
	fnc_dim = dim;
	args->dim=dim;
	args->eps = 1e-08;
	args->dim_primg_prml=0;
	args->dim_primg_dl=0;
	args->primalV2dualH = NULL;
	args->dualV2primalH = &cone_polar;

	
	return;
}

void poly__initialise (poly_args *args)
{
	size_t k;

	args->primal.dim=args->dim;
	args->primal.dim_primg=args->dim_primg_prml;

	args->dual.dim=args->dim;
	args->dual.dim_primg=args->dim_primg_dl;

	args->dual.dim=args->dim;
	args->primal.blcks=1;
	args->dual.blcks=1;
	poly__poly_init (&args->primal);
	poly__poly_init (&args->dual);

	args->primal.dual=&args->dual;
	args->dual.dual=&args->primal;


	args->val = (double *)malloc(args->dim*sizeof(double));
	args->val_primg_prml=(double *)malloc(args->dim_primg_prml*sizeof(double));
	args->val_primg_dl=(double *)malloc(args->dim_primg_dl*sizeof(double));

	args->dual.v2h = args->dualV2primalH;
	args->primal.v2h=args->primalV2dualH;


	add_vrtx(&args->dual);
	k=args->dim-1;
	*(args->dual.data+k)=-1.0;
	while (k--)
		*(args->dual.data+k)=0;
	k=args->dim_primg_dl;
	while (k--)
		*(args->dual.data_primg+k)=0;

	ST_BT(args->dual.ideal,0);
	args->init_data.H = (double *)malloc(args->dim*args->dim*sizeof(double));
	args->init_data.R = (double *)malloc(args->dim*(args->dim+1)/2*sizeof(double));
	args->init_data.alph = (double *)malloc(args->dim*sizeof(double));
	poly__list_init (&args->init_data.queue);
	poly__list_init (&args->init_data.gnrtrs);
	args->init_data.intlsd=0;


	return;
}

int poly__add_vrtx (poly_args *args)
{
	double *val,sum,*hp = (double *)malloc((args->dual.dim+1)*sizeof(double));
	size_t j;

	add_vrtx (&args->dual);

	if(args->ideal)
		ST_BT(args->dual.ideal,(args->dual.cnt-1));
	for (val=args->val;val<args->val+args->dim;val++)
		*(args->dual.data+args->dual.dim*(args->dual.cnt-1)+(val-args->val))=*val;
	for (val=args->val_primg_dl;val<args->val_primg_dl+args->dim_primg_dl;val++)
		*(args->dual.data_primg+args->dual.dim_primg*(args->dual.cnt-1)+(val-args->val_primg_dl))=*val;


	args->dualV2primalH (args->val,args->ideal,hp);
	if (args->init_data.intlsd){
		for (args->idx=0;args->idx<args->primal.cnt;args->idx++)
			if (IS_ELEM(args->primal.used,args->idx)){
				sum=0;
				for (j=0;j<args->dim;j++)
					sum+=*(hp+j)**(VRTX_VAL(&args->primal,args->idx)+j);
				if (sum<(IS_ELEM(args->primal.ideal,args->idx)?0:*(hp+args->dim))-POLY_EPS)
					break;
			}

		if (args->idx<args->primal.cnt)
			poly__cut (&args->primal, args->idx, hp);
		else{
			UNST_BT (args->dual.used,(args->dual.cnt-1));
			free (hp);
			return EXIT_FAILURE;
		}

		for (size_t *k = (args->dual.incidence+args->dual.cnt-1)->data;k<(args->dual.incidence+args->dual.cnt-1)->data+(args->dual.incidence+args->dual.cnt-1)->cnt;k++)
			for (size_t *j = (args->dual.incidence+args->dual.cnt-1)->data;j<k;j++)
				if (edge_test(args->dual.dual,*k,*j)){
					add_lst_elem(args->dual.dual->adjacence+*k,*j);
					add_lst_elem(args->dual.dual->adjacence+*j,*k);
				}
	}else
		add_lst_elem (&args->init_data.queue,args->dual.cnt-1);

	free (hp);


	return EXIT_SUCCESS;
}

int poly__intl_apprx (poly_args *poly)
{
	double max_nrm, *hp;
	size_t max_idx,*perm;

	if (poly->init_data.queue.cnt<poly->dim)
		return EXIT_FAILURE;

	hp = (double *)malloc((poly->dim+1)*poly->init_data.queue.cnt*sizeof(double));
	for (size_t *dual_vrtx=poly->init_data.queue.data;dual_vrtx<poly->init_data.queue.data+poly->init_data.queue.cnt;dual_vrtx++)
		poly->dualV2primalH (poly->dual.data+*dual_vrtx*poly->dim,IS_ELEM(poly->dual.ideal,*dual_vrtx),hp+(dual_vrtx-poly->init_data.queue.data)*(poly->dim+1));

	perm = (size_t *)malloc((poly->dim+1)*sizeof(size_t));
	*perm=0;
	while (poly->init_data.gnrtrs.cnt<poly->dim){
		max_nrm=0;
		for (size_t k=0;k<poly->init_data.queue.cnt;k++)
			if (max_nrm<bslv__normalise(hp+k*(poly->dim+1),poly->init_data.H,poly->init_data.R,poly->init_data.gnrtrs.cnt,poly->dim)){
				max_nrm=bslv__normalise(hp+k*(poly->dim+1),poly->init_data.H,poly->init_data.R,poly->init_data.gnrtrs.cnt,poly->dim);
				max_idx=k;
			}
		if (max_nrm<1.0e-10){
			free (hp);
			return EXIT_FAILURE;
		}
		bslv__normalise (hp+max_idx*(poly->dim+1),poly->init_data.H,poly->init_data.R,poly->init_data.gnrtrs.cnt,poly->dim);
		*(poly->init_data.alph+poly->init_data.gnrtrs.cnt)=*(hp+max_idx*(poly->dim+1)+poly->dim);
		add_lst_elem (&poly->init_data.gnrtrs,*(poly->init_data.queue.data+max_idx));
		*(perm+poly->init_data.gnrtrs.cnt)=*(poly->init_data.queue.data+max_idx);
		for (double *hp_val=hp+max_idx*(poly->dim+1);hp_val<hp+(poly->dim+1)*(max_idx+1);hp_val++)
			*hp_val=*(hp_val+(poly->dim+1)*(poly->init_data.queue.cnt-max_idx-1));
		*(poly->init_data.queue.data+max_idx)=*(poly->init_data.queue.data+--poly->init_data.queue.cnt);
	}
	poly__poly_initialise (&poly->primal,poly->init_data.H,poly->init_data.R,poly->init_data.alph,perm);
	free (perm);
	poly->init_data.intlsd=1;

	for (size_t *k_ptr=poly->init_data.queue.data;k_ptr<poly->init_data.queue.data+poly->init_data.queue.cnt;k_ptr++)
		UNST_BT(poly->dual.used,*k_ptr);
	for (size_t k=0;k<poly->init_data.queue.cnt;k++){
		for (size_t j=0;j<poly->dim;j++)
			*(poly->val+j)=*(poly->dual.data+*(poly->init_data.queue.data+k)*poly->dim+j);
		poly->ideal=IS_ELEM(poly->dual.ideal,(*(poly->init_data.queue.data+k)));
		poly__add_vrtx (poly);
	}

	free (hp);
	free (poly->init_data.queue.data);
	free (poly->init_data.gnrtrs.data);
	free (poly->init_data.H);
	free (poly->init_data.R);
	free (poly->init_data.alph);


	return EXIT_SUCCESS;
}

int poly__get_vrtx (poly_args *args)
{
	size_t k;

	for (args->idx=0; args->idx<args->primal.cnt; args->idx++)
		if (IS_ELEM(args->primal.used,args->idx)&&!IS_ELEM(args->primal.sltn,args->idx))
			break;
	if(args->idx==args->primal.cnt)
		return EXIT_FAILURE;

	args->ideal = IS_ELEM (args->primal.ideal, args->idx);
	for (k=0; k<args->primal.dim; k++)
		args->val[k]=*(VRTX_VAL(&args->primal,args->idx)+k);


	return EXIT_SUCCESS;
}

void poly__poly_init(polytope *poly)
{
	poly->cnt = 0;
	poly->data = (double *) malloc (poly->blcks*VRTXBLCK*poly->dim*sizeof(double));
	poly->data_primg = (double *) malloc (poly->blcks*VRTXBLCK*poly->dim_primg*sizeof(double));
	poly->ip = (double *)malloc(poly->dim*sizeof(double));
	poly->adjacence = (poly_list *) malloc (poly->blcks*VRTXBLCK*sizeof(poly_list));
	for (poly_list *lst=poly->adjacence; lst<poly->adjacence+poly->blcks*VRTXBLCK; lst++)
		poly__list_init(lst);
	poly->incidence = (poly_list *) malloc (poly->blcks*VRTXBLCK*sizeof(poly_list));
	for (poly_list *lst=poly->incidence; lst<poly->incidence+poly->blcks*VRTXBLCK; lst++)
		poly__list_init(lst);
	poly->used = (vrtx_strg *) calloc (poly->blcks*ALLOCFCTR,sizeof(vrtx_strg));
	poly->ideal = (vrtx_strg *) calloc (poly->blcks*ALLOCFCTR,sizeof(vrtx_strg));
	poly->sltn = (vrtx_strg *) calloc (poly->blcks*ALLOCFCTR,sizeof(vrtx_strg));


	return;
}

void poly__list_init (poly_list *list)
{
	list->cnt=0;
	list->blcks=1;
	list->data=(size_t *)malloc(LSTBLCK*sizeof(size_t));


	return;
}

void poly__kill (poly_args *args)
{
	poly__poly_kill (&args->primal);
	poly__poly_kill (&args->dual);
	free (args->val);
	if(!args->init_data.intlsd){
		free (args->init_data.H);
		free (args->init_data.R);
		free (args->init_data.alph);
		free (args->init_data.queue.data);
		free (args->init_data.gnrtrs.data);
	}
	free (args->val_primg_prml);
	free (args->val_primg_dl);


	return;
}

void poly__poly_kill (polytope *poly)
{
	free (poly->data);
	free (poly->data_primg);
	free (poly->ip);
	free (poly->used);
	free (poly->ideal);
	free (poly->sltn);
	for (size_t k=0; k<poly->blcks*VRTXBLCK; k++){
		free ((poly->adjacence+k)->data);
		free ((poly->incidence+k)->data);
	}
	free (poly->adjacence);
	free (poly->incidence);


	return;
}

void poly_defrag (polytope *poly)
{
	while (!IS_ELEM(poly->used,(poly->cnt-1)))
		poly->cnt--;
	for (size_t k=0; k<poly->cnt-1; k++){
		if(!IS_ELEM(poly->used,k)){
			ST_BT(poly->used,k);
			UNST_BT(poly->used,poly->cnt);
			vrtx_cpy (poly,--poly->cnt,k);
			while (!IS_ELEM(poly->used,(poly->cnt-1)))
				poly->cnt--;
		}
	}


	return;
}

void poly__initialise_permutation (polytope *poly, permutation *prm)
{
	size_t idx;

	prm->cnt=0;
	prm->data=(size_t *)malloc(poly->cnt*sizeof(size_t));
	prm->inv = (size_t *)malloc(poly->cnt*sizeof(size_t));
	for (idx=0;idx<poly->cnt;idx++)
		if (IS_ELEM(poly->used,idx)){
			*(prm->data+prm->cnt)=idx;
			*(prm->inv+idx)=prm->cnt++;
		}
	prm->data=(size_t *)realloc((void *)prm->data,prm->cnt*sizeof(size_t));


	return;
}

void poly__kill_permutation (permutation *prm)
{
	free (prm->data);
	free (prm->inv);


	return;
}

void poly__vrtx2file (polytope *poly, permutation *prm, const char *fname, const char *frmt)
{
	size_t *idx;
	double *val;

	FILE *stream = fname?fopen(fname,"w"):stdout;

	for (idx=prm->data;idx<prm->data+prm->cnt;idx++){
		fprintf (stream, "%-1.1d ", 1-(int)IS_ELEM(poly->ideal,*idx));
		for (val=poly->data+*idx*poly->dim;val<poly->data+(*idx+1)*poly->dim;val++)
			fprintf(stream, frmt?frmt:"%g ", *val);
		fseek (stream, -sizeof(char), SEEK_CUR);
		fprintf (stream, "\n");
	}
	if(fname)
		fclose(stream);


	return;
}

void poly__primg2file (polytope *poly, permutation *prm, const char *fname, const char *frmt)
{
	size_t *idx;
	double *val;
	FILE *stream = fname?fopen(fname,"w"):stdout;

	for (idx=prm->data;idx<prm->data+prm->cnt;idx++)
		if (IS_ELEM(poly->sltn,*idx)){
			for (val=poly->data_primg+*idx*poly->dim_primg;val<poly->data_primg+(*idx+1)*poly->dim_primg;val++)
				fprintf(stream, frmt?frmt:"%g ", *val);
			fseek (stream, -sizeof(char), SEEK_CUR);
			fprintf (stream, "\n");
		}
	if(fname)
		fclose(stream);


	return;
}

void poly__adj2file (polytope *poly, permutation *prm, const char *fname, const char *frmt)
{
	size_t *vrtx, *nghbr;
	FILE *stream = fname?fopen(fname,"w"):stdout;

	for (vrtx=prm->data;vrtx<prm->data+prm->cnt;vrtx++){
		for (nghbr=(poly->adjacence+*vrtx)->data;nghbr<(poly->adjacence+*vrtx)->data+(poly->adjacence+*vrtx)->cnt;nghbr++)
			fprintf (stream, frmt?frmt:"%u ", (unsigned int)(*(prm->inv+*nghbr)));
		fseek (stream, -sizeof(char), SEEK_CUR);
		fprintf (stream, "\n");
	}
	if(fname)
		fclose(stream);

	return;
}

void poly__inc2file (polytope *poly, permutation *prm, permutation *prm_dual, const char *fname, const char *frmt)
{
	size_t *fct,*vrtx;
	FILE *stream = fname?fopen(fname,"w"):stdout;

	for (fct=prm_dual->data;fct<prm_dual->data+prm_dual->cnt;fct++){
		for (vrtx=(poly->dual->incidence+*fct)->data;vrtx<(poly->dual->incidence+*fct)->data+(poly->dual->incidence+*fct)->cnt;vrtx++)
			fprintf (stream, frmt?frmt:"%u ",(unsigned int)(*(prm->inv+*vrtx)));
		fseek (stream, -sizeof(char), SEEK_CUR);
		fprintf(stream, "\n");
	}
	if(fname)
		fclose(stream);

	return;
}

void add_vrtx (polytope *poly)
{
	if (++poly->cnt==poly->blcks*VRTXBLCK){
		poly->data = (double *) realloc (poly->data,++poly->blcks*poly->dim*VRTXBLCK*sizeof(double));
		poly->data_primg = (double *) realloc (poly->data_primg,poly->blcks*poly->dim_primg*VRTXBLCK*sizeof(double));
		poly->adjacence = (poly_list *) realloc (poly->adjacence, poly->blcks*VRTXBLCK*sizeof(poly_list));
		for (poly_list *lst=poly->adjacence+(poly->blcks-1)*VRTXBLCK;lst<poly->adjacence+poly->blcks*VRTXBLCK;lst++){
			lst->cnt=0;
			lst->blcks=1;
			lst->data = (size_t *) malloc (LSTBLCK*sizeof(size_t));
		}
		poly->incidence = (poly_list *) realloc (poly->incidence, poly->blcks*VRTXBLCK*sizeof(poly_list));
		for (poly_list *lst=poly->incidence+(poly->blcks-1)*VRTXBLCK;lst<poly->incidence+poly->blcks*VRTXBLCK;lst++){
			lst->cnt=0;
			lst->blcks=1;
			lst->data = (size_t *) malloc (LSTBLCK*sizeof(size_t));
		}
		poly->used = (vrtx_strg *) realloc (poly->used, poly->blcks*ALLOCFCTR*sizeof(vrtx_strg));
		for (size_t k=0;k<ALLOCFCTR;k++)
			*(poly->used+(poly->blcks-1)*ALLOCFCTR+k) = 0;
		poly->ideal = (vrtx_strg *) realloc (poly->ideal, poly->blcks*ALLOCFCTR*sizeof(vrtx_strg));
		for (size_t k=0;k<ALLOCFCTR;k++)
			*(poly->ideal+(poly->blcks-1)*ALLOCFCTR+k) = 0;
		poly->sltn = (vrtx_strg *) realloc (poly->sltn, poly->blcks*ALLOCFCTR*sizeof(vrtx_strg));
		for (size_t k=0;k<ALLOCFCTR;k++)
			*(poly->sltn+(poly->blcks-1)*ALLOCFCTR+k) = 0;
	}
	ST_BT(poly->used,(poly->cnt-1));


	return;
}

void add_lst_elem (poly_list *lst, size_t elem)
{
	if (lst->cnt++==lst->blcks*LSTBLCK)
		lst->data = (size_t *) realloc (lst->data, ++lst->blcks*LSTBLCK*sizeof(size_t));
	*(lst->data+lst->cnt-1) = elem;


	return;
}

void delete_lst_elem (poly_list *lst, size_t idx)
{
	*(lst->data+idx)=*(lst->data+--lst->cnt);


	return;
}

int edge_test (polytope *poly, size_t v1, size_t v2)
{
	int ret;
	size_t *k,*j,*fct;
	poly_list mutual_facets={0,1,NULL},mutual_vertices={0,1,NULL};
	mutual_facets.data=(size_t *)malloc(LSTBLCK*sizeof(size_t));
	mutual_vertices.data=(size_t *)malloc(LSTBLCK*sizeof(size_t));

	for (k=(poly->incidence+v1)->data;k<(poly->incidence+v1)->data+(poly->incidence+v1)->cnt;k++)
		for (j=(poly->incidence+v2)->data;j<(poly->incidence+v2)->data+(poly->incidence+v2)->cnt;j++)
			if(*k==*j){
				add_lst_elem (&mutual_facets,*k);
				break;
			}

	if (poly->dim==1)
		ret=1;
	else if (mutual_facets.cnt<poly->dim-1)
		ret=0;
	else{
		for (k=(poly->dual->incidence+*mutual_facets.data)->data;k<(poly->dual->incidence+*mutual_facets.data)->data+(poly->dual->incidence+*mutual_facets.data)->cnt;k++)
			if (*k==v1||*k==v2)
				continue;
			else
				add_lst_elem (&mutual_vertices,*k);
		*mutual_facets.data=*(mutual_facets.data+--mutual_facets.cnt);
		for (fct=mutual_facets.data;fct<mutual_facets.data+mutual_facets.cnt;fct++)
			for (j=mutual_vertices.data;j<mutual_vertices.data+mutual_vertices.cnt;j++){
				for (k=(poly->incidence+*j)->data;k<(poly->incidence+*j)->data+(poly->incidence+*j)->cnt;k++)
					if (*fct==*k)
						break;
				if (k==(poly->incidence+*j)->data+(poly->incidence+*j)->cnt)
					*(j--)=*(mutual_vertices.data+--mutual_vertices.cnt);
			}

		if (mutual_vertices.cnt>0)
			ret=0;
		else
			ret=1;
	}


	free (mutual_facets.data);
	free (mutual_vertices.data);
	return ret;
}

void vrtx_cpy (polytope *poly, size_t idx, size_t cpy)
{
	size_t *k,*j;
	double *src,*dst;
	UNST_BT (poly->ideal,cpy);
	UNST_BT (poly->sltn,cpy);
	if (IS_ELEM(poly->sltn,idx)){
		ST_BT(poly->sltn,cpy);
		UNST_BT(poly->sltn,idx);
	}
	if (IS_ELEM(poly->ideal,idx)){
		ST_BT(poly->ideal,cpy);
		UNST_BT(poly->ideal,idx);
	}
	UNST_BT(poly->used,idx);
	(poly->adjacence+cpy)->cnt=0;
	(poly->incidence+cpy)->cnt=0;
	src=poly->data+poly->dim*idx;
	for (dst=poly->data+poly->dim*cpy;dst<poly->data+poly->dim*(cpy+1);dst++)
		*dst=*(src++);

	src=poly->data_primg+poly->dim_primg*idx;
	for (dst=poly->data_primg+poly->dim_primg*cpy;dst<poly->data_primg+poly->dim_primg*(cpy+1);dst++)
		*dst=*(src++);

	for (k=(poly->adjacence+idx) ->data; k< (poly->adjacence+idx) ->data+(poly->adjacence+idx) ->cnt; k++){
		add_lst_elem(poly->adjacence+cpy, *k);
		for (j=(poly->adjacence+*k) ->data; j< (poly->adjacence+*k) ->data+(poly->adjacence+*k) ->cnt; j++)
			if (*j==idx){
				*j=cpy;
				break;
			}
	}

	for (k=(poly->incidence+idx)->data; k<(poly->incidence+idx)->data+(poly->incidence+idx)->cnt; k++){
		add_lst_elem(poly->incidence+cpy,*k);
		for (j=(poly->dual->incidence+*k)->data; j<(poly->dual->incidence+*k)->data+(poly->dual->incidence+*k)->cnt; j++)
			if (*j==idx){
				*j=cpy;
				break;
			}
	}


	return;
}


void poly__cut (polytope *poly, size_t v, double *ct_hp)
{
	size_t *k,*j,*l,*js,*ks,*ls,v_out;
	double alpha=*(ct_hp+poly->dim),mu,*vrtx=(double *)malloc(poly->dim*sizeof(double)),*drctn=(double *)malloc(poly->dim*sizeof(double)),*val;

	int smpl=1;
	UNST_BT (poly->used, v);
	double tmp=0;
	for (size_t kk=0;kk<poly->dim;kk++)
		tmp+=ct_hp[kk]**(VRTX_VAL(poly,v)+kk);

	if (tmp>((IS_ELEM(poly->ideal,v))?0:alpha)-/* 1.0e-8 */POLY_EPS){
		smpl=0;
		v_out=poly->cnt;
		add_vrtx(poly);
		add_lst_elem (poly->dual->incidence+poly->dual->cnt-1,v_out);
		add_lst_elem (poly->incidence+v_out, poly->dual->cnt-1);
		for (val=poly->data+v_out*poly->dim;val<poly->data+(v_out+1)*poly->dim;val++)
			*val=*(val-poly->dim*(v_out-v));
		if(IS_ELEM(poly->ideal,v))
			ST_BT(poly->ideal,v_out);
		if(IS_ELEM(poly->sltn,v)){
			ST_BT(poly->sltn,v_out);
			for (val=poly->data_primg+v_out*poly->dim_primg;val<poly->data_primg+(v_out+1)*poly->dim_primg;val++)
				*val=*(val-poly->dim_primg*(v_out-v));
		}
	}

	for (k=(poly->adjacence+v)->data;k<(poly->adjacence+v)->data+(poly->adjacence+v)->cnt;k++){
		if (!IS_ELEM(poly->used,*k))
			continue;
		tmp=0;
		for (size_t kk=0;kk<poly->dim;kk++)
			tmp+=ct_hp[kk]**(VRTX_VAL(poly,*k)+kk);
		if (tmp>(IS_ELEM(poly->ideal,*k)?0:alpha)+POLY_EPS){
			if (smpl){
				v_out=poly->cnt;
				add_vrtx (poly);
				for (val=VRTX_VAL(poly,v_out);val<VRTX_VAL(poly,(v_out+1));val++)
					*val=*(VRTX_VAL(poly,(IS_ELEM(poly->ideal,*k)?v:*k))+(val-(VRTX_VAL(poly,v_out))));
				for (val=drctn;val<drctn+poly->dim;val++)
					*val=*(VRTX_VAL(poly,(IS_ELEM(poly->ideal,*k)?*k:v))+(val-drctn));
				if (IS_ELEM(poly->ideal,*k)&&IS_ELEM(poly->ideal,v)){
					for (val=drctn;val<drctn+poly->dim;val++)
						*val-=*(VRTX_VAL(poly,v)+(val-drctn));
					ST_BT(poly->ideal,v_out);
				}else if (!IS_ELEM(poly->ideal, *k)&&!IS_ELEM(poly->ideal,v))
					for (val=drctn;val<drctn+poly->dim;val++)
						*val-=*(VRTX_VAL(poly,*k)+(val-drctn));
				tmp=0;
				for (val=VRTX_VAL(poly,v_out);val<VRTX_VAL(poly,(v_out+1));val++,ct_hp++)
					tmp+=*ct_hp**val;
				ct_hp-=poly->dim;
				
				mu = ((IS_ELEM(poly->ideal,*k)&&IS_ELEM(poly->ideal,v))?0:alpha)-tmp;
				tmp=0;
				for (val=drctn;val<drctn+poly->dim;val++,ct_hp++)
					tmp+=*ct_hp**val;
				ct_hp-=val-drctn;
				mu/=tmp;
				for (val=VRTX_VAL(poly,v_out);val<VRTX_VAL(poly,(v_out+1));val++,drctn++)
					*val+=mu**drctn;
				drctn-=poly->dim;
				add_lst_elem (poly->dual->incidence+poly->dual->cnt-1,v_out);
				add_lst_elem (poly->incidence+v_out,poly->dual->cnt-1);
			}
			for (j=(poly->adjacence+*k)->data;j<(poly->adjacence+*k)->data+(poly->adjacence+*k)->cnt;j++)
				if (*j==v){
					*j=v_out;
					break;
				}
			add_lst_elem(poly->adjacence+v_out,*k);
			for (ks=(poly->incidence+*k)->data;ks<(poly->incidence+*k)->data+(poly->incidence+*k)->cnt;ks++)
				for (js=(poly->incidence+v)->data;js<(poly->incidence+v)->data+(poly->incidence+v)->cnt;js++)
					if (*js==*ks){
						if (!smpl){
							for (ls=(poly->incidence+v_out)->data;ls<(poly->incidence+v_out)->data+(poly->incidence+v_out)->cnt;ls++){
								if (*ls==*ks)
									break;
							}
							if(ls==(poly->incidence+v_out)->data+(poly->incidence+v_out)->cnt){
								add_lst_elem(poly->incidence+v_out,*ks);
								for (l = (poly->dual->incidence+*ks)->data;l<(poly->dual->incidence+*ks)->data+(poly->dual->incidence+*ks)->cnt;l++){
									if (*l==v){
										*l=v_out;
										break;
									}
								}
								if (l==(poly->dual->incidence+*ks)->data+(poly->dual->incidence+*ks)->cnt)
									add_lst_elem(poly->dual->incidence+*ks,v_out);
							}
						}else{
							add_lst_elem(poly->incidence+v_out,*ks);
							for (l = (poly->dual->incidence+*ks)->data;l<(poly->dual->incidence+*ks)->data+(poly->dual->incidence+*ks)->cnt;l++){
								if (*l==v){
									*l=v_out;
									break;
								}
							}
							if (l==(poly->dual->incidence+*ks)->data+(poly->dual->incidence+*ks)->cnt)
								add_lst_elem(poly->dual->incidence+*ks,v_out);
						}
						break;
					}
		}else if (tmp > (IS_ELEM(poly->ideal,*k)?0:alpha)+1.0e-2*POLY_EPS){
			mu=tmp-(IS_ELEM(poly->ideal,*k)?0:alpha);
			tmp=0;
			for (val=ct_hp;val<ct_hp+poly->dim;val++)
				tmp+=*val**val;
			mu /= tmp;
			for (val=VRTX_VAL(poly,*k);val<VRTX_VAL(poly,(*k+1));val++,ct_hp++)
				*val-=mu**ct_hp;
			poly__cut (poly, *k, ct_hp-=poly->dim);
		}else{
			for (j=(poly->adjacence+*k)->data;j<(poly->adjacence+*k)->data+(poly->adjacence+*k)->cnt;j++)
				if (*j==v)
					delete_lst_elem(poly->adjacence+*k,j-(poly->adjacence+*k)->data);
			for (ks=(poly->incidence+*k)->data;ks<(poly->incidence+*k)->data+(poly->incidence+*k)->cnt;ks++){
				for (j=(poly->dual->incidence+*ks)->data;j<(poly->dual->incidence+*ks)->data+(poly->dual->incidence+*ks)->cnt;j++){
					if (*j==v){
						delete_lst_elem(poly->dual->incidence+*ks,j-(poly->dual->incidence+*ks)->data);
						break;
					}
				}
				if ((poly->dual->incidence+*ks)->cnt==0)
					UNST_BT(poly->dual->used,*ks);
			}
			if (IS_ELEM (poly->used,*k))
				poly__cut (poly, *k, ct_hp);
			else
				continue;
		}
	}
	free(vrtx);
	free(drctn);
	for (size_t *fct=(poly->incidence+v)->data;fct<(poly->incidence+v)->data+(poly->incidence+v)->cnt;fct++)
		if ((poly->dual->incidence+*fct)->cnt){
			for (ks=(poly->dual->incidence+*fct)->data;ks<(poly->dual->incidence+*fct)->data+(poly->dual->incidence+*fct)->cnt;ks++)
				if (*ks==v){
					delete_lst_elem(poly->dual->incidence+*fct,ks-(poly->dual->incidence+*fct)->data);
					break;
				}
		}else
			UNST_BT (poly->dual->used,*fct);


	return;
}

void poly__poly_initialise (polytope *poly, double *H, double *R, double *alph, size_t *perm)
{
	double *tmp = (double *)calloc(((poly->dim+2)*(poly->dim+2)-3)/2,sizeof(double));
	double *tmp_dir = (double *)malloc(poly->dim*poly->dim*sizeof(double));
	double *wrk = (double *)calloc(poly->dim,sizeof(double));
	double sum=0,*val,*h;
	for (size_t l=0,j=0;l<poly->dim;l++)
		for (size_t k=0;k<l+1;k++,j++)
			wrk[k]+=R[j];

	for (val=poly->dual->ip;val<poly->dual->ip+poly->dim;val++){
		*val=0;
		for (h=H+(val-poly->dual->ip);h<H+poly->dim*poly->dim;h+=poly->dim)
			*val+=*h**(wrk+(h-(H+(val-poly->dual->ip)))/poly->dim);
		sum+=*val**(poly->ip+(val-poly->dual->ip));

	}
	for (val=poly->dual->ip;val<poly->dual->ip+poly->dim;val++)
		*val/=sum;

	add_vrtx (poly);
	for (size_t k=0;k<poly->dim;k++){
		tmp[k]=alph[k];
		for (size_t l=0;l<poly->dim;l++)
			tmp_dir[l*poly->dim+k]=0;
		tmp_dir[k*poly->dim+k]=1.0;
		add_vrtx (poly);
		ST_BT (poly->ideal,(k+1));
		for (size_t j=0;j<poly->dim;j++)
			*(VRTX_VAL(poly,(k+1))+j)=0;
		*(VRTX_VAL(poly,(k+1))+k)=1.0;
		tmp[poly->dim+(k*(k+1))/2+k]=1.0;
		for (size_t j=0;j<k;j++){
			tmp[k]-=tmp[j]*R[(k*(k+1))/2+j];
			for (size_t l=0;l<poly->dim;l++)
				tmp_dir[l*poly->dim+k]-=tmp_dir[l*poly->dim+j]*R[(k*(k+1))/2+j];
			*(VRTX_VAL(poly,(j+1))+k)-=*(VRTX_VAL(poly,(j+1))+j)*R[(k*(k+1))/2+j];
			tmp[poly->dim+(k*(k+1))/2+j]-=tmp[poly->dim+(j*(j+1))/2+j]*R[(k*(k+1))/2+j];
		}
		tmp[k]/=R[(k*(k+01))/2+k];
		for (size_t l=0;l<poly->dim;l++)
			tmp_dir[l*poly->dim+k]/=R[(k*(k+1))/2+k];
		for (size_t j=0;j<k+1;j++){
			*(VRTX_VAL(poly,(j+1))+k)/=R[(k*(k+1))/2+k];
			tmp[poly->dim+(k*(k+1))/2+j]/=R[(k*(k+1))/2+k];
		}
	}
	for (size_t k=0;k<poly->dim;k++){
		poly->data[k]=0;
		for (size_t j=0;j<poly->dim;j++)
			poly->data[k]+=tmp[j]*H[j*poly->dim+k];
		for (size_t j=0;j<poly->dim;j++){
			*(VRTX_VAL(poly,(k+1))+j)=0;
			for (size_t l=0;l<poly->dim;l++)
				*(VRTX_VAL(poly,(k+1))+j)+=H[l*poly->dim+j]*tmp_dir[k*poly->dim+l];
		}
	}

	for (size_t k=0;k<poly->dim+1;k++)
		{
			for (size_t j=0;j<poly->dim+1;j++)
				{
					if (j!=k)
						{
							add_lst_elem (poly->dual->incidence+*(perm+k),j);
							add_lst_elem (poly->incidence+j,*(perm+k));
							add_lst_elem (poly->adjacence+k,j);
						}
				}
		}
	free(wrk);
	free(tmp);
	free(tmp_dir);


	return;
}

void poly__list_cpy (poly_list *orig,poly_list *dest)
{
	size_t *data;
	dest->cnt=0;
	for (data=orig->data;data<orig->data+orig->cnt;data++)
		add_lst_elem(dest,*data);


	return;
}

void poly__poly_cpy (polytope *orig, polytope *dest)
{
	size_t idx;
	double *val;

	dest->dim=orig->dim;
	dest->blcks=orig->blcks;
	poly__poly_init (dest);
	dest->cnt=orig->cnt;
	for (idx=0;idx<orig->cnt;idx++)
		if (IS_ELEM(orig->used,idx)){
			ST_BT (dest->used,idx);
			for (val=orig->data+idx*orig->dim;val<orig->data+(idx+1)*orig->dim;val++)
				*(dest->data+idx*dest->dim+(val-(orig->data+idx*orig->dim)))=*val;
			if (IS_ELEM(orig->ideal,idx))
				ST_BT (dest->ideal,idx);
			poly__list_cpy (orig->adjacence+idx,dest->adjacence+idx);
			poly__list_cpy (orig->incidence+idx,dest->incidence+idx);
		}


	return;
}


void trnsfrm_polar (double *dual, int is_dir, double *hp)
{
	for (double *val=dual;val<dual+3;val++,hp++)
		*hp=*val;
	*hp=is_dir?0:-1.0;


	return;
}


void poly__swap (poly_args *poly_in,poly_args *poly_out)
{
	double *val;
	size_t *fct;

	for (poly_in->idx=0;poly_in->idx<poly_in->dual.cnt;poly_in->idx++)
		if (IS_ELEM(poly_in->dual.used,poly_in->idx)&&!IS_ELEM(poly_in->dual.ideal,poly_in->idx)){
			for (fct=(poly_in->dual.incidence+poly_in->idx)->data;fct<(poly_in->dual.incidence+poly_in->idx)->data+(poly_in->dual.incidence+poly_in->idx)->cnt;fct++){
				for (val=poly_in->primal.data+*fct*poly_in->dim;val<poly_in->primal.data+(*fct+1)*poly_in->dim;val++,poly_out->val++)
					*poly_out->val=*val;
				poly_out->val-=poly_in->dim;
				poly_out->ideal=IS_ELEM(poly_in->primal.ideal,*fct);
				poly__add_vrtx (poly_out);
			}
			break;
		}

	poly__intl_apprx (poly_out);

	for (poly_in->idx=0;poly_in->idx<poly_in->primal.cnt;poly_in->idx++)
		if (IS_ELEM(poly_in->primal.used,poly_in->idx)){
			for (val=poly_in->primal.data+poly_in->idx*poly_in->dim;val<poly_in->primal.data+(poly_in->idx+1)*poly_in->dim;val++,poly_out->val++)
				*poly_out->val=*val;
			poly_out->val-=poly_in->dim;
			poly_out->ideal=IS_ELEM(poly_in->primal.ideal,poly_in->idx);
			poly__add_vrtx(poly_out);
		}


	return;
}

void poly__plot (polytope *poly, const char *fname)
{
	permutation perm,dual_perm;
	size_t idx,*vrtx,*nghbr,*prm_idx;

	poly_list fct={0,1,NULL};
	FILE *strm;

	poly__initialise_permutation (poly, &perm);
	poly__initialise_permutation (poly->dual, &dual_perm);

	for (vrtx=dual_perm.data;vrtx<dual_perm.data+dual_perm.cnt;vrtx++)
		for (nghbr=vrtx+1;nghbr<dual_perm.data+dual_perm.cnt;nghbr++)
			if (edge_test(poly->dual,*vrtx,*nghbr)){
				add_lst_elem (poly->dual->adjacence+*vrtx,*nghbr);
				add_lst_elem (poly->dual->adjacence+*nghbr,*vrtx);
			}

	strm=fopen (fname,"w");

	fprintf (strm, "OFF\n%zu %zu 0\n\n", perm.cnt, dual_perm.cnt);


	fprintf (strm, "#vertices:\n");
	for (prm_idx=perm.data;prm_idx<perm.data+perm.cnt;prm_idx++){
		for (idx=0;idx<poly->dim;idx++)
			fprintf(strm, "%g ", *(poly->data+*prm_idx*poly->dim+idx));
		fseek (strm, -sizeof(char), SEEK_CUR);
		fprintf (strm, "\n");
	}

	fprintf(strm, "\n#facets:\n");
	fct.data=(size_t *)malloc(LSTBLCK*sizeof(size_t));
	for (prm_idx=dual_perm.data;prm_idx<dual_perm.data+dual_perm.cnt;prm_idx++){
		for (vrtx=(poly->dual->incidence+*prm_idx)->data;vrtx<(poly->dual->incidence+*prm_idx)->data+(poly->dual->incidence+*prm_idx)->cnt;vrtx++)
			add_lst_elem (&fct,*vrtx);
		fprintf (strm, "%zu\t", fct.cnt);
		while (1){
			fprintf (strm, "%zu ", *(perm.inv+*fct.data));
			if (fct.cnt>1){
				for (nghbr=fct.data+1;nghbr<fct.data+fct.cnt;nghbr++){
					for (idx=0;idx<(poly->adjacence+*fct.data)->cnt;idx++)
						if (*((poly->adjacence+*fct.data)->data+idx)==*nghbr)
							break;
					if (idx<(poly->adjacence+*fct.data)->cnt)
						break;
				}
				if (nghbr<fct.data+fct.cnt){
					*fct.data=*nghbr;
					*(fct.data+(nghbr-fct.data))=*(fct.data+--fct.cnt);
				}else{
					fprintf(stderr,"Error: Fault in plot.. exiting\n");
					goto exit;
				}
			}else{
				fct.cnt--;
				break;
			}
		}
		fseek (strm, -sizeof(char), SEEK_CUR);
		fprintf (strm,"\n");
	}

 exit:
	free (fct.data);
	fclose (strm);
	poly__kill_permutation (&perm);
	poly__kill_permutation (&dual_perm);

	return;
}

void poly__polyck (poly_args *poly)
{
	double *hp = (double *)malloc((poly->dim+1)*sizeof(double));
	double scprd,alph,val,eps=1.0e-6;
	size_t *inc,idx,*dual_inc,*nghbr,*nghbr_nghbr,k;

	for (poly->idx=0;poly->idx<poly->dual.cnt;poly->idx++)
		if (IS_ELEM(poly->dual.used,poly->idx)){
			poly->dualV2primalH(poly->dual.data+poly->idx*poly->dim, IS_ELEM(poly->dual.ideal,poly->idx), hp);
			for (inc=(poly->dual.incidence+poly->idx)->data;inc<(poly->dual.incidence+poly->idx)->data+(poly->dual.incidence+poly->idx)->cnt;inc++){
				scprd=0;
				for (idx=0;idx<poly->dim;idx++)
					scprd+=*(hp+idx)**(poly->primal.data+*inc*poly->dim+idx);
				alph = IS_ELEM(poly->primal.ideal,*inc)?0:*(hp+poly->dim);
				val = (scprd-alph<0)?alph-scprd:scprd-alph;
				if (val>eps)
					fprintf(stderr,"Error:\tHyperplane %zu does not contain vertex %zu.\n",poly->idx,*inc);
				for (dual_inc=(poly->primal.incidence+*inc)->data;dual_inc<(poly->primal.incidence+*inc)->data+(poly->primal.incidence+*inc)->cnt;dual_inc++)
					if (*dual_inc==poly->idx)
						break;
				if (dual_inc==(poly->primal.incidence+*inc)->data+(poly->primal.incidence+*inc)->cnt)
					fprintf(stderr,"Error:\tHyperplane %zu, Vertex %zu.\n",poly->idx, *inc);
			}
		}
	for (poly->idx=0;poly->idx<poly->primal.cnt;poly->idx++)
		if (IS_ELEM(poly->primal.used,poly->idx)){
			for (nghbr=(poly->primal.adjacence+poly->idx)->data;nghbr<(poly->primal.adjacence+poly->idx)->data+(poly->primal.adjacence+poly->idx)->cnt;nghbr++){
				for (nghbr_nghbr=(poly->primal.adjacence+*nghbr)->data;nghbr_nghbr<(poly->primal.adjacence+*nghbr)->data+(poly->primal.adjacence+*nghbr)->cnt;nghbr_nghbr++)
					if (*nghbr_nghbr==poly->idx)
						break;
				if (nghbr_nghbr==(poly->primal.adjacence+*nghbr)->data+(poly->primal.adjacence+*nghbr)->cnt)
					fprintf(stderr,"Error:\tVertex %zu appears in vertex' %zu adjacence-list, but not vice versa.\n", *nghbr, poly->idx);
			}
		}
	for (poly->idx=0;poly->idx<poly->primal.cnt;poly->idx++)
		if (IS_ELEM(poly->primal.used,poly->idx))
			for (k=0;k<poly->idx;k++)
				if (IS_ELEM(poly->primal.used,k))
					if (edge_test(&poly->primal,poly->idx,k)){
						for (nghbr=(poly->primal.adjacence+poly->idx)->data;nghbr<(poly->primal.adjacence+poly->idx)->data+(poly->primal.adjacence+poly->idx)->cnt;nghbr++)
							if (*nghbr==k)
								break;
						if (nghbr==(poly->primal.adjacence+poly->idx)->data+(poly->primal.adjacence+poly->idx)->cnt)
							fprintf (stderr,"Error:\tVertices %zu and %zu are adjacent (due to their incidence relation), but vertex %zu does not apper in vertex' %zu adjacence-list.\n", poly->idx, k, k, poly->idx);
					}

	free(hp);

	
	return;
}

void poly__update_adjacence (polytope *poly)
{
	permutation perm;
	size_t *vrtx,*nghbr;


	poly__initialise_permutation (poly, &perm);

	for (vrtx=perm.data;vrtx<perm.data+perm.cnt;vrtx++)
		for (nghbr=vrtx+1;nghbr<perm.data+perm.cnt;nghbr++)
			if (edge_test(poly,*vrtx,*nghbr)){
				add_lst_elem (poly->adjacence+*vrtx,*nghbr);
				add_lst_elem (poly->adjacence+*nghbr,*vrtx);
			}
	poly__kill_permutation (&perm);


	return;
}

void bslv__cpy (double *x, double *y, size_t d)
{
	size_t k;
	for (k=0;k<d;k++,x++,y++)
		*y=*x;
	x-=d;
	y-=d;
}

double bslv__nrm (double *x, size_t d)
{
	double nrm=0,*y;
	for (y=x;y<x+d;y++)
		nrm+=pow(*y,2);

	return sqrt(nrm);
}

double bslv__normalise (double *x, double *X, double *R, size_t k, size_t n)
{
	size_t j,l;
	double sc_prod,scl,*y,nrm_in;
	nrm_in=bslv__nrm(x,n);
	bslv__cpy (x,X+k*n,n);
	for (j=0;j<k;j++){
		scl=0;
		for (l=0;l<n;l++)
			scl += *(X+j*n+l)**(X+k*n+l);
		for (l=0;l<n;l++){
			*(X+k*n+l)-=scl**(X+j*n+l);
		}
	}
	scl = bslv__nrm(X+k*n,n);
	if (scl<1.0e-6)
		return 0;

	for (y=X+k*n;y<X+(k+1)*n;y++)
		*y /= scl;
	scl=scl/nrm_in;
	for (j=0;j<k+1;j++){
		sc_prod = 0;
		for (l=0;l<n;l++)
			sc_prod += *(X+j*n+l)**(x+l);
		*(R+(k*(k+1))/2+j)=sc_prod;
	}


	return scl;
}
