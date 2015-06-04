/*
This file is part of BENSOLVE - VLP solver

Copyright © 2014-2015 Andreas Löhne and Benjamin Weißing

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

#ifndef BSLV_POLY_H
#define BSLV_POLY_H


#include <limits.h>
#include <stdio.h>


#ifndef ALLOCFCTR
#define ALLOCFCTR 1
#endif
#ifndef VRTXBLCK
#define VRTXBLCK (ALLOCFCTR*BTCNT)
#endif
#ifndef LSTBLCK
#define LSTBLCK 1
#endif

#define VRTX_VAL(poly,idx) (poly)->data+idx*(poly)->dim

typedef size_t btstrg;
typedef btstrg vrtx_strg;
#define BTCNT (CHAR_BIT*sizeof(btstrg))
#define ST_BT(lst,idx) (*(lst+idx/BTCNT)|=(btstrg)1<<idx%BTCNT)
#define UNST_BT(lst,idx) (*(lst+idx/BTCNT)&=~((btstrg)1<<idx%BTCNT))
#define IS_ELEM(lst,idx) ((btstrg)1U&(*(lst+idx/BTCNT)>>idx%BTCNT))

#define POLY_EPS 1e-9

typedef struct poly_list_strct{
	size_t cnt;
	size_t blcks;
	size_t *data;
} poly_list;

typedef struct polytope_strct{
	size_t dim,dim_primg;
	size_t cnt;
	size_t blcks;
	double *ip;
	double *data;
	double *data_primg;
	poly_list *adjacence;
	poly_list *incidence;
	vrtx_strg *ideal;
	vrtx_strg *used;
	vrtx_strg *sltn;
	struct polytope_strct *dual;
	void (*v2h)(double *, int, double *);
} polytope;

typedef struct{
	size_t dim,dim_primg_prml,dim_primg_dl;
	unsigned int ideal:1;
	size_t idx;
	double *val,*val_primg_prml,*val_primg_dl;
	double eps;
	polytope primal;
	polytope dual;
	void (*primalV2dualH)();
	void (*dualV2primalH)();
	struct{double *H,*R,*alph;poly_list queue,gnrtrs;unsigned int intlsd:1;} init_data;
} poly_args;

typedef struct {
	size_t cnt;
	size_t *data;
	size_t *inv;
} permutation;

void poly__initialise_permutation (polytope *, permutation *);
void poly__kill_permutation (permutation *);
void poly__vrtx2file (polytope *, permutation *, const char *, const char *);
void poly__primg2file (polytope *, permutation *, const char *, const char *);
void poly__adj2file (polytope *, permutation *, const char *, const char *);
void poly__inc2file (polytope *, permutation *, permutation *, const char *, const char *);

void poly__set_default_args (poly_args *args, size_t dim);
void poly__initialise (poly_args *);
void poly__kill (poly_args *);
void poly__cut (polytope *, size_t, double *);
void poly__poly_initialise (polytope *, double *,double *,double *,size_t *);
int poly__add_vrtx (poly_args *);
int poly__get_vrtx (poly_args *);

void poly__poly_init (polytope *);
void poly__poly_kill (polytope *);
void poly__list_init (poly_list *);
void add_vrtx (polytope *);
int poly__intl_apprx (poly_args *);
void add_lst_elem (poly_list *, size_t);
int edge_test (polytope *,size_t,size_t);
void poly__update_adjacence (polytope *);
void vrtx_cpy (polytope *, size_t,size_t);
void poly__swap (poly_args *, poly_args *);
void poly__plot (polytope *, const char *);
void poly__polyck (poly_args *poly);

double bslv__normalise (double *, double *, double *, size_t, size_t);

#endif
