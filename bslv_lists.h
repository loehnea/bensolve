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

#ifndef BSLV_LISTS_H
#define BSLV_LISTS_H

#include "bslv_main.h"
#include "bslv_poly.h"

typedef struct
{
	lp_idx size;	// number of entries
	lp_idx *idx;	// vector of indices (smallest index is 1)
	double *data;	// vector of coefficients
} list1d;

typedef struct
{
	size_t size;	// number of entries
	lp_idx *idx1;	// vector of first (row) indices (smallest index is 1)
	lp_idx *idx2;	// vector of second (column) indices (smallest index is 1)
	double *data;	// vector of coefficients
} list2d;

typedef struct
{ 
	lp_idx size;	// number of entries	
	lp_idx *idx;	// vector of indices (smallest index is 1)
	double *lb;		// vector of lower bound entries
	double *ub;		// vector of upper bound entries
	char *type;		// vector of GLPK types: 'f', 'l', 'u', 'd', 's'
} boundlist;

void string_fprint(const char *filename, const char *string);
void matrix_fprint (double *mat_arr, int m, int n, int tda, char *filename, const char *format);
void matrix_print (double *mat_arr, int m, int n, const char *format);
int string_to_int(char *str, char *error_msg);
int string_to_positive_int(char *str, char *error_msg);
double string_to_positive_double(char *str, char *error_msg);
void orthogonal_vector(double * mat_arr, int dim, int cidx);
int is_equal(const lp_idx size, const double *vec1, const double *vec2, const double tol);
int is_zero(const lp_idx size, const double *vec, const double tol);

// list1d

list1d *list1d_alloc(lp_idx size);
list1d *list1d_calloc(lp_idx size);
void list1d_init_idx(list1d *list, lp_idx firstidx);
void list1d_free(list1d *list);
void vector_to_list1d(list1d *const list, const double *vec_arr, int n);
void list1d_print (list1d const *list, int size);

// list2d

list2d *list2d_alloc(size_t size);
list2d *list2d_calloc(size_t size);
void list2d_init_idx(list2d *list, lp_idx nrows, lp_idx ncols);
void list2d_free(list2d *list);
void list2d_print (list2d const *list);

// boundlist

boundlist *boundlist_alloc(lp_idx size);
boundlist *boundlist_calloc(lp_idx size, char type);
void boundlist_init_idx(boundlist *list, lp_idx firstidx);
void boundlist_free(boundlist *list);
void boundlist_print (boundlist const *list);

#endif
