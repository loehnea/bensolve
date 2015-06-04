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

#include <string.h> // memcpy
#include <limits.h> // INT_MIN, INT_MAX
#include <assert.h> // assert
#include <math.h>
#include <stdlib.h>

#include "bslv_lists.h"

void string_fprint(const char *filename, const char *str)
{
	FILE *fp = fopen (filename, "w+");
	if (fp == NULL)	
	{
		printf("unable to open file %s", filename);
		exit(1);
	}
	fprintf(fp, "%s", str);
	fclose(fp);
}

int string_to_positive_int(char *str, char *error_msg)
{
	long val;
	char *end_ptr; 
	val = strtol(str, &end_ptr, 10);
	if (end_ptr == str || val > INT_MAX || val < 1 || '\0' != *end_ptr)
	{
		printf("%s", error_msg);
		exit(1);
	}
	else
		return (int) val;
}

int string_to_int(char *str, char *error_msg)
{
	long val;
	char *end_ptr;
	val = strtol(str, &end_ptr, 10);
	if (end_ptr == str || val > INT_MAX || val < INT_MIN || '\0' != *end_ptr)
	{	
		printf("%s", error_msg);
		exit(1);
	}
	else
		return (int) val;
}

double string_to_positive_double(char *str, char *error_msg)
{
	double val;
	char *end_ptr;
	val = strtod(str, &end_ptr);
	if (end_ptr == str || '\0' != *end_ptr || val <= 0)
	{
		printf("%s", error_msg);
		exit(1);
	}
	else
		return val;
}

void matrix_fprint (double *mat, int m, int n, int tda, char *filename, const char *format)
{
	FILE *fp;
	fp = fopen (filename, "w+");
	if (fp == NULL)	
	{
		printf("unable to open file %s", filename);
		exit(1);
	}
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < n; j++)
			fprintf(fp,format?format:"%g ", ABS(mat[i * tda + j])<1e-14?0:mat[i * tda + j]);
		fseek (fp, -sizeof(char), SEEK_CUR);
		fprintf(fp,"\n");
	}
	fclose(fp);
}

void matrix_print (double *mat_arr, int m, int n, const char *format)
{
	for (size_t i = 0; i < m; i++)
	{
		for (size_t j=0; j< n; j++)
			printf(format?format:"%g ", ABS(mat_arr[i * n + j])<1e-14?0:mat_arr[i * n + j]);
		printf("\n");
	}
}

void orthogonal_vector(double *C, int dim, int i)
{
	// replace column i by a vector which is orthogonal to (orthogonal) columns 0,...,i-1
	double scl, scl1;

	for (size_t r=0; r<dim;r++)
	{
		for (size_t k=0; k<dim; k++)
			C[k*dim+i]=0;
		C[((i+r)%dim)*dim+i]=1; // try all unit vectors, start with e^i

		for (size_t j=0; j<i; j++)
		{
			scl=0;
			scl1=0;
			for (size_t k=0; k<dim; k++)
				scl += C[k*dim+j]*C[k*dim+i];
			for (size_t k=0; k<dim; k++)
				scl1 += C[k*dim+j]*C[k*dim+j];
			for (size_t k=0; k<dim; k++)
				C[k*dim+i]-=scl/scl1*C[k*dim+j];
		}
		scl=0;
		for (size_t k=0; k<dim; k++)
			scl += C[k*dim+i]*C[k*dim+i];
		if (scl > 1e-3)
			break;
	}
	for (size_t k=0; k<dim; k++)
		C[k*dim+i] *= pow(scl,-0.5);
}

int is_equal(const lp_idx size, const double *vec1, const double *vec2, const double tol)
{
	double tmp = 0;
	double tol_sq = tol*tol;
	for (lp_idx i = 0; i < size; i++)
	{
		tmp += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
		if (tmp > tol_sq)
			return 0;
	}
	return 1;
}

int is_zero(const lp_idx size, const double *vec, const double tol)
{
	double tmp = 0;
	double tol_sq = tol*tol;
	for (lp_idx i = 0; i < size; i++)
	{
		tmp += vec[i] * vec[i];
		if (tmp > tol_sq)
			return 0;
	}
	return 1;
}

// list1d

list1d *list1d_alloc(lp_idx size)
{
	list1d *list;
    list = malloc(sizeof(list1d));
	list->size = size;
	list->idx = (lp_idx*) malloc (sizeof(lp_idx)*size);
	list->data = (double*) malloc (sizeof(double)*size);
	return list;
}

list1d *list1d_calloc(lp_idx size)
{
	list1d *list;
	list = malloc(sizeof(list1d));
	list->size = size;
	list->idx = (lp_idx*) calloc(size, sizeof(lp_idx));
	list->data = (double*) calloc(size, sizeof(double));
	return list;
}

void list1d_init_idx(list1d *list, lp_idx firstidx)
{
	for(lp_idx k = 0; k < list->size; k++)
			list->idx[k] = k + firstidx;
}

void vector_to_list1d(list1d *const list, const double *vec_arr, int n)
{
	for(lp_idx i = 0; i < n; i++)
		list->data[i] = vec_arr[i];
}	

void list1d_free(list1d *list)
{
	free(list->idx);
	free(list->data);
	free(list);
}

// list2d

list2d *list2d_alloc(size_t size)
{
	list2d *list;
	list = malloc(sizeof(list2d));
	list->size = size;
	list->idx1 = (lp_idx*) malloc (sizeof(lp_idx)*size);
	list->idx2 = (lp_idx*) malloc (sizeof(lp_idx)*size);
	list->data = (double*) malloc (sizeof(double)*size);
	return list;
}

list2d *list2d_calloc(size_t size)
{
	list2d *list;
	list = malloc(sizeof(list2d));
	list->size = size;
	list->idx1 = (lp_idx*) calloc(size, sizeof(lp_idx));
	list->idx2 = (lp_idx*) calloc(size, sizeof(lp_idx));
	list->data = (double*) calloc(size, sizeof(double));
	return list;
}

void list2d_init_idx(list2d *list, lp_idx nrows, lp_idx ncols)
{
	if (list->size < nrows * ncols) {
		printf(	"Error:"
				" list2d_init_idx(list2d *list, int nrows, int ncols):"
				" index out of scope\n");
	}		
	for(lp_idx i = 0; i < nrows; i++){
		for(lp_idx j = 0; j < ncols; j++){
			list->idx1[j + ncols * i] = i+1;
			list->idx2[j + ncols * i] = j+1;
		}
	}
}

void list2d_free(list2d *list)
{
	free(list->idx1);
	free(list->idx2);
	free(list->data);
	free(list);
}

void list2d_print (list2d const *list)
{
	printf("Size: %zd\n",list->size);
	for (size_t k = 0; k < list->size; k++)
		printf("%3d %3d %7.5g\n", list->idx1[k], list->idx2[k], list->data[k]);
}

// boundlist

boundlist *boundlist_alloc(lp_idx size)
{
	boundlist *list;
	list = malloc(sizeof(boundlist));
	list->size = size;
	list->idx = (lp_idx*) malloc (sizeof(lp_idx)*size);
	list->lb = (double*) malloc (sizeof(double)*size);
	list->ub = (double*) malloc (sizeof(double)*size);
	list->type = (char*) malloc (sizeof(char)*size);
	return list;
}

boundlist *boundlist_calloc(lp_idx size, char type)
{
	boundlist *list;
	list = malloc(sizeof(boundlist));
	list->size = size;
	list->idx = (lp_idx*) calloc(size, sizeof(lp_idx));
	list->lb = (double*) calloc(size, sizeof(double));
	list->ub = (double*) calloc(size, sizeof(double));
	list->type = (char*) malloc (sizeof(char)*size);
	for (lp_idx k=0; k<size; k++) list->type[k] = type;
	return list;
}

void boundlist_init_idx(boundlist *list, lp_idx firstidx)
{
	for(lp_idx k = 0; k < list->size; k++)
			list->idx[k] = k+firstidx;
}

void boundlist_free(boundlist *list)
{
	free(list->idx);
	free(list->lb);
	free(list->ub);
	free(list->type);
	free(list);
}
