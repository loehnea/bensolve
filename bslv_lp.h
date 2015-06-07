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

#ifndef BSLV_LP_H
#define BSLV_LP_H

#include "bslv_main.h"
#include "bslv_lists.h"

void lp_write (size_t i);

double lp_write_sol(size_t i);

// initialize lp
void lp_init (int, int, int, lp_idx *, lp_idx *, double *);
/*
 - for nonapprearing coefficients in obj and A, zero is assumed
 - for non-appearing rows, GLPK standard type 'f' is assumed
 - for non-appearing colums, GLPK standard type 's' is assumed

 - obj index range: 0 ... ncols, where 0 stands for a 'shift'
 - A index range: (1,1) ... (nrows,ncols)
 - rows index range: 1 ... nrows
 - cols index range: 1 ... ncols
*/

// set lp options

typedef enum {PRIMAL_SIMPLEX, DUAL_SIMPLEX, DUAL_PRIMAL_SIMPLEX, LP_METHOD_AUTO} lp_method_type;
typedef enum {LP_INFEASIBLE, LP_UNBOUNDED, LP_UNEXPECTED_STATUS, LP_UNDEFINED_STATUS, LP_OPTIMAL} lp_status_type;
typedef enum {HOMOGENEOUS, INHOMOGENEOUS} lp_hom_type;

struct lp_opt {
  lp_method_type method_phase0,method_phase1,method_phase2;
  int message_level;
};

void lp_set_options(const struct lp_opt *opt, phase_type phase);

// create a copy lp[dest]] of lp[src]
void lp_copy(size_t dest, size_t src);

// delete extra rows and set num new extra rows
void lp_update_extra_coeffs (lp_idx n_rows, lp_idx n_cols);

// set (replace) row bounds
void lp_set_rows (size_t i, boundlist const *rows);

// set (replace) row bounds of homogeneous problem
void lp_set_rows_hom (size_t i, boundlist const *rows);

// set (replace) column bounds
void lp_set_cols (size_t i, boundlist const *cols);

// set (replace) column bounds of homogeneous problem
void lp_set_cols_hom (size_t i, boundlist const *cols);

// set (replace) constraint coefficients
int lp_set_mat (size_t i, list2d const *A);

// set (replace) constraint coefficients row
void lp_set_mat_row (size_t i, list1d *list, lp_idx ridx);

// set all objective coefficients to zero
void lp_clear_obj_coeffs (size_t i);

// set (replace) (a subset of) objective coefficients
void lp_set_obj_coeffs (size_t i, list1d const *list);

// solve problem, return
lp_status_type lp_solve(size_t i);

// retrieve solutions, x and v need to be allocated before functions are called */
void lp_primal_solution_rows(size_t i, double *const x, lp_idx firstidx, lp_idx size, double sign);
void lp_primal_solution_cols(size_t i, double *const x, lp_idx firstidx, lp_idx size, double sign);
void lp_dual_solution_rows(size_t i, double *const u, lp_idx firstidx, lp_idx size, double sign);
void lp_dual_solution_cols(size_t i, double *const u, lp_idx firstidx, lp_idx size, double sign);

// return (optimal) objective value
double lp_obj_val(size_t i);

// return CPU time of lp solver in seconds
double lp_get_time (size_t i);

// return number of LPs (type i) solved
int lp_get_num (size_t i);

void lp_free(size_t i);

#endif
