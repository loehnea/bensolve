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

#ifndef BSLV_VLP_H
#define BSLV_VLP_H

#include <setjmp.h> // jmp_buf

#include "bslv_main.h"		/* lp_idx */
#include "bslv_lists.h"
#include "bslv_lp.h"


#define PROBLEM_DESIGNATOR "vlp"
#define STOP_AT_WARNING 0	// set 1 to compile with stop at warning

typedef struct
{
	jmp_buf jump;		// label for go to in case of error
	const char *fname;	// name of input text file
	FILE *fp;			// stream assigned to input text file */
	int count;			// line count
	int c;				// current character
	char field[255+1];	// data field
	int empty;			// warning 'empty line ignored' was printed
	int nonint;			// warning 'non-integer data detected' was printed
	char msg[255+1];	// error and warning messages
	int error;			// error flag
	int warning;		// warning flag
} csatype; 

typedef struct 
{
	list2d *A_ext;					// non-zero constraint and objective coefficients: A_ext = (A,0;P -I)
	boundlist *rows;				// non-standard row bounds (standard is 'f')
	boundlist *cols;				// non-standard column bounds (standard is 's')
	int optdir;						// 1 for minimization, -1 for maximization
	cone_gen_type cone_gen;			// type of ordering cone generators CONE: CONE | DUALCONE | DEFAULT
	double *gen;					// generators of ordering cone (primal or dual)
	double *c;						// duality parameter vector (given data, unscaled)
	long int nz;					// number of non-zero entries of A
	long int nzobj;					// number of non-zero entries of P
	lp_idx n;						// number of variables (cols)
	lp_idx m;						// number of constraints (rows)
	lp_idx q;						// number of objectives
	lp_idx n_gen;					// number of generators of ordering cone (primal or dual)
} vlptype;

typedef struct 
{
	lp_idx m;						// number of rows (constraints)
	lp_idx n;						// number of cols (varibles)
	lp_idx q;						// number of objectives
	lp_idx o;						// number of generators of ordering cone (after vertex enumeration and scaling)
	lp_idx p;						// number of generators of dual of ordering cone (after vertex enumeration and scaling)
	lp_idx r;						// number of generators of dual cone of recession cone of upper image
	lp_idx h;						// number of generators of recession cone of upper image
	double *eta;					// result of phase0
	double *Y;						// generators of ordering cone as columns of Y (non-redundant and scaled columns)
	double *Z;						// generators of dual cone of C as columns in matrix Z (non-redundant and scaled columns that that Z' * c == (1,...,1)')
	double *c;						// geometric duality parameter vector (scaled such that c_q=1)
	double *R;						// result of phase1: columns are generators of dual cone of recession cone of upper image
	double *H;						// result of phase1: columns are generators of recession cone of upper image
	sol_status_type status;			// solution status of VLP
	c_dir_type c_dir;				// type of duality parameter vector
	size_t pp;						// number of vertices of upper image
	size_t dd;						// number of vertices of lower image
	size_t pp_dir;					// number of extreme directions of upper image
	size_t dd_dir;					// number of extreme directions of lower image
} soltype;

typedef struct 
{
	int bounded;
	int plot;
	char filename[255+1];
	pre_img_type solution;					// PRE_IMG_OFF - PRE_IMG_ON
	format_type format;						// SHORT - LONG - AUTO
	int message_level;						// 0 - 1 - 2 - 3
	alg_type alg_phase1;					// PRIMAL - DUAL
	alg_type alg_phase2;					// PRIMAL - DUAL
	double eps_phase0;						// Epsilon used in Phase 0
	double eps_phase1;						// Epsilon used in Phase 1
	double eps_benson_phase1;				// Epsilon of Benson algorithm in Phase 1
	double eps_benson_phase2;				// Epsilon of Benson algorithm in Phase 2
    struct lp_opt lp;
} opttype;

int vlp_init(char const *filename, vlptype *vlp, const opttype *opt);
void vlp_free(vlptype *vlp);

void sol_init(soltype *sol, const vlptype *vlp, const opttype *opt);
void sol_free(soltype *sol);

void set_default_opt(opttype *opt);
#endif
