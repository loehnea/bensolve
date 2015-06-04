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

#ifndef BSLV_ALGS_H
#define BSLV_ALGS_H

#include "bslv_vlp.h"

int cone_vertenum(double **prim, lp_idx *n_prim, double **dual, lp_idx *n_dual, double *prim_in, const size_t n_prim_in, const size_t dim, const opttype *opt, cone_out_type output, swap_type swap);

void phase0(soltype *const sol, const vlptype *vlp, const opttype *opt);
void phase1_primal(soltype *const sol, const vlptype *vlp, const opttype *opt);
void phase2_primal(soltype *const sol, const vlptype *vlp, const opttype *opt);
void phase1_dual(soltype *const sol, const vlptype *vlp, const opttype *opt);
void phase2_dual(soltype *const sol, const vlptype *vlp, const opttype *opt);
void phase2_init(soltype *sol, const vlptype *vlp);

#endif
