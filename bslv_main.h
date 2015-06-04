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

#ifndef BSLV_MAIN_H
#define BSLV_MAIN_H

typedef int lp_idx;

#define ABS(x) (((x) < 0 ? -(x) : (x)))
#define MAX(x,y) (((x) < (y) ? (y) : (x)))

#define THISVERSION "version 2.0.1"
#define WELCOME "BENSOLVE: VLP Solver, %s\nCopyright (C) 2014-2015 Andreas Löhne and Benjamin Weißing\nThis is free software; see the source code for copying conditions.\nThere is ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY or\nFITNESS FOR A PARTICULAR PURPOSE.\n"
#define USAGE "\nUsage: bensolve file [options]\n\n"

#ifndef POLY_TEST
#define POLY_TEST 0
#endif

#ifndef EPS_C
#define EPS_C 1e-7
#endif

#ifndef EPS_POLY
#define EPS_POLY 1e-9
#endif

#ifndef EPS_OUTPUT_CHOP
#define EPS_OUTPUT_CHOP 1e-10
#endif

#ifndef PRIMAL_PLOT_CUT_SHIFT
#define PRIMAL_PLOT_CUT_SHIFT 0.5
#endif

#ifndef DUAL_PLOT_CUT_SHIFT
#define DUAL_PLOT_CUT_SHIFT 0.2
#endif

#ifndef FORMAT_SHORT_STR
#define FORMAT_SHORT_STR "%10.4g "
#endif

#ifndef FORMAT_LONG_STR
#define FORMAT_LONG_STR "%.14g "
#endif

#define DEFAULT_MESSAGE_LEVEL 1
#define DEFAULT_MESSAGE_LEVEL_STR "1"

#define DEFAULT_LP_MESSAGE_LEVEL		1
#define DEFAULT_LP_MESSAGE_LEVEL_STR	"1"

#define DEFAULT_EPS_PHASE0 1e-8
#define DEFAULT_EPS_PHASE1 1e-8

#define DEFAULT_EPS_BENSON_PHASE1 1e-7
#define DEFAULT_EPS_BENSON_PHASE1_STR "1e-7"

#define DEFAULT_EPS_BENSON_PHASE2 1e-7
#define DEFAULT_EPS_BENSON_PHASE2_STR "1e-7"

#define MIN_P_STR "Upper image of primal problem:\n"
#define MIN_D_STR "Lower image of dual problem:\n"
#define MAX_P_STR "Lower image of primal problem:\n"
#define MAX_D_STR "Upper image of dual problem:\n"
#define CONE_P_STR "Ordering cone:\n"
#define CONE_D_STR "Dual of ordering cone:\n"
#define CONE_ENDING_STR ".cone"
#define SOL_ENDING_STR ".sol"

#define MAX_STR_LNGTH 10 // adopt string length when changing these strings
#define PRE_IMG_P_STR "_pre_img_p"
#define PRE_IMG_D_STR "_pre_img_d"
#define IMG_P_STR "_img_p"
#define IMG_D_STR "_img_d"
#define ADJ_P_STR "_adj_p"
#define ADJ_D_STR "_adj_d"
#define INC_P_STR "_inc_p"
#define INC_D_STR "_inc_d"

enum _alg_type {PRIMAL_BENSON, DUAL_BENSON};
enum _lp_method_type {PRIMAL_SIMPLEX, DUAL_SIMPLEX, DUAL_PRIMAL_SIMPLEX, LP_METHOD_AUTO};
enum _lp_hom_type {HOMOGENEOUS, INHOMOGENEOUS};
enum _cone_out_type {CONE_OUT_OFF, CONE_OUT_ON};
enum _phase_type {PHASE0, PHASE1_PRIMAL, PHASE1_DUAL, PHASE2_PRIMAL, PHASE2_DUAL};
enum _format_type {FORMAT_SHORT, FORMAT_LONG, FORMAT_AUTO};
enum _lp_status_type {LP_INFEASIBLE, LP_UNBOUNDED, LP_UNEXPECTED_STATUS, LP_UNDEFINED_STATUS, LP_OPTIMAL};
enum _sol_status_type {VLP_NOSTATUS, VLP_INFEASIBLE, VLP_UNBOUNDED, VLP_NOVERTEX, VLP_OPTIMAL, VLP_INPUTERROR};
enum _cone_gen_type {CONE, DUALCONE, DEFAULT};
enum _c_dir_type{C_DIR_POS, C_DIR_NEG};
enum _swap_type{SWAP, NO_SWAP};
enum _pre_img_type{PRE_IMG_OFF, PRE_IMG_ON}; 

typedef enum _alg_type alg_type;
typedef enum _lp_method_type lp_method_type;
typedef enum _lp_hom_type lp_hom_type;
typedef enum _cone_out_type cone_out_type;
typedef enum _phase_type phase_type;
typedef enum _format_type format_type;
typedef enum _lp_status_type lp_status_type;
typedef enum _sol_status_type sol_status_type;
typedef enum _cone_gen_type cone_gen_type;
typedef enum _c_dir_type c_dir_type;
typedef enum _swap_type swap_type;
typedef enum _pre_img_type pre_img_type;

#endif
