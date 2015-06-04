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

#include <sys/time.h>	// for gettimeofday()
#include <getopt.h>
#include <string.h>
#include <assert.h>
//#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "bslv_vlp.h"
#include "bslv_lp.h"
#include "bslv_lists.h"
#include "bslv_algs.h"

struct timeval t_start,t_end;

int main(int argc, char **argv)
{
	const char *options_string = 
	"  --help, -h             Print short help message and exit.\n"
	"  --bounded, -b          Assume that the problem is bounded. Skip phases 0 and 1.\n"
	"  --plot, -p             Generate an OFF graphics file of upper and lower images.\n"
	"  --test, -t             Run integrity tests for polytopes.\n"
	"  --solution, -s         Write primal and dual solution to files\n"
	"  --format, -f           Choose a output format, args: long - short (default auto: short on screen, long in files)\n"
	"  --output_filename, -o  Use alternative filename for output, arg: [filename]\n"
	"  --lp_method_phase0, -k Choose lp method for phase 0, ...\n"
	"  --lp_method_phase1, -L Choose lp method for phase 1, ...\n"
	"  --lp_method_phase2, -l Choose lp method for phase 2, ...\n"
	"                         ... args: primal_simplex - dual_simplex - dual_primal_simplex (default: auto)\n"	
	"  --message_level, -m    Display less or more messages, args: 0 - 1 - 2 - 3 (default: " DEFAULT_MESSAGE_LEVEL_STR ")\n"
	"  --lp_message_level, -M Display less or more messages of the LP solver, args: 0 - 1 - 2 - 3 (default: " DEFAULT_LP_MESSAGE_LEVEL_STR ")\n"
	"  --alg_phase1, -A       Choose algorithm type for phase 1, args: primal (default) - dual\n"
	"  --alg_phase2, -a       Choose algorithm type for phase 2, args: primal (default) - dual\n"
	"  --eps_phase1, -E       Determine epsilon used in phase 1, example: '-E 0.01', default: '-e " DEFAULT_EPS_BENSON_PHASE1_STR "'\n"
	"  --eps_phase2, -e       Determine epsilon used in phase 2, example: '-e 0.01', default: '-e " DEFAULT_EPS_BENSON_PHASE1_STR "'\n"
	" \n"
	" For additional information, see the reference manual: doc/manual.pdf\n\n"
	" To contact the authors, visit http://www.bensolve.org\n\n";

	char *filename = *(argv+1);

	if(argc<=1 || filename[0]=='-')
	{
		printf(WELCOME, THISVERSION);
		printf(USAGE);
		printf("%s",options_string);
		return 1;
	}

	opttype _opt, *opt = &_opt;
	set_default_opt(opt);
	/*
	 *
	 * set options
	 *
	 */
	int c;
	for(;;)
	{
		static struct option long_options[] =
		{
			{"help", no_argument, 0, 'h'},
			{"bounded", no_argument, 0, 'b'},
			{"plot", no_argument, 0, 'p'},
			{"solution", no_argument, 0, 's'},
			{"format", required_argument, 0, 'f'},				// short - long - mixed
			{"output_filename", required_argument, 0, 'o'},		// [filename]
			{"lp_method_phase0", required_argument, 0, 'k'},	// primal_simplex - dual_simplex - dual_primal_simplex
			{"lp_method_phase1", required_argument, 0, 'L'},	// primal_simplex - dual_simplex - dual_primal_simplex
			{"lp_method_phase2", required_argument, 0, 'l'},	// primal_simplex - dual_simplex - dual_primal_simplex
			{"message_level", required_argument, 0, 'm'},		// 0 - 1 - 2 - 3
			{"lp_message_level", required_argument, 0, 'M'},	// 0 - 1 - 2 - 3
			{"alg_phase1", required_argument, 0, 'A'},			// primal - dual
			{"alg_phase2", required_argument, 0, 'a'},			// primal - dual
			{"eps_phase1", required_argument, 0, 'E'},			// [epsilon]
			{"eps_phase2", required_argument, 0, 'e'},			// [epsilon]
			{0, 0, 0, 0}
		};
		int option_index = 0;
		c = getopt_long (argc, argv, "hbptsf:o:k:L:l:m:M:A:a:E:e:", long_options, &option_index);
		if (c == -1) break; // detect end of options
		switch (c)
		{
			case 'h':
				printf(WELCOME, THISVERSION);
				printf(USAGE);
				printf("%s",options_string);
				return 1;
			case 'b':
				opt->bounded = 1;
				break;
			case 'p':
				opt->plot = 1;
				break;
			case 's':
				opt->solution = PRE_IMG_ON;
				break;
			case 'f':
				if (strcmp(optarg,"auto") == 0) opt->format = FORMAT_AUTO;
				else if (strcmp(optarg,"long") == 0) opt->format = FORMAT_LONG;
				else if (strcmp(optarg,"short") == 0) opt->format = FORMAT_SHORT;
				else
				{
					printf("option --format (-f): invalid argument\n");
					exit(1);
				}
				break;
			case 'o':
				strcpy(opt->filename, optarg);
				break;
			case 'k':
				if (strcmp(optarg,"primal_simplex") == 0) opt->lp_method_phase0 = PRIMAL_SIMPLEX;
				else if (strcmp(optarg,"dual_simplex") == 0) opt->lp_method_phase0 = DUAL_SIMPLEX;
				else if (strcmp(optarg,"dual_primal_simplex") == 0) opt->lp_method_phase0 = DUAL_PRIMAL_SIMPLEX;
				else
				{
					printf("option --lp_method_phase0 (-k): invalid argument\n");
					exit(1);
				}
				break;
			case 'L':
				if (strcmp(optarg,"primal_simplex") == 0) opt->lp_method_phase1 = PRIMAL_SIMPLEX;
				else if (strcmp(optarg,"dual_simplex") == 0) opt->lp_method_phase1 = DUAL_SIMPLEX;
				else if (strcmp(optarg,"dual_primal_simplex") == 0) opt->lp_method_phase1 = DUAL_PRIMAL_SIMPLEX;
				else if (strcmp(optarg,"auto") == 0) opt->lp_method_phase1 = LP_METHOD_AUTO;
				else
				{
					printf("option --lp_method_phase1 (-L): invalid argument\n");
					exit(1);
				}
				break;
			case 'l':
				if (strcmp(optarg,"primal_simplex") == 0) opt->lp_method_phase2 = PRIMAL_SIMPLEX;
				else if (strcmp(optarg,"dual_simplex") == 0) opt->lp_method_phase2 = DUAL_SIMPLEX;
				else if (strcmp(optarg,"dual_primal_simplex") == 0) opt->lp_method_phase2 = DUAL_PRIMAL_SIMPLEX;
				else if (strcmp(optarg,"auto") == 0) opt->lp_method_phase2 = LP_METHOD_AUTO;
				else
				{
					printf("option --lp_method_phase2 (-l): invalid argument\n");
					exit(1);
				}
				break;
			case 'M':
				opt->lp_message_level = string_to_int(optarg, "option --lp_message_level (-M): invalid argument\n");
				if (opt->lp_message_level>3 || opt->lp_message_level<0)
				{
 					printf("option --lp_message_level (-M): invalid argument\n");
 					exit(1);
				}
				break; 
			case 'm':
				opt->message_level = string_to_int(optarg, "option --message_level (-m): invalid argument\n");
				if (opt->message_level>3 || opt->message_level<0)
				{
					printf("option --message_level (-m): invalid argument\n\n");
					exit(1);
				}
				break;
			case 'A':
				if (strcmp(optarg,"primal") == 0) opt->alg_phase1 = PRIMAL_BENSON;
				else if (strcmp(optarg,"dual") == 0) opt->alg_phase1 = DUAL_BENSON;
				else
				{
					printf("option --alg_phase1 (-A): invalid argument\n");
					exit(1);
				}
				break;
			case 'a':
				if (strcmp(optarg,"primal") == 0) opt->alg_phase2 = PRIMAL_BENSON;
				else if (strcmp(optarg,"dual") == 0) opt->alg_phase2 = DUAL_BENSON;
				else
				{
					printf("option --alg_phase2 (-a): invalid argument\n");
					exit(1);
				}
				break;
			case 'E':
				opt->eps_benson_phase1 = string_to_positive_double(optarg, "option --eps_benson_phase1 (-E): invalid argument\n");
				break;
			case 'e':
				opt->eps_benson_phase2 = string_to_positive_double(optarg, "option --eps_benson_phase2 (-e): invalid argument\n");
				break;
			case '?':
				// getopt_long prints error message.
				return 1;
				default:
				printf("invalid option -%c\n", c);
				return 1;
			}
	}
	
	// generate default output filename
	if (strlen(opt->filename) == 0)
	{
		strcpy(opt->filename, filename);
		strtok(opt->filename,".");
	}
	
	if (opt->message_level >= 1)
		 printf(WELCOME, THISVERSION);

	/*
	 *  read problem from file 
	 */
	vlptype _vlp, *vlp = &_vlp;
	if (opt->message_level >= 1)
		 printf("loading ... \n");
	if (vlp_init(filename, vlp, opt)) 
	{	
		vlp_free(vlp);
		return 1;
	}
	if (opt->message_level >= 1) printf("done: %d rows, %d columns, %zd non-zero matrix coefficients\n", vlp->m, vlp->n, vlp->nz);

	// begin of computations - start timer
	gettimeofday(&t_start, NULL);
	
	soltype _sol, *sol = &_sol;
	sol_init(sol,vlp,opt);
	if (sol->status == VLP_INPUTERROR)
	{
		vlp_free(vlp);
		sol_free(sol);
		return 1;
	}
	
	if (opt->plot)
	{
		if (vlp->q!=3)
		{
			printf("OFF file generation for problem with 3 objectives only - try again without option -p\n");
			vlp_free(vlp);
			sol_free(sol);
			return 1;
		}
	}

	lp_init(vlp);
	
	if (opt->message_level == 1)
		 printf("running ... \n");
	
	if (opt->bounded)
		phase2_init(sol,vlp);
	else
	{
		/*
		 * phase 0
		 */
		if (opt->message_level >= 3) printf("\n********* PHASE 0 *********\n");
		phase0(sol, vlp, opt);
		if (sol->status == VLP_UNBOUNDED)
		{
			printf("VLP is totally unbounded, there is no solution\n");
			lp_free(0);
			vlp_free(vlp);
			sol_free(sol);
			return 1;
		}
		if (sol->status == VLP_NOVERTEX)
		{
			printf("upper image of VLP has no vertex (this case is not covered by this version)\n");
			lp_free(0);
			vlp_free(vlp);
			sol_free(sol);
			return 1;
		}
		if (opt->message_level >= 2) { printf("Result of phase 0: eta = \n  "); matrix_print(sol->eta, 1, vlp->q, opt->format==FORMAT_LONG?FORMAT_LONG_STR:FORMAT_SHORT_STR);}	  

		/*
		 * phase 1 
		 */
		{
			if (opt->alg_phase1 == PRIMAL_BENSON)
			{
				if (opt->message_level >= 3) printf("\n********* PHASE 1 -- PRIMAL ALGORITHM *********\n");
				phase1_primal(sol, vlp, opt);
			}
			else
			{
				assert(opt->alg_phase1 == DUAL_BENSON);
				if (opt->message_level >= 3) printf("\n********* PHASE 1 -- DUAL ALGORITHM *********\n");
				phase1_dual(sol,vlp,opt);
			}
		}
	}

	/*
	 * phase 2
	 */
	if (opt->alg_phase2 == PRIMAL_BENSON)
	{
		if (opt->message_level >= 3) printf("\n********* PHASE 2 -- PRIMAL ALGORITHM *********\n");
		phase2_primal(sol, vlp, opt);
	}
	else
	{
		if (opt->message_level >= 3) printf("\n********* PHASE 2 -- DUAL ALGORITHM *********\n");
		phase2_dual(sol, vlp, opt);
	}

	if (sol->status == VLP_INFEASIBLE)
	{
		printf("VLP is infeasible\n");
		lp_free(0);
		vlp_free(vlp);
		sol_free(sol);
		return 1;
	}
	if (sol->status == VLP_UNBOUNDED)
	{
		if (opt->bounded==1)
			printf("VLP is not bounded, re-run without option -b\n");
		else
			printf("LP in phase 2 is not bounded, probably by inaccuracy in phase 1 \n");
		lp_free(0);
		vlp_free(vlp);
		sol_free(sol);
		return 1;
	}
	
	double elapsedTime = (t_end.tv_sec - t_start.tv_sec) * 1000.0;	// sec to ms
	elapsedTime += (t_end.tv_usec - t_start.tv_usec) / 1000.0;		// us to ms

	// write log file
	{
		FILE *log_fp;
		char filename[strlen(opt->filename)+4+1];
		strcpy(filename,opt->filename);
		strcat(filename, ".log");
		log_fp = fopen (filename, "w+");
		if (log_fp == NULL)	
		{
			printf("unable to open file %s", filename);
			exit(1);
		}
	
		fprintf(log_fp, "BENSOLVE: VLP solver, %s\n", THISVERSION);
#ifdef LOG_HOST_NAME
		{
			char hostname[64]={0};
			gethostname(hostname, 63);
			fprintf(log_fp, "  host name:         %s\n", hostname);
		}
#endif		
		fprintf(log_fp, "Problem parameters\n");
		fprintf(log_fp, "  problem file:      %s\n", filename);
		fprintf(log_fp, "  problem rows:      %7d\n", vlp->m);
		fprintf(log_fp, "  problem columns:   %7d\n", vlp->n);
		fprintf(log_fp, "  matrix non-zeros:  %7lu\n", vlp->nz);
		fprintf(log_fp, "  primal generators: %7d\n", sol->o);
		fprintf(log_fp, "  dual generators:   %7d\n", sol->p);
		fprintf(log_fp, "Options\n");
		fprintf(log_fp, "  bounded:            %s\n", opt->bounded ? "yes (run phase 2 only)" : "no (run phases 0 to 2)");
		fprintf(log_fp, "  solution:           %s\n", opt->solution == PRE_IMG_OFF ? "off (no solution output)":"on (solutions (pre-image) written to files)"); 
		fprintf(log_fp, "  format:             %s\n", opt->format == FORMAT_AUTO ? "auto" : opt->format == FORMAT_LONG ? "long": "short");
		fprintf(log_fp, "  lp_method_phase0:   %s\n", opt->lp_method_phase0 == PRIMAL_SIMPLEX ? "primal_simplex" : opt->lp_method_phase0 == DUAL_SIMPLEX ? "dual_simplex" : "dual_primal_simplex (dual simplex, if not succesful, primal simplex)");
		fprintf(log_fp, "  lp_method_phase1:   %s\n", opt->lp_method_phase1 == PRIMAL_SIMPLEX ? "primal_simplex" : opt->lp_method_phase1 == DUAL_SIMPLEX ? "dual_simplex" : opt->lp_method_phase1 == DUAL_PRIMAL_SIMPLEX ? "dual_primal_simplex (dual simplex, if not succesful, primal simplex)" : "auto");
		fprintf(log_fp, "  lp_method_phase2:   %s\n", opt->lp_method_phase2 == PRIMAL_SIMPLEX ? "primal_simplex" : opt->lp_method_phase2 == DUAL_SIMPLEX ? "dual_simplex" : opt->lp_method_phase2 == DUAL_PRIMAL_SIMPLEX ? "dual_primal_simplex (dual simplex, if not succesful, primal simplex)" : "auto");
		fprintf(log_fp, "  message_level:      %d\n", opt->message_level);
		fprintf(log_fp, "  lp_message_level:   %d\n", opt->lp_message_level);
		fprintf(log_fp, "  alg_phase1:         %s\n", opt->alg_phase1 == PRIMAL_BENSON ? "primal" : "dual");
		fprintf(log_fp, "  alg_phase2:         %s\n", opt->alg_phase2 == PRIMAL_BENSON ? "primal" : "dual");
		fprintf(log_fp, "  eps_benson_phase1:  %g\n", opt->eps_benson_phase1);
		fprintf(log_fp, "  eps_benson_phase2:  %g\n", opt->eps_benson_phase2);
		fprintf(log_fp, "  eps_phase0:         %g\n", opt->eps_phase0);
		fprintf(log_fp, "  eps_phase1:         %g\n", opt->eps_phase1);
		fprintf(log_fp, "Computational results\n");
		fprintf(log_fp, "  CPU time (ms):      %g\n", elapsedTime);
		fprintf(log_fp, "  # LPs:              %d\n", lp_get_num(0));
		fprintf(log_fp, "Solution properties\n");
		fprintf(log_fp, "  # primal solution points:     %7zu\n", sol->pp);
		fprintf(log_fp, "  # primal solution directions: %7zu\n", sol->pp_dir);
		fprintf(log_fp, "  # dual solution points:       %7zu\n", sol->dd);
		fprintf(log_fp, "  # dual solution directions:   %7zu\n", sol->dd_dir);
		fclose(log_fp);
	}

	if (opt->message_level >= 1)
	{
		printf("CPU time            : %.4g %s.\n", elapsedTime >= 1000 ? elapsedTime / 1000 : elapsedTime, elapsedTime >= 1000 ? "s" : "ms");
		printf("Number of LPs solved: %d.\n", lp_get_num(0));
		printf("\n");
	}

	lp_free(0);
	vlp_free(vlp);
	sol_free(sol);
}
