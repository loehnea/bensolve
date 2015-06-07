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

#include <setjmp.h>
#include <assert.h> // assert,
#include <stdlib.h> // strtod,
#include <ctype.h>  // isspace, iscntrl, isdigit,
#include <limits.h> // INT_MIN, INT_MAX
#include <float.h>  // DBL_MIN, DBL_MAX
#include <string.h> // strcmp
#include <stdio.h>  // fopen, fgetc, ferror, fclose
#include <math.h>

#include "bslv_vlp.h"
#include "bslv_algs.h"		/* cone_vertenum */


static void error(csatype *csa, char const *msg)
{
	csa->error = 1;
	for(int i = 0; i<sizeof(csa->msg)-1; i++)
	{
		csa->msg[i] = msg[i];
		if (csa->msg[i] == '\0')
			break;
	}
	longjmp(csa->jump, 1);
}

static void warning(csatype *csa, char const *msg)
{
	csa->warning = 1;
	for(int i = 0; i < sizeof(csa->msg)-1; i++)
	{
		csa->msg[i] = msg[i];
		if (csa->msg[i] == '\0')
			break;
	}
	if (STOP_AT_WARNING)
		longjmp(csa->jump, 1);
}

static void read_char(csatype *csa)
{	// read character from input text file
	int c;
	if (csa->c == '\n')
		csa->count++;
	c = fgetc(csa->fp);
	if (c < 0)
	{
		if (ferror(csa->fp))
			error(csa,"reading error");
		else if (csa->c == '\n')
			error(csa,"unexpected end of file");
		else
		{
			// missing final end of line
			c = '\n';
		}
	}
	else if (c == '\n')
		;
	else if (isspace(c))
		c = ' ';
	else if (iscntrl(c))
		error(csa,"invalid control character");
	csa->c = c;
}

static void read_designator(csatype *csa)
{   // read one-character line designator
	assert(csa->c == '\n');
	read_char(csa);
	for (;;)
	{
		while (csa->c == ' ') // skip preceding white-space characters
			read_char(csa);
		if (csa->c == '\n')
		{
			if (!csa->empty)
			{
				warning(csa,"empty line ignored");
				csa->empty = 1;
			}
			read_char(csa);
		}
		else if (csa->c == 'c')
		{ // skip comment line
			while (csa->c != '\n')
				read_char(csa);
			read_char(csa);
		}
		else
		{   // candidate for a line designator
			csa->field[0] = (char)csa->c, csa->field[1] = '\0';
			// check that it is followed by a white-space character
			read_char(csa);
			if (!(csa->c == ' ' || csa->c == '\n'))
				error(csa,"line designator missing or invalid");
			break;
		}
	}
}

static void read_field(csatype *csa)
{
	// read data field
	int len = 0;
	// skip preceding white-space characters
	while (csa->c == ' ')
		read_char(csa);
	// scan data field 
	if (csa->c == '\n')
		error(csa,"unexpected end of line");
	while (!(csa->c == ' ' || csa->c == '\n'))
	{
		if (len == sizeof(csa->field)-1)
			error(csa,"data field to long");
		csa->field[len++] = (char)csa->c;
		read_char(csa);
	}
	csa->field[len] = '\0';
}

static int is_field_read(csatype *csa)
{
	// read data field
	int len = 0;
	// skip preceding white-space characters
	while (csa->c == ' ')
		read_char(csa);
	// scan data field 
	if (csa->c == '\n')
		return 0;
	while (!(csa->c == ' ' || csa->c == '\n'))
	{
		if (len == sizeof(csa->field)-1)
			error(csa,"data field to long");
		csa->field[len++] = (char)csa->c;
		read_char(csa);
	}
	csa->field[len] = '\0';
	return 1;
}

static void end_of_line(csatype *csa)
{
	// skip white-space characters until end of line 
	while (csa->c == ' ')
		read_char(csa);
	if (csa->c != '\n')
		error(csa,"too many data fields specified");
}

static int getint(csatype *csa)
{
	int d, k, s, val = 0;
	// scan optional sign
	if (csa->field[0] == '+')
		s = +1, k = 1;
	else if (csa->field[0] == '-')
		s = -1, k = 1;
	else
		s = +1, k = 0;
	// check for the first digit
	if (!isdigit((unsigned char)csa->field[k]))
		csa->error = 1;
	// scan digits
	while (isdigit((unsigned char)csa->field[k]))
	{
		d = csa->field[k++] - '0';
		if (s > 0)
		{
			if (val > INT_MAX / 10)
				csa->error = 1;
			val *= 10;
			if (val > INT_MAX - d)
				csa->error = 1;
			val += d;
		}
		else
		{
			if (val < INT_MIN / 10)
				csa->error = 1;
			val *= 10;
			if (val < INT_MIN + d)
				csa->error = 1;
			val -= d;
		}
	}
	// check for terminator
	if (csa->field[k] != '\0')
		csa->error = 1;
	// conversion has been done
	if (csa->error == 1)
		return 0;
	return val;
}

static double getnum(csatype *csa)
{
	int k;
	double val;
	// scan optional sign
	k = (csa->field[0] == '+' || csa->field[0] == '-' ? 1 : 0);
	// check for decimal point
	if (csa->field[k] == '.')
	{
		k++;
		// a digit should follow it
		if (!isdigit((unsigned char)csa->field[k])) 
			csa->error = 1;
		k++;
		goto frac;
	}
	// integer part should start with a digit
	if (!isdigit((unsigned char)csa->field[k])) 
		csa->error = 1;
	// scan integer part
	while (isdigit((unsigned char)csa->field[k])) k++;
	// check for decimal point
	if (csa->field[k] == '.') k++;
	frac: // scan optional fraction part
	while (isdigit((unsigned char)csa->field[k])) k++;
	// check for decimal exponent
	if (csa->field[k] == 'E' || csa->field[k] == 'e')
	{
		k++;
		// scan optional sign
		if (csa->field[k] == '+' || csa->field[k] == '-') k++;
		// a digit should follow E, E+ or E-
		if (!isdigit((unsigned char)csa->field[k]))
			csa->error = 1;
	}
	// scan optional exponent part
	while (isdigit((unsigned char)csa->field[k])) k++;
	// check for terminator
	if (csa->field[k] != '\0') 
		csa->error = 1;
	// perform conversion
	{
		char *endptr;
		val = strtod(csa->field, &endptr);
		if (*endptr != '\0')
			csa->error = 1;
	}
	// check for overflow
	if (!(-DBL_MAX <= val && val <= +DBL_MAX)) 
		csa->error = 1;
	// check for underflow
	if (-DBL_MIN < val && val < +DBL_MIN) val = 0.0;
	// conversion has been done
	if (csa->error == 1)
		return 0;
	return val;
}

int vlp_init(char const *filename, vlptype *vlp, const opttype *opt)
{
	// initialize struct vlp
	{
		vlp->A_ext = NULL;
		vlp->rows = NULL;
		vlp->cols= NULL;
		vlp->optdir = 0;
		vlp->cone_gen = 0;
		vlp->gen = NULL;
		vlp->c = NULL;
		vlp-> nz = 0;
		vlp-> nzobj = 0;
		vlp-> n = 0;
		vlp-> m = 0;
		vlp-> q = 0;
		vlp-> n_gen = 0;
	}

	csatype _csa, *csa = &_csa;
	if (setjmp(csa->jump))
		goto done;
	csa->fname = filename;
	csa->fp = NULL;
	csa->count = 0;
	csa->c = '\n';
	csa->field[0] = '\0';
	csa->empty = csa->nonint = 0;
	csa->fp = fopen(csa->fname, "r");
	csa->error = 0;
	csa->warning = 0;
	if (csa->fp == NULL)
		error(csa, "file not found or unable to open");
	// read problem line
	read_designator(csa);
	if (strcmp(csa->field, "p") != 0)
		error(csa,"problem line missing or invalid");
	read_field(csa);
	if (strcmp(csa->field, PROBLEM_DESIGNATOR) != 0)
		error(csa, "wrong problem designator");
	read_field(csa);
	if (strcmp(csa->field, "min") == 0)
		vlp->optdir=1;
	else if (strcmp(csa->field, "max") == 0)
		vlp->optdir=-1;
	else
		error(csa,"objective sense missing or invalid");
	read_field(csa);
	vlp->m=getint(csa);
	if (csa->error || vlp->m < 0)
		error(csa,"number of rows missing or invalid");
	read_field(csa);
	vlp->n=getint(csa);
	if (csa->error || vlp->n < 0)
		error(csa,"number of columns missing or invalid");
	read_field(csa);
	vlp->nz=getint(csa);
	if (csa->error || vlp->nz < 0)
		error(csa,"number of nonzeros missing or invalid");
	read_field(csa);
	vlp->q=getint(csa);
	if (csa->error || vlp->q < 1)
		error(csa,"number of objectives missing or invalid");
	read_field(csa);
	vlp->nzobj=getint(csa);
	if (csa->error || vlp->nzobj < 0)
		error(csa,"number of objective matrix nonzeros missing or invalid");
	if (is_field_read(csa))
	{
		if (strcmp(csa->field, "cone") == 0)
			vlp->cone_gen=CONE;
		else if (strcmp(csa->field, "dualcone") == 0){
			vlp->cone_gen=DUALCONE;
		}
		else
			error(csa,"type of cone generators missing or invalid");
	}
	else 
		vlp->cone_gen=DEFAULT;

	int nzgen = 0;
	if (vlp->cone_gen==CONE || vlp->cone_gen==DUALCONE)
	{
		read_field(csa);
		vlp->n_gen=getint(csa);
		if (csa->error || vlp->n_gen < 0)
			error(csa,"number of cone generating vectors missing or invalid");
		read_field(csa);
		nzgen=getint(csa);
		if (csa->error || nzgen < 0)
			error(csa,"number of cone generator non-zeros missing or invalid");	
	}
	else
	{
		vlp->n_gen = 0; // number of generators
	}
	end_of_line(csa);

	// read other lines

	// extended coefficient matrix A_ext = (A,0; -P, I)
	vlp->A_ext = list2d_calloc(vlp->nz + vlp->nzobj + vlp->q); 
	for (int i = 0; i < vlp->q; i++) // set -I
	{
		vlp->A_ext->idx1[vlp->nz + vlp->nzobj + i] = vlp->m + i + 1;
		vlp->A_ext->idx2[vlp->nz + vlp->nzobj + i] = vlp->n + i + 1;
		vlp->A_ext->data[vlp->nz + vlp->nzobj + i] = 1.0;
	}

	// row descriptors (all entries are stored because number is unknown)
	vlp->rows = boundlist_calloc(vlp->m, 'x');	// initialize with invalid type
	boundlist_init_idx(vlp->rows, 1);			// initialize indices: 1,...,rows->size

	// columns descriptors (all entries are stored because number is unknown)
	vlp->cols = boundlist_calloc(vlp->n, 'x');
	boundlist_init_idx(vlp->cols, 1); // initialize indices: 1,...,cols->size

	// geometric duality parameter vector c (all entries are stored)
	// must be inizialized with zero
	vlp->c = (double *) calloc (vlp->q, sizeof(double)); 

	// generating vectors of orderung cone (vlp->cone_gen == CONE) or
	// generating vectors of the dual of the ordering cone (vlp->cone_gen == DUALCONE)
	if (vlp->cone_gen == CONE || vlp->cone_gen == DUALCONE)
		vlp->gen = calloc(vlp->q*vlp->n_gen, sizeof(double));
	else
		vlp->gen = NULL;

	int ridx = 0, cidx = 0, na = 0, no = 0, nk = 0;
	char type;
	for(;;)
	{
		read_designator(csa);
		// coefficient descriptor
		if (strcmp(csa->field, "a") == 0)
		{
			if (na == vlp->nz)
				error(csa,"too many constraint coefficient descriptors");
			read_field(csa);
			ridx = getint(csa);
			if (csa->error)
				error(csa,"constraint coefficient row number missing or invalid");
			if (1 > ridx || ridx > vlp->m)
				error(csa,"constraint coefficient row number out of range");
			read_field(csa);
			cidx = getint(csa);
			if (csa->error)
				error(csa,"constraint coefficient column number  missing or invalid");
			if (1 > cidx || cidx > vlp->n)
				error(csa,"constraint coefficient column number out of range");
			vlp->A_ext->idx1[na] = ridx;
			vlp->A_ext->idx2[na] = cidx;
			read_field(csa);
			vlp->A_ext->data[na] = getnum(csa);
			if (csa->error)
				error(csa,"constraint coefficient missing or invalid");
			na++;
		}
		// objective matrix descriptor
		else if (strcmp(csa->field, "o") == 0)
		{
			if (no == vlp->nzobj)
				error(csa,"too many objective coefficient descriptors");
			read_field(csa);
			ridx = getint(csa);
			if (csa->error)
				error(csa,"objective coefficient row number missing or invalid");
			if (1 > ridx || ridx > vlp->q)
				error(csa,"objective coefficient row number out of range");
			read_field(csa);
			cidx = getint(csa);
			if (csa->error)
				error(csa,"objective coefficient column number missing or invalid");
			if (1 > cidx || cidx > vlp->n)
				error(csa,"objective coefficient column number out of range");
			vlp->A_ext->idx1[vlp->nz + no] = vlp->m + ridx;
			vlp->A_ext->idx2[vlp->nz + no] = cidx;
			read_field(csa);
			vlp->A_ext->data[vlp->nz + no] = -getnum(csa);
			if (csa->error)
				error(csa,"objective coefficient missing or invalid");
			no++;
		}
		// cone generator matrix descriptor
		else if (strcmp(csa->field, "k") == 0)
		{
			if(vlp->cone_gen != CONE && vlp->cone_gen != DUALCONE)
				error(csa,"invalid designator k");
			read_field(csa);
			ridx = getint(csa);
			if (csa->error)
				error(csa,"cone generator coefficient row number missing or invalid");
			if (1 > ridx || ridx > vlp->n_gen)
				error(csa,"cone generator coefficient row number out of range");
			read_field(csa);
			cidx = getint(csa);
			if (csa->error)
				error(csa,"cone generator coefficient column number missing or invalid");
			if (0 > cidx || cidx > vlp->n_gen)
				error(csa,"cone generator coefficient column number out of range");
			read_field(csa);
			if (cidx==0) // store geometric duality parameter vector c
				vlp->c[ridx-1] = getnum(csa);
			else // store cone data
			{
				if (nk == nzgen)
					error(csa,"too many cone generator coefficient descriptors");
				vlp->gen[vlp->n_gen*(ridx-1)+(cidx-1)] = getnum(csa);
				nk++;
			}
			if (csa->error)
				error(csa,"cone generator coefficient missing or invalid");
		}

		// row descriptor
		else if (strcmp(csa->field, "i") == 0)
		{
			read_field(csa);
			ridx = getint(csa);
			if (csa->error)
				error(csa,"row number missing or invalid");
			if (1 > ridx || ridx > vlp->m)
				error(csa,"row number out of range");
			if (vlp->rows->type[ridx-1] != 'x')
				error(csa,"duplicate row descriptor");
			read_field(csa);
			if (strcmp(csa->field, "f") == 0) type = 'f';
			else if (strcmp(csa->field, "l") == 0)	type = 'l';
			else if (strcmp(csa->field, "u") == 0)	type = 'u';
			else if (strcmp(csa->field, "d") == 0)	type = 'd';
			else if (strcmp(csa->field, "s") == 0)	type = 's';
			else
				error(csa,"row type missing or invalid");
			if (type == 'l' || type == 'd' || type == 's')
			{
				read_field(csa);
				vlp->rows->lb[ridx-1] = getnum(csa);
				if (csa->error)
					error(csa,"row lower bound missing or invalid");
			}
			if (type == 'u' || type == 'd')
			{
				read_field(csa);
				vlp->rows->ub[ridx-1] = getnum(csa);
				if (csa->error)
					error(csa,"row upper bound missing or invalid");
			}
			vlp->rows->type[ridx-1] = type;
		}
		// column descriptor
		else if (strcmp(csa->field, "j") == 0)
		{
			read_field(csa);
			cidx = getint(csa);
			if (csa->error)
				error(csa,"column number missing or invalid");
			if (1 > cidx || cidx > vlp->n)
				error(csa,"column descriptor out of range");
			if (vlp->cols->type[cidx-1] != 'x')
				error(csa,"duplicate column descriptor");
			read_field(csa);
			if      (strcmp(csa->field, "f") == 0) type = 'f';
			else if (strcmp(csa->field, "l") == 0) type = 'l';
			else if (strcmp(csa->field, "u") == 0) type = 'u';
			else if (strcmp(csa->field, "d") == 0) type = 'd';
			else if (strcmp(csa->field, "s") == 0) type = 's';
			else
				error(csa,"column type missing or invalid");
			if (type == 'l' || type == 'd' || type == 's')
			{
				read_field(csa);
				vlp->cols->lb[cidx-1] = getnum(csa);
				if (csa->error)
					error(csa,"column lower bound missing or invalid");
			}
			if (type == 'u' || type == 'd')
			{
				read_field(csa);
				vlp->cols->ub[cidx-1] = getnum(csa);
				if (csa->error)
					error(csa,"column upper bound missing or invalid");
			}
			vlp->cols->type[cidx-1] = type;
		}
		else if (strcmp(csa->field, "e") == 0)
			break;
		else
			error(csa,"line designator missing or invalid");
		end_of_line(csa);
	}

	// replace undefined rows by default bounds
	for (int i = 0; i < vlp->m; i++)
		if (vlp->rows->type[i] == 'x')
			vlp->rows->type[i] = 'f';

	// replace undefined columns by default bounds
	for (int j = 0; j < vlp->n; j++)
		if (vlp->cols->type[j] == 'x')
			vlp->cols->type[j] = 's';
	done:
	if (csa->warning)
	{
		if (STOP_AT_WARNING)
			printf("Warning occurred while reading %s: stopped at line %d: %s\n", filename, csa->count, csa->msg);
		else
			printf("Warning occurred while reading %s: %s\n", filename, csa->msg);
	}
	if (csa->error)
		printf("Error while reading %s: line %d: %s\n", filename, csa->count, csa->msg);
	if (csa->fp != NULL) fclose(csa->fp);
	if (csa->error) return 1;
	return 0;
}

void vlp_free(vlptype *vlp)
{
	if (vlp->A_ext != NULL)	{list2d_free(vlp->A_ext); vlp->A_ext=NULL;}
	if (vlp->rows != NULL) {boundlist_free(vlp->rows); vlp->rows=NULL;}
	if (vlp->cols != NULL) {boundlist_free(vlp->cols); vlp->cols=NULL;}
	if (vlp->gen != NULL) {free(vlp->gen); vlp->gen=NULL;}
	if (vlp->c != NULL) {free(vlp->c); vlp->c=NULL;}
}

void sol_init(soltype *sol, const vlptype *vlp, const opttype *opt)
{
	// initialize sol to zero
	{
		sol->m = 0;
		sol->n = 0;
		sol->q = 0;
		sol->o = 0;
		sol->p = 0;
		sol->r = 0;
		sol->h = 0;
		sol->eta = NULL;
		sol->Y = NULL;
		sol->Z = NULL;
		sol->c = NULL;
		sol->R = NULL;
		sol->H = NULL;
		sol->status = 0;
		sol->c_dir = 0;
		sol->pp=0;
		sol->dd=0;
		sol->pp_dir=0;
		sol->dd_dir=0;
	}

	sol->m = vlp->m;
	sol->n = vlp->n;
	sol->q = vlp->q;
	sol->eta = calloc(sol->q, sizeof(double));

	if (vlp->cone_gen == CONE) // generators of C are given
	{
		int flag=cone_vertenum(&sol->Y,&sol->o,&sol->Z,&sol->p,vlp->gen,vlp->n_gen,vlp->q,opt,CONE_OUT_ON,SWAP);
		if (flag==EXIT_FAILURE) // cone is not pointed
		{
			sol->status=VLP_INPUTERROR;
			printf("Input error: Ordering cone has empty interior (1)\n");
			return;
		}
		else if (sol->p < vlp->q || sol->o < vlp->q) // test does not cover all cases which may orrur
		{
			sol->status=VLP_INPUTERROR;
			printf("Input error: Ordering cone is not pointed (2)\n");
			return;
		}
	}
	else if(vlp->cone_gen == DUALCONE) // generators of C^* are given
	{
		int flag=cone_vertenum(&sol->Z,&sol->p,&sol->Y,&sol->o,vlp->gen,vlp->n_gen,vlp->q,opt,CONE_OUT_ON,NO_SWAP);
		if (flag==EXIT_FAILURE) 
		{
			sol->status=VLP_INPUTERROR;
			printf("Input error: Ordering cone is not pointed (1)\n");
			return;
		}
		else if (sol->p < vlp->q || sol->o < vlp->q) // test does not cover all cases which may orrur
		{
			sol->status=VLP_INPUTERROR;
			printf("Input error: Ordering cone has empty interior (2)\n");
			return;
		}
	}
	else  // standard cone R^q_+
	{
		assert(vlp->cone_gen == DEFAULT);
		sol->Y=calloc(vlp->q*vlp->q,sizeof(double));
		sol->Z=calloc(vlp->q*vlp->q,sizeof(double));
		for (size_t k=0; k<vlp->q; k++){
			sol->Y[k*vlp->q+k]=1.0;
			sol->Z[k*vlp->q+k]=1.0;
		}
		sol->p=vlp->q;
		sol->o=vlp->q;
	}
	
	// copy c and scale c such that c_q == 1
	{
		sol->c=calloc(vlp->q,sizeof(double));
		if (vlp->cone_gen == DEFAULT)
		{
			for (size_t j=0; j<vlp->q; j++)
				sol->c[j] = 1.0;
			sol->c_dir=C_DIR_POS;
		}
		else
		{
			assert(vlp->cone_gen == CONE || vlp->cone_gen == DUALCONE);
			// scale columns of sol->Y (2-norm == 1)
			{
				for (size_t k=0;k<sol->o;k++)
				{
					double tmp=0;
					for (size_t j=0;j<vlp->q;j++)
						tmp+=sol->Y[k+j*sol->o]*sol->Y[k+j*sol->o];
					for (size_t j=0;j<vlp->q;j++)
						sol->Y[k+j*sol->o]/=sqrt(tmp);
				}
			}

			if (ABS(vlp->c[vlp->q-1]) > EPS_C)
			{
				for (size_t i=0; i<vlp->q; i++)
					sol->c[i] = vlp->c[i] / ABS(vlp->c[vlp->q-1]);
				sol->c_dir=(vlp->c[vlp->q-1]>0)?C_DIR_POS:C_DIR_NEG;
			}
			else
			{
				double tmp1[vlp->q];
				double tmp2[vlp->q];
				double max=0, min=0;
				size_t k1=0,k2=0;
				for (size_t j=0; j<vlp->q; j++)
				{
					tmp1[j] = 0;
					tmp2[j] = 0;
				}
				for (size_t i=0; i<sol->o; i++)
				{
					if (sol->Y[(vlp->q-1)*sol->o+i]>0) // collect vectors with positive last component
					{
						if (sol->Y[(vlp->q-1)*sol->o+i]>max)
							max=sol->Y[(vlp->q-1)*sol->o+i];
						for (size_t j=0; j<vlp->q; j++)
							tmp1[j] += sol->Y[j*sol->o+i];
						k1++;
					}
					else // and vectors with non-positive last component
					{
						if (sol->Y[(vlp->q-1)*sol->o+i]<min)
							min=sol->Y[(vlp->q-1)*sol->o+i];
						for (size_t j=0; j<vlp->q; j++)
							tmp2[j] += sol->Y[j*sol->o+i];
						k2++;
					}
				}
				if (k1 == 0 && min<EPS_C) // no vector with positive component and at least one with negative
				{
					sol->c_dir=C_DIR_NEG;
					for (size_t i=0; i<vlp->q; i++) 
						sol->c[i] = tmp2[i] / ABS(tmp2[vlp->q-1]);
				}
				else if (k2 == 0 && max>EPS_C) // no vector with non-positive component and at least one with positive
				{
					sol->c_dir=C_DIR_POS;
					for (size_t i=0; i<vlp->q; i++) 
						sol->c[i] = tmp1[i] / ABS(tmp1[vlp->q-1]);
				}
				else if (min<-EPS_C || max>EPS_C)
				{
					double lambda;
					if (-min>max)
					{
						sol->c_dir=C_DIR_NEG;
						lambda=0.2*(-min/(max-min));
					}
					else
					{
						sol->c_dir=C_DIR_POS;
						lambda=0.8-0.2*min/(max-min);
					}
					for (size_t i=0; i<vlp->q; i++)
						sol->c[i] = lambda*tmp1[i]/k1 +(1-lambda)*tmp2[i]/k2;
					for (size_t i=0; i<vlp->q; i++)
						sol->c[i] = sol->c[i] / ABS(sol->c[vlp->q-1]);
				}
				else
				{
					printf("Input error: ordering cone is not solid (3)\n");
					sol->status=VLP_INPUTERROR;
					return;
				}
				if (opt->message_level >= 1) printf("Warning: geometric duality parameter vector c was generated\n");
			}
		}
	}

	// scale columns of sol->Z (Z' * c == (1,...,1)')
	{
		double tmp;
		for (size_t k=0; k<sol->p; k++)
		{
			tmp=0;
			for (size_t j=0; j<vlp->q; j++)
				tmp+=sol->Z[k+j*sol->p]*sol->c[j];
			if(tmp < 1e-8)
			{
				printf("Input error: c does not belong to interior of ordering cone\n");
				sol->status=VLP_INPUTERROR;
				return;
			}	
			for (size_t j=0; j<vlp->q; j++)
				sol->Z[k+j*sol->p]/=tmp;
		}
	} 

	// further tests whether C is pointed and solid
	if (vlp->cone_gen == CONE || vlp->cone_gen == DUALCONE)
	{
		double *sum_Y=calloc(sol->q,sizeof(double));
		double *sum_Z=calloc(sol->q,sizeof(double));
		for(size_t j=0; j<sol->q; j++)
			for(size_t k=0; k<sol->o; k++)
				sum_Y[j]+=sol->Y[sol->o*j+k];
		for(size_t j=0; j<sol->q; j++)
			for(size_t k=0; k<sol->p; k++)
				sum_Z[j]+=sol->Z[sol->p*j+k];
		for(size_t k=0; k<sol->p; k++)
		{
			double tmp = 0;
			for(size_t j=0; j<sol->q; j++)
				tmp+=sol->Z[sol->p*j+k]*sum_Y[j];
			if(tmp < 1e-8)
			{
				printf("Input error: ordering cone is not solid (4)\n");
				sol->status=VLP_INPUTERROR;
				return;
			}
		}
		for(size_t k=0; k<sol->o; k++)
		{
			double tmp = 0;
			for(size_t j=0; j<sol->q; j++)
				tmp+=sol->Y[sol->o*j+k]*sum_Z[j];
			if(tmp < 1e-8)
			{
				printf("Input error: ordering cone is not pointed (4)\n");
				sol->status=VLP_INPUTERROR;
				return;
			}
		}
		free(sum_Y);
		free(sum_Z);
	}
	
	// write c to file and stdout
	{
		double c[vlp->q];
		for(size_t k=0; k<vlp->q; k++)
			c[k]=sol->c[k];
		char filename[strlen(opt->filename)+6+1];
		strcpy(filename,opt->filename);
		strcat(filename, "_c.sol");
		matrix_fprint (c, 1, vlp->q, 1, filename, opt->format==FORMAT_SHORT?FORMAT_SHORT_STR:FORMAT_LONG_STR);
		if (opt->message_level >= 2) { printf("Duality parameter vector c = \n  "); matrix_print(sol->c, 1, vlp->q, opt->format==FORMAT_LONG?FORMAT_LONG_STR:FORMAT_SHORT_STR);}
	}

	// invert C and c in case of c_q<0 in order to obtain standard problem of type c_q>0
	if (sol->c_dir == C_DIR_NEG)
	{
		for (size_t k=0; k<sol->o*vlp->q; k++)
			sol->Y[k]=-sol->Y[k];
		for (size_t k=0; k<sol->p*vlp->q; k++)
			sol->Z[k]=-sol->Z[k];
		for (size_t k=0; k<vlp->q; k++)
			sol->c[k]=-sol->c[k];
	}

	// invert P in cases min/c_q<0 or max/c_q>0 in order to get standard problem of type min/c_q>0
	if ((sol->c_dir == C_DIR_NEG && vlp->optdir == 1) || (sol->c_dir == C_DIR_POS && vlp->optdir == -1))
	{
		for (lp_idx i = 0; i<vlp->nzobj; i++)
			vlp->A_ext->data[vlp->nz + i] = -vlp->A_ext->data[vlp->nz + i];
	}

	sol->status=VLP_NOSTATUS;
}

void sol_free(soltype *sol)
{
	if (sol->Z != NULL) {free(sol->Z); sol->Z=NULL;}
	if (sol->Y != NULL) {free(sol->Y); sol->Y=NULL;}
	if (sol->c != NULL) {free(sol->c); sol->c=NULL;}
	if (sol->R != NULL) {free(sol->R); sol->R=NULL;}
	if (sol->H != NULL) {free(sol->H); sol->H=NULL;}
	if (sol->eta != NULL) {free(sol->eta); sol->eta=NULL;}
}

void set_default_opt(opttype *opt)
{
	opt->bounded = 0;
	opt->plot = 0;
	opt->filename[0]='\0';
	opt->solution = PRE_IMG_OFF;
	opt->format = FORMAT_AUTO;
	opt->lp_method_phase0 = PRIMAL_SIMPLEX;
	opt->lp_method_phase1 = LP_METHOD_AUTO;
	opt->lp_method_phase2 = LP_METHOD_AUTO;
	opt->message_level = DEFAULT_MESSAGE_LEVEL;
	opt->lp_message_level = DEFAULT_LP_MESSAGE_LEVEL;
	opt->alg_phase1 = PRIMAL_BENSON;
	opt->alg_phase2 = PRIMAL_BENSON;
	opt->eps_phase0 = DEFAULT_EPS_PHASE0;
	opt->eps_phase1 = DEFAULT_EPS_PHASE1;	
	opt->eps_benson_phase1 = DEFAULT_EPS_BENSON_PHASE1;
	opt->eps_benson_phase2 = DEFAULT_EPS_BENSON_PHASE2;
}
