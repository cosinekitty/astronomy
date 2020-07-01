/*
    top2013.c  -  Don Cross  -  2020-06-30
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "codegen.h"
#include "top2013.h"


void TopInitModel(top_model_t *model)
{
    memset(model, 0, sizeof(top_model_t));
}


void TopFreeModel(top_model_t *model)
{
    int f, s;
    for (f=0; f < TOP_NCOORDS; ++f)
        for (s=0; s < model->formula[f].nseries_total; ++s)
            free(model->formula[f].series[s].terms);

    TopInitModel(model);
}


int TopLoadModel(top_model_t *model, const char *filename, int planet)
{
    int error = 1;
    FILE *infile = NULL;
    int nterms_remaining, lnum, nscanned, check_planet, check_var, tpower;
    top_formula_t *formula = NULL;
    top_series_t *series = NULL;
    top_term_t *term = NULL;
    char line[100];

    TopInitModel(model);
    infile = fopen(filename, "rt");
    if (infile == NULL)
        FAIL("TopLoadModel: cannot open file: %s\n", filename);

    lnum = 0;
    nterms_remaining = 0;
    check_planet = 0;
    check_var = 0;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;
        if (nterms_remaining == 0)
        {
            nscanned = sscanf(line,
                " TOP2013ELL PLANET %d VARIABLE %d T**%d %d term(s)",
                &check_planet, &check_var, &tpower, &nterms_remaining);

            if (nscanned != 4)
                FAIL("TopLoadModel(%s line %d): invalid data format.\n", filename, lnum);

            if (check_var < 1 || check_var > TOP_NCOORDS)
                FAIL("TopLoadModel(%s line %d): invalid variable number %d\n", filename, lnum, check_var);

            --check_var;    /* convert one-based to zero-based indexing */

            if (check_planet == planet)
            {
                if (series != NULL)
                {
                    /* make sure the previos series is complete */
                    if (series->nterms_calc != series->nterms_total)
                        FAIL("TopLoadModel(%s line %d): previous series has %d terms; expected %d\n", filename, lnum, series->nterms_calc, series->nterms_total);
                }

                /* Allocate room for the new terms. */
                formula = &model->formula[check_var];
                if (formula->nseries_total != tpower)
                    FAIL("TopLoadModel(%s line %d): expected tpower %d, found %d\n", filename, lnum, formula->nseries_total, tpower);
                if (formula->nseries_total >= TOP_MAX_SERIES)
                    FAIL("TopLoadModel(%s line %d): too many series.\n", filename, lnum);
                series = &formula->series[formula->nseries_total++];
                series->nterms_total = nterms_remaining;
                series->nterms_calc = 0;
                series->terms = calloc(sizeof(top_term_t), nterms_remaining);
                if (series->terms == NULL)
                    FAIL("TopLoadModel(%s line %d): out of memory\n", filename, lnum);
            }
        }
        else
        {
            --nterms_remaining;

            if (check_planet == planet)
            {
                if (series == NULL)
                    FAIL("TopLoadModel(%s line %d): series == NULL\n", filename, lnum);

                if (series->nterms_calc >= series->nterms_total)
                    FAIL("TopLoadModel(%s line %d): too many terms\n", filename, lnum);

                term = &series->terms[series->nterms_calc++];

                if (strlen(line) < 61)
                    FAIL("TopLoadModel(%s line %d) line is too short.\n", filename, lnum);

                /* patch in 'e' to make numbers in scientific notation. */
                if (line[31] != ' ' || line[57] != ' ')
                    FAIL("TopLoadModel(%s line %d): expected spaces between mantissas and exponents.\n", filename, lnum);

                line[31] = line[57] = 'e';
                nscanned = sscanf(line, "%lf %lf %lf %lf", &term->k, &term->c, &term->s, &term->p);
                if (nscanned == 3)
                    term->p = 0.0;
                else if (nscanned != 4)
                    FAIL("TopLoadModel(%s line %d): invalid term data format.\n", filename, lnum);
            }
        }
    }
    if (nterms_remaining != 0)
        FAIL("TopLoadModel(%s): missing %d terms at the end.\n", filename, nterms_remaining);

    error = 0;
fail:
    if (infile != NULL) fclose(infile);
    if (error) TopFreeModel(model);
    return error;
}
