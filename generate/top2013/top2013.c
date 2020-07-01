/*
    top2013.c  -  Don Cross  -  2020-06-30
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "codegen.h"
#include "top2013.h"

static int FormatTermLine(int lnum, char *line, size_t size, const top_term_t *term);


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


static int RoundingAdjustment(char original, char regen, int *diff)
{
    switch (original - regen)
    {
    case 0:
        *diff = 0;
        return 0;

    case -1:
    case +9:
        *diff = -1;
        return 0;

    case +1:
    case -9:
        *diff = +1;
        return 0;

    default:
        fprintf(stderr, "RoundingAdjustment: original=%c, regen=%c\n", original, regen);
        return 1;
    }
}


int TopLoadModel(top_model_t *model, const char *filename, int planet)
{
    int error = 1;
    FILE *infile = NULL;
    int nterms_remaining, lnum, nscanned, check_planet, check_var, tpower;
    int count = 0;      /* total number of terms we found for this planet */
    top_formula_t *formula = NULL;
    top_series_t *series = NULL;
    top_term_t *term = NULL;
    char line[100];
    char gline[100];

    TopInitModel(model);
    model->planet = planet;

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
                    /* make sure the previous series is complete */
                    if (series->nterms_calc != series->nterms_total)
                        FAIL("TopLoadModel(%s line %d): previous series has %d terms; expected %d\n", filename, lnum, series->nterms_calc, series->nterms_total);
                }

                /* Allocate room for the new terms. */
                formula = &model->formula[check_var];
                if (formula->nseries_total >= TOP_MAX_SERIES)
                    FAIL("TopLoadModel(%s line %d): too many series.\n", filename, lnum);
                series = &formula->series[formula->nseries_total++];
                formula->nseries_calc = formula->nseries_total;
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

                /*
                    Super ugly hack: the weird exponential notation used in TOP2013.dat
                    makes it really difficult to regenerate exactly in all cases.
                    So I generate each line back as text, and remember what I will have to do later
                    to reconstruct the input exactly.
                */
                CHECK(FormatTermLine(lnum, gline, sizeof(gline), term));
                CHECK(RoundingAdjustment(line[30], gline[30], &term->rc));
                CHECK(RoundingAdjustment(line[56], gline[56], &term->rs));

                /* Sanity check the rounding hack. */
                CHECK(FormatTermLine(lnum, gline, sizeof(gline), term));
                line[31] = line[57] = ' ';
                if (strcmp(line, gline))
                {
                    fprintf(stderr, "INPUT:%s", line);
                    fprintf(stderr, "REGEN:%s", gline);
                    FAIL("TopLoadModel(%s line %d): unable to reconstruct identical term line.\n", filename, lnum);
                }

                ++count;
            }
        }
    }

    if (nterms_remaining != 0)
        FAIL("TopLoadModel(%s): missing %d terms at the end.\n", filename, nterms_remaining);

    if (count == 0)
        FAIL("TopLoadModel(%s): could not find any terms for planet %d.\n", filename, planet);

    error = 0;
fail:
    if (infile != NULL) fclose(infile);
    if (error) TopFreeModel(model);
    return error;
}


static int AppendTrigCoeff(char *line, int lnum, double x, int rounding_adjust)
{
    int error = 1;
    int m, length, exponent;
    char polarity;
    char buffer[40];

    if (rounding_adjust < -1 || rounding_adjust > +1)
        FAIL("AppendTrigCoeff(%d): invalid rounding_adjust = %d\n", lnum, rounding_adjust);

    /*
        The data format for TOP2013 has a weird form of scientific notation.
        It always has a 0 to the left of the decimal point, which wastes a digit of precision.
        "-5.2026032025158849e+00"       <== C scientific notation
        "   -0.5202603202515885 +01"    <== output
    */

    snprintf(buffer, sizeof(buffer), "%23.16le", x);
    length = strlen(buffer);
    if (length != 23)
        FAIL("AppendTrigCoeff(%d): output string '%s' has incorrect length %d.\n", lnum, buffer, length);

    if (buffer[19] != 'e')
        FAIL("AppendTrigCoeff(%d): expected 'e' at index 19 in string '%s'\n", lnum, buffer);

    buffer[19] = '\0';      /* truncate the string at the 'e' */

    if (1 != sscanf(buffer+20, "%d", &exponent))
        FAIL("AppendTrigCoeff(%d): cannot scan exponent from '%s'\n", lnum, buffer+20);

    ++exponent;
    if (x == 0.0)
        exponent = 0;
    if (exponent >= 0)
    {
        polarity = '+';
    }
    else
    {
        polarity = '-';
        exponent = -exponent;
    }

    /* Copy digits and shift decimal point */
    buffer[22] = '\0';
    for(m = 17; m >= 0 && buffer[m] != '.'; --m)
        buffer[m+4] = buffer[m];

    if (m != 2 || buffer[2] != '.')
        FAIL("AppendTrigCoeff(%d): decimal point is in the wrong place: '%s'\n", lnum, buffer);

    buffer[6] = buffer[1];
    buffer[5] = '.';
    buffer[4] = '0';
    buffer[3] = buffer[0];
    buffer[2] = ' ';
    buffer[1] = ' ';
    buffer[0] = ' ';

    if (rounding_adjust != 0)
    {
        for (m=21; m >= 0; --m)
        {
            if (buffer[m] != '.')
            {
                if (buffer[m] < '0' || buffer[m] > '9')
                    FAIL("AppendTrigCoeff(%d): rounding failure\n", lnum);
                buffer[m] += rounding_adjust;
                if (buffer[m] < '0')
                    buffer[m] = '9';
                else if (buffer[m] > '9')
                    buffer[m] = '0';
                else
                    break;
            }
        }
    }

    sprintf(buffer+22, " %c%02d", polarity, exponent);
    length = strlen(buffer);
    if (length != 26)
        FAIL("AppendTrigCoeff(%d): generated incorrect length %d in string '%s' for x=%lg\n", lnum, length, buffer, x);

    strcat(line, buffer);
    error = 0;
fail:
    return error;
}


static int FormatTermLine(int lnum, char *line, size_t size, const top_term_t *term)
{
    int error = 1;
    int length;

    snprintf(line, size, "%9.0lf", term->k);
    CHECK(AppendTrigCoeff(line, lnum, term->c, term->rc));
    CHECK(AppendTrigCoeff(line, lnum, term->s, term->rs));
    length = strlen(line);
    if (61 != length)
        FAIL("FormatTermLine(%d): incorrect output line length = %d.\n", lnum, length);
    if (term->k != 0)
        snprintf(line + length, size - (size_t)length, " %11.6lf", term->p);
    strcat(line, "\n");
    error = 0;
fail:
    return error;
}


int TopWriteModel(const top_model_t *model, FILE *outfile)
{
    int error = 1;
    int f, s, t;
    int lnum = 0;
    char line[100];

    for (f=0; f < TOP_NCOORDS; ++f)
    {
        const top_formula_t *formula = &model->formula[f];
        for (s=0; s < formula->nseries_calc; ++s)
        {
            const top_series_t *series = &formula->series[s];

            ++lnum;
            if (0 > fprintf(outfile, " TOP2013ELL    PLANET %d    VARIABLE %d    T**%02d %7d term(s)\n", model->planet, f+1, s, series->nterms_calc))
                FAIL("TopWriteModel(%d): error writing header record to output stream.\n", lnum);

            for (t=0; t < series->nterms_calc; ++t)
            {
                const top_term_t *term = &series->terms[t];
                ++lnum;
                CHECK(FormatTermLine(lnum, line, sizeof(line), term));
                if (0 > fprintf(outfile, "%s", line))
                    FAIL("TopWriteModel(%d): error writing term record to output stream.\n", lnum);
            }
        }
    }
    error = 0;
fail:
    return error;
}


int TopSaveModel(const top_model_t *model, const char *filename)
{
    int error;
    FILE *outfile = fopen(filename, "wt");
    if (outfile == NULL)
        FAIL("TopSaveModel: Cannot open output file: %s\n", filename);
    CHECK(TopWriteModel(model, outfile));
fail:
    if (outfile) fclose(outfile);
    return error;
}
