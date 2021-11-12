/*
    MIT License

    Copyright (c) 2019-2021 Don Cross <cosinekitty@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "novas.h"
#include "codegen.h"
#include "top2013.h"
#include "ephfile.h"
#include "gravsim/pluto_gravsim.h"

#define CG_MAX_LINE_LENGTH  200
#define MAX_DATA_PER_LINE    20
#define IAU_DATA_PER_ROW     11
#define IAU_DATA_NUM_ROWS    77
#define ADDSOL_DATA_PER_ROW   8

typedef struct
{
    FILE *outfile;
    cg_language_t language;
    const char *datapath;
    const char *verb;
    const char *args;
    const char *inFileName;
    int lnum;
}
cg_context_t;

typedef int (* cg_directive_func) (cg_context_t *);

typedef struct
{
    const char *verb;
    cg_directive_func func;
}
cg_directive_entry;

static int ScanDirective(
    char *line,
    const char **verb_text,
    const char **args_text,
    const char **tail_text);

static int ProcessDirective(cg_context_t *context);
static int LogError(const cg_context_t *context, const char *format, ...);
static int ParseVsopBodyName(const cg_context_t *context, const char *name, vsop_body_t *body);

int GenerateCode(
    cg_language_t language,
    const char *outCodeFileName,
    const char *inTemplateFileName,
    const char *dataPath)
{
    int error = 1;
    cg_context_t context;
    FILE *infile = NULL;
    char line[CG_MAX_LINE_LENGTH];
    const char *tail;

    memset(&context, 0, sizeof(cg_context_t));
    context.language = language;
    context.inFileName = inTemplateFileName;
    context.datapath = dataPath;

    infile = fopen(inTemplateFileName, "rt");
    if (infile == NULL)
    {
        fprintf(stderr, "GenerateCode: Cannot open input template file: %s\n", inTemplateFileName);
        goto fail;
    }

    context.outfile = fopen(outCodeFileName, "wt");
    if (context.outfile == NULL)
    {
        fprintf(stderr, "GenerateCode: Cannot open output file: %s\n", outCodeFileName);
        goto fail;
    }

    for (context.lnum=1; fgets(line, sizeof(line), infile); ++context.lnum)
    {
        /* Look for $ASTRO_...(...)$ directives */
        if (ScanDirective(line, &context.verb, &context.args, &tail))
        {
            fputs(line, context.outfile);
            error = ProcessDirective(&context);
            if (error) goto fail;
            fputs(tail, context.outfile);
        }
        else
        {
            /* Not a directive, so copy to output verbatim. */
            fputs(line, context.outfile);
        }
    }
    error = 0;

fail:
    if (infile != NULL) fclose(infile);
    if (context.outfile != NULL) fclose(context.outfile);
    if (error) remove(outCodeFileName);
    return error;
}

static int IndexOf(const char *text, int offset, const char *pattern)
{
    int i, j;
    for (i=offset; text[i]; ++i)
    {
        for (j=0; pattern[j] && (pattern[j] == text[i+j]); ++j);
        if (!pattern[j]) return i;
    }
    return -1;
}

static int ScanDirective(
    char *line,
    const char **verb_text,
    const char **args_text,
    const char **tail_text)
{
    static const char prefix[] = "$ASTRO_";
    static const char delim[] = "(";
    static const char suffix[] = ")";
    int prefix_length = (int)strlen(prefix);
    int delim_length = (int)strlen(delim);
    int suffix_length = (int)strlen(suffix);
    int prefix_index, delim_index, suffix_index;

    *verb_text = NULL;
    *args_text = NULL;
    *tail_text = NULL;

    /*
        Search for a directive anywhere in the line that looks like this:

        $ASTRO_...(...)$

        If that pattern is found, set *verb_text to point to the
        reset of the verb after the prefix, set *args_text to the
        string after the delimiter, and return 1.
    */

    /* Search for the front of the directive in the line. */
    prefix_index = IndexOf(line, 0, prefix);
    if (prefix_index < 0)
        return 0;    /* Not found, so bail out. */

    /* Search for the delimiter between the verb and the arguments. */
    delim_index = IndexOf(line, prefix_index + prefix_length, delim);
    if (delim_index < 0)
        return 0;   /* Not a pattern match, so bail out. */

    /* Search for the arguments terminator. */
    suffix_index = IndexOf(line, delim_index + delim_length, suffix);
    if (suffix_index < 0)
        return 0;   /* Not a pattern match, so bail out. */

    /* We have found a directive, so now we are allowed to change the contents of the line. */

    line[prefix_index] = '\0';     /* terminate any leading text in the line */
    line[delim_index] = '\0';      /* terminate the prefix string */
    line[suffix_index] = '\0';     /* terminate the argument string */

    /* skip over the "$ASTRO_" before the verb */
    *verb_text = &line[prefix_index + prefix_length];

    /* skip over the delimiter between the verb and the arguments */
    *args_text = &line[delim_index + delim_length];

    /* return a pointer to the part of the input line after the directive */
    *tail_text = &line[suffix_index + suffix_length];

    return 1;
}

static int ParseVsopBodyName(const cg_context_t *context, const char *name, vsop_body_t *body)
{
    if (!strcmp(name, "Sun"))       { *body = VSOP_SUN;     return 0; }
    if (!strcmp(name, "Mercury"))   { *body = VSOP_MERCURY; return 0; }
    if (!strcmp(name, "Venus"))     { *body = VSOP_VENUS;   return 0; }
    if (!strcmp(name, "EMB"))       { *body = VSOP_EMB;     return 0; }
    if (!strcmp(name, "Earth"))     { *body = VSOP_EARTH;   return 0; }
    if (!strcmp(name, "Mars"))      { *body = VSOP_MARS;    return 0; }
    if (!strcmp(name, "Jupiter"))   { *body = VSOP_JUPITER; return 0; }
    if (!strcmp(name, "Saturn"))    { *body = VSOP_SATURN;  return 0; }
    if (!strcmp(name, "Uranus"))    { *body = VSOP_URANUS;  return 0; }
    if (!strcmp(name, "Neptune"))   { *body = VSOP_NEPTUNE; return 0; }

    *body = VSOP_INVALID_BODY;
    return LogError(context, "Unknown VSOP body name '%s'", name);
}

static int ListVsop(cg_context_t *context)
{
    int error;
    const char *name;
    vsop_body_t body;
    vsop_model_t model;
    int check_length;
    int k, s, i;
    char filename[100];

    VsopInit(&model);

    name = context->args;
    error = ParseVsopBodyName(context, name, &body);
    if (error) goto fail;

    check_length = snprintf(filename, sizeof(filename), "%s/vsop_%d.txt", context->datapath, (int)body);
    if (check_length < 0 || check_length != (int)strlen(filename))
    {
        error = LogError(context, "VSOP model filename is too long!");
        goto fail;
    }

    error = VsopReadTrunc(&model, filename);
    if (error) goto fail;

    /* Represent the VSOP model using JavaScript syntax. */
    fprintf(context->outfile, "[\n");
    for (k=0; k < model.ncoords; ++k)
    {
        const vsop_formula_t *formula = &model.formula[k];

        fprintf(context->outfile, "  [\n");
        for (s=0; s < formula->nseries_total; ++s)
        {
            const vsop_series_t *series = &formula->series[s];
            fprintf(context->outfile, "    [\n");
            for (i=0; i < series->nterms_total; ++i)
            {
                const vsop_term_t *term = &series->term[i];
                fprintf(context->outfile, "      [%0.11lf, %0.11lf, %0.11lf]%s\n",
                    term->amplitude,
                    term->phase,
                    term->frequency,
                    (i+1 < series->nterms_total) ? "," : "");
            }
            fprintf(context->outfile, "    ]%s\n", (s+1 < formula->nseries_total) ? "," : "");
        }
        fprintf(context->outfile, "  ]%s\n", (k+1 < model.ncoords) ? "," : "");
    }
    fprintf(context->outfile, "]");

fail:
    VsopFreeModel(&model);
    return error;
}

static int CVsop_Series(cg_context_t *context, const vsop_series_t *series, const char *varprefix, int s)
{
    int i;

    if (series->nterms_total > 0)
    {
        fprintf(context->outfile, "static const vsop_term_t %s_%d[] =\n{\n", varprefix, s);
        for (i = 0; i < series->nterms_total; ++i)
        {
            const vsop_term_t *term = &series->term[i];

            fprintf(context->outfile, "    { %0.11lf, %0.11lf, %0.11lf }%s\n",
                term->amplitude,
                term->phase,
                term->frequency,
                (i + 1 < series->nterms_total) ? "," : "");
        }
        fprintf(context->outfile, "};\n\n");
    }

    return 0;
}

static int CVsop_Formula(cg_context_t *context, const vsop_formula_t *formula, const char *coord_name, const char *body_name)
{
    int error = 0;
    int s;
    char varprefix[100];
    char sname[120];

    snprintf(varprefix, sizeof(varprefix), "vsop_%s_%s", coord_name, body_name);

    for (s=0; s < formula->nseries_total; ++s)
        CHECK(CVsop_Series(context, &formula->series[s], varprefix, s));

    fprintf(context->outfile, "static const vsop_series_t %s[] =\n{\n", varprefix);
    for (s=0; s < formula->nseries_total; ++s)
    {
        if (formula->series[s].nterms_total == 0)
            strcpy(sname, "NULL");
        else
            snprintf(sname, sizeof(sname), "%s_%d", varprefix, s);

        fprintf(context->outfile, "    { %d, %s }%s\n",
            formula->series[s].nterms_total,
            sname,
            (s+1 < formula->nseries_total) ? "," : "");
    }
    fprintf(context->outfile, "};\n\n");

fail:
    return error;
}

static int CVsop(cg_context_t *context)
{
    int error;
    const char *name;
    vsop_body_t body;
    vsop_model_t model;
    int check_length;
    int k;
    char filename[100];
    const char *coord_name[3] = { "lon", "lat", "rad" };

    VsopInit(&model);

    name = context->args;
    CHECK(ParseVsopBodyName(context, name, &body));

    check_length = snprintf(filename, sizeof(filename), "%s/vsop_%d.txt", context->datapath, (int)body);
    if (check_length < 0 || check_length != (int)strlen(filename))
        CHECK(LogError(context, "VSOP model filename is too long!"));

    CHECK(VsopReadTrunc(&model, filename));

    for (k=0; k < model.ncoords; ++k)
        CHECK(CVsop_Formula(context, &model.formula[k], coord_name[k], name));

fail:
    VsopFreeModel(&model);
    return error;
}

static int CsharpVsop_Series(cg_context_t *context, const vsop_series_t *series, const char *varprefix, int s)
{
    int i;

    fprintf(context->outfile, "        private static readonly vsop_term_t[] %s_%d = new vsop_term_t[]\n        {\n", varprefix, s);
    for (i = 0; i < series->nterms_total; ++i)
    {
        const vsop_term_t *term = &series->term[i];

        fprintf(context->outfile, "            new vsop_term_t(%0.11lf, %0.11lf, %0.11lf)%s\n",
            term->amplitude,
            term->phase,
            term->frequency,
            (i + 1 < series->nterms_total) ? "," : "");
    }
    fprintf(context->outfile, "        };\n\n");

    return 0;
}

static int CsharpVsop_Formula(cg_context_t *context, const vsop_formula_t *formula, const char *coord_name, const char *body_name)
{
    int error = 0;
    int s;
    char varprefix[100];
    char sname[120];

    snprintf(varprefix, sizeof(varprefix), "vsop_%s_%s", coord_name, body_name);

    for (s=0; s < formula->nseries_total; ++s)
        CHECK(CsharpVsop_Series(context, &formula->series[s], varprefix, s));

    fprintf(context->outfile, "        private static readonly vsop_series_t[] %s = new vsop_series_t[]\n        {\n", varprefix);
    for (s=0; s < formula->nseries_total; ++s)
    {
        snprintf(sname, sizeof(sname), "%s_%d", varprefix, s);

        fprintf(context->outfile, "            new vsop_series_t(%s)%s\n",
            sname,
            (s+1 < formula->nseries_total) ? "," : "");
    }
    fprintf(context->outfile, "        };\n\n");

fail:
    return error;
}

static int CsharpVsop(cg_context_t *context)
{
    int error;
    const char *name;
    vsop_body_t body;
    vsop_model_t model;
    int check_length;
    int k;
    char filename[100];
    const char *coord_name[3] = { "lon", "lat", "rad" };

    VsopInit(&model);

    name = context->args;
    CHECK(ParseVsopBodyName(context, name, &body));

    check_length = snprintf(filename, sizeof(filename), "%s/vsop_%d.txt", context->datapath, (int)body);
    if (check_length < 0 || check_length != (int)strlen(filename))
        CHECK(LogError(context, "VSOP model filename is too long!"));

    CHECK(VsopReadTrunc(&model, filename));

    for (k=0; k < model.ncoords; ++k)
        CHECK(CsharpVsop_Formula(context, &model.formula[k], coord_name[k], name));

fail:
    VsopFreeModel(&model);
    return error;
}

#define UDEF(expr)  ((u = (expr)), (u2 = u*u), (u3 = u*u2), (u4 = u2*u2), (u5 = u2*u3), (u6 = u3*u3), (u7 = u3*u4))

double ExtrapolatedDeltaT(int year)
{
    double dt=NAN, y;
    double u, u2, u3, u4, u5, u6, u7;

    /*
        Fred Espenak writes about Delta-T generically here:
        https://eclipse.gsfc.nasa.gov/SEhelp/deltaT.html
        https://eclipse.gsfc.nasa.gov/SEhelp/deltat2004.html

        He provides polynomial approximations for distant years here:
        https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
    */

    y = (double)year + (0.5/12.0);      /* Espenak's "decimal year" value y for January of the given year. */

    if (year < -500)
    {
        u = (year - 1820.0) / 100.0;
        dt = -20.0 + (32.0 * u*u);
    }
    else if (year < 500)
    {
        UDEF(y / 100.0);
        dt = 10583.6 - 1014.41*u + 33.78311*u2 - 5.952053*u3
		    - 0.1798452*u4 + 0.022174192*u5 + 0.0090316521*u6;
    }
    else if (year < 1600)
    {
        UDEF((y - 1000.0) / 100.0);
        dt = 1574.2 - 556.01*u + 71.23472*u2 + 0.319781*u3
		     - 0.8503463*u4 - 0.005050998*u5 + 0.0083572073*u6;
    }
    else if (year < 1700)
    {
        UDEF(y - 1600);
        dt = 120 - 0.9808*u - 0.01532*u2 + u3/7129.0;
    }
    else if (year < 1800)
    {
        UDEF(y - 1700);
        dt = 8.83 + 0.1603*u - 0.0059285*u2 + 0.00013336*u3 - u4/1174000;
    }
    else if (year < 1860)
    {
        UDEF(y - 1800);
        dt = 13.72 - 0.332447*u + 0.0068612*u2 + 0.0041116*u3 - 0.00037436*u4
		     + 0.0000121272*u5 - 0.0000001699*u6 + 0.000000000875*u7;
    }
    else if (year < 1900)
    {
        UDEF(y - 1860);
        dt = 7.62 + 0.5737*u - 0.251754*u2 + 0.01680668*u3
		     -0.0004473624*u4 + u5/233174;
    }
    else if (year < 1920)
    {
        UDEF(y - 1900);
        dt = -2.79 + 1.494119*u - 0.0598939*u2 + 0.0061966*u3 - 0.000197*u4;
    }
    else if (year < 1941)
    {
        UDEF(y - 1920);
        dt = 21.20 + 0.84493*u - 0.076100*u2 + 0.0020936*u3;
    }
    else if (year < 1961)
    {
        UDEF(y - 1950);
        dt = 29.07 + 0.407*u - u2/233 + u3/2547;
    }
    else if (year < 1986)
    {
        UDEF(y - 1975);
        dt = 45.45 + 1.067*u - u2/260 - u3/718;
    }
    else if (year < 2005)
    {
        UDEF(y - 2000);
        dt = 63.86 + 0.3345*u - 0.060374*u2 + 0.0017275*u3 + 0.000651814*u4
		     + 0.00002373599*u5;
    }
    else if (year < 2050)
    {
        UDEF(y - 2000);
        dt = 62.92 + 0.32217*u + 0.005589*u2;
    }
    else if (year < 2150)
    {
        u = (y-1820)/100;
        dt = -20.0 + 32.0*u*u - 0.5628*(2150.0 - y);
    }
    else    /* all years after 2150 */
    {
        u = (year - 1820.0) / 100.0;
        dt = -20.0 + 32.0*u*u;
    }

    return dt;
}


static int ScanRealArray(
    cg_context_t *context,
    const char *filename,
    int lnum,
    char *line,
    int numExpected,
    double *data)
{
    int i, t, len, inspace;
    char *token[MAX_DATA_PER_LINE];

    if (numExpected < 1 || numExpected > MAX_DATA_PER_LINE)
        return LogError(context, "Invalid value for numExpected=%d\n", numExpected);

    /* Split the line into space delimited tokens. */
    len = (int)strlen(line);
    inspace = 1;
    t = 0;
    for (i=0; i < len; ++i)
    {
        if (line[i] == ' ' || line[i] == '\t' || line[i] == '\r' || line[i] == '\n')
        {
            if (!inspace)
            {
                line[i] = '\0';     /* terminate the previous token */
                inspace = 1;
            }
        }
        else
        {
            if (inspace)
            {
                /* we just found the beginning of a new token. */
                if (t < IAU_DATA_PER_ROW)
                {
                    token[t++] = &line[i];
                    inspace = 0;
                }
                else
                    return LogError(context, "ScanRealArray(%s %d): too many data on line.", filename, lnum);
            }
        }
    }

    /* Verify there are the correct number of tokens. */
    if (t != numExpected)
        return LogError(context, "ScanRealArray(%s %d): found %d data, but expected %d\n", filename, lnum, t, numExpected);

    /* Parse each token as a floating point number. */
    for (t=0; t < numExpected; ++t)
        if (1 != sscanf(token[t], "%lf", &data[t]))
            return LogError(context, "ScanRealArray(%s %d): invalid floating point token '%s'\n", filename, lnum, token[t]);

    return 0;   /* successful parse */
}

static int OptimizeConst(cg_context_t *context, char *buffer, size_t size, double c, const char *v)
{
    int nprinted;
    const char *op;

    if (c == 0.0)
    {
        buffer[0] = '\0';
        return 0;
    }

    if (c < 0.0)
    {
        op = " - ";
        c *= -1.0;
    }
    else
        op = " + ";

    if (c == 1.0)
        nprinted = snprintf(buffer, size, "%s%s", op, v);
    else
        nprinted = snprintf(buffer, size, "%s%0.1lf*%s", op, c, v);

    if (nprinted >= (int)size)
        return LogError(context, "OptimizeConst: print buffer overflowed.");

    return 0;
}

static int OptimizeLinear(cg_context_t *context, char *buffer, size_t size, double a, double b)
{
    int nprinted;

    if (b == 0.0)
        nprinted = snprintf(buffer, size, "%0.1lf", a);
    else if (a == 0.0)
        nprinted = snprintf(buffer, size, "%0.1lf*t", b);
    else if (b < 0.0)
        nprinted = snprintf(buffer, size, "%0.1lf - %0.1lf*t", a, -b);
    else
        nprinted = snprintf(buffer, size, "%0.1lf + %0.1lf*t", a, b);

    if ((size_t)nprinted >= size)
        return LogError(context, "OptimizeLinear: print buffer overflowed.");

    return 0;
}

static int OptIauPython(cg_context_t *context, const double *data)
{
    static const char * const nv[] = {"el", "elp", "f", "d", "om"};    /* variable names */
    int first;
    double n;
    int i;
    const char *op;
    char dotprod[200];
    char term[40];
    char linear[100];
    char cpart[100];
    int nprinted;
    int lonevar;

    fprintf(context->outfile, "\n");    /* must start on new line to maintain correct indentation in Python source. */

    /* Optimize dot product of data[0]..data[4] with nv[]. */
    first = 1;
    dotprod[0] = '\0';
    for (i=0; i<5; ++i)
    {
        n = data[i];
        if (n != 0.0)
        {
            if (n < 0.0)
            {
                n *= -1.0;
                op = first ? "-" : " - ";
            }
            else
                op = first ? "" : " + ";

            if (n == 1.0)
                nprinted = snprintf(term, sizeof(term), "%s%s", op, nv[i]);
            else
                nprinted = snprintf(term, sizeof(term), "%s%0.1lf*%s", op, n, nv[i]);

            if (nprinted >= sizeof(term))
                return LogError(context, "Truncated iau2000b term.");

            if (nprinted + strlen(dotprod) >= sizeof(dotprod))
                return LogError(context, "Dot product overflow in iau2000b formula.");

            strcat(dotprod, term);
            first = 0;
        }
    }

    /* Did we print a lone variable, e.g. "elp"? */
    lonevar = 0;
    for (i=0; i<5 && !lonevar; ++i)
        if (!strcmp(dotprod, nv[i]))
            lonevar = 1;

    if (lonevar)
    {
        fprintf(context->outfile, "        sarg = math.sin(%s)\n", dotprod);
        fprintf(context->outfile, "        carg = math.cos(%s)\n", dotprod);
    }
    else
    {
        fprintf(context->outfile, "        arg = %s\n", dotprod);
        fprintf(context->outfile, "        sarg = math.sin(arg)\n");
        fprintf(context->outfile, "        carg = math.cos(arg)\n");
    }

    if (OptimizeLinear(context, linear, sizeof(linear), data[5], data[6])) return 1;
    if (OptimizeConst(context, cpart, sizeof(cpart), data[7], "carg")) return 1;
    fprintf(context->outfile, "        dp += (%s)*sarg%s\n", linear, cpart);

    if (OptimizeLinear(context, linear, sizeof(linear), data[8], data[9])) return 1;
    if (OptimizeConst(context, cpart, sizeof(cpart), data[10], "sarg")) return 1;
    fprintf(context->outfile, "        de += (%s)*carg%s\n", linear, cpart);

    fprintf(context->outfile, "\n");
    return 0;
}

static int OptIauC(cg_context_t *context, int lnum, const double *data)
{
    int i;

    fprintf(context->outfile, "        { { %2.0lf", data[0]);

    for (i=1; i < 5; ++i)
        fprintf(context->outfile, ", %2.0lf", data[i]);

    fprintf(context->outfile, " }, { %12.0lf", data[i]);
    for (++i; i < 11; ++i)
        fprintf(context->outfile, ", %12.0lf", data[i]);

    fprintf(context->outfile, " } }");
    if (lnum != IAU_DATA_NUM_ROWS)
        fprintf(context->outfile, ",");
    fprintf(context->outfile, "\n");

    return 0;
}

static int OptIauCsharp(cg_context_t *context, int lnum, const double *data)
{
    int i;

    fprintf(context->outfile, "        new iau_row_t { nals0 = %2.0lf", data[0]);

    for (i=1; i < 5; ++i)
        fprintf(context->outfile, ", nals%d = %2.0lf", i, data[i]);

    fprintf(context->outfile, " , cls%d = %12.0lf", i-5, data[i]);
    for (++i; i < 11; ++i)
        fprintf(context->outfile, ", cls%d = %12.0lf", i-5, data[i]);

    fprintf(context->outfile, " }");
    if (lnum != IAU_DATA_NUM_ROWS)
        fprintf(context->outfile, ",");
    fprintf(context->outfile, "\n");

    return 0;
}

static int OptIauJS(cg_context_t *context, int lnum, const double *data)
{
    int i;

    fprintf(context->outfile, "    [ [ %2.0lf", data[0]);

    for (i=1; i < 5; ++i)
        fprintf(context->outfile, ", %2.0lf", data[i]);

    fprintf(context->outfile, " ], [ %12.0lf", data[i]);
    for (++i; i < 11; ++i)
        fprintf(context->outfile, ", %12.0lf", data[i]);

    fprintf(context->outfile, " ] ]");
    if (lnum != IAU_DATA_NUM_ROWS)
        fprintf(context->outfile, ",");
    fprintf(context->outfile, "\n");

    return 0;
}

static int OptIauData(cg_context_t *context)
{
    int error = 1;
    int lnum;
    FILE *infile;
    const char *filename;
    char line[100];
    double data[IAU_DATA_PER_ROW];

    filename = "model_data/iau2000b.txt";
    infile = fopen(filename, "rt");
    if (infile == NULL) goto fail;

    fprintf(context->outfile, "\n");

    lnum = 0;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;
        CHECK(ScanRealArray(context, filename, lnum, line, IAU_DATA_PER_ROW, data));
        switch (context->language)
        {
        case CODEGEN_LANGUAGE_PYTHON:
            CHECK(OptIauPython(context, data));
            break;

        case CODEGEN_LANGUAGE_C:
            CHECK(OptIauC(context, lnum, data));
            break;

        case CODEGEN_LANGUAGE_CSHARP:
            CHECK(OptIauCsharp(context, lnum, data));
            break;

        case CODEGEN_LANGUAGE_JS:
            CHECK(OptIauJS(context, lnum, data));
            break;

        default:
            CHECK(LogError(context, "OptIauData: Unsupported language %d", context->language));
        }
    }

    if (lnum != IAU_DATA_NUM_ROWS)
        CHECK(LogError(context, "OptIauData: expected %d rows, found %d\n", IAU_DATA_NUM_ROWS, lnum));

    error = 0;
fail:
    if (infile == NULL)
        error = LogError(context, "Cannot open input file: %s", filename);
    else
        fclose(infile);

    return error;
}

static int OptAddSolPython(
    cg_context_t *context,
    double cl, double cs, double cg, double cp, double p, double q, double r, double s)
{
    const char *op;

    fprintf(context->outfile, "\n    # AddSol(%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf)\n", cl, cs, cg, cp, p, q, r, s);

    op = "";
    fprintf(context->outfile, "    z = ");
    if (p != 0.0)
    {
        fprintf(context->outfile, "ex[%0.0lf][1]", p);
        op = " * ";
    }
    if (q != 0.0)
    {
        fprintf(context->outfile, "%sex[%0.0lf][2]", op, q);
        op = " * ";
    }
    if (r != 0.0)
    {
        fprintf(context->outfile, "%sex[%0.0lf][3]", op, r);
        op = " * ";
    }
    if (s != 0.0)
    {
        fprintf(context->outfile, "%sex[%0.0lf][4]", op, s);
    }
    fprintf(context->outfile, "\n");

    if (cl != 0.0)
        fprintf(context->outfile, "    DLAM  += %0.3lf * z.imag\n", cl);

    if (cs != 0.0)
        fprintf(context->outfile, "    DS    += %0.2lf * z.imag\n", cs);

    if (cg != 0.0)
        fprintf(context->outfile, "    GAM1C += %0.3lf * z.real\n", cg);

    if (cp != 0.0)
        fprintf(context->outfile, "    SINPI += %0.4lf * z.real\n", cp);

    return 0;
}

static int OptAddSolC(cg_context_t *context, const double *data)
{
    int i;

    fprintf(context->outfile, "    AddSol(ctx");

    for (i=0; i < 4; ++i)
        fprintf(context->outfile, ",%11.4lf", data[i]);

    for(; i < 8; ++i)
        fprintf(context->outfile, ",%2.0lf", data[i]);

    fprintf(context->outfile, ");\n");
    return 0;
}

static int OptAddSolCsharp(cg_context_t *context, const double *data)
{
    int i;

    fprintf(context->outfile, "            AddSol(");

    for (i=0; i < 4; ++i)
        fprintf(context->outfile, "%11.4lf,", data[i]);

    for(; i < 8; ++i)
        fprintf(context->outfile, "%2.0lf%s", data[i], (i < 7) ? "," : ");\n");

    return 0;
}

static int OptAddSolJS(cg_context_t *context, const double *data)
{
    int i;

    fprintf(context->outfile, "    AddSol(%11.4lf", data[0]);

    for (i=1; i < 4; ++i)
        fprintf(context->outfile, ",%11.4lf", data[i]);

    for(; i < 8; ++i)
        fprintf(context->outfile, ",%2.0lf", data[i]);

    fprintf(context->outfile, ");\n");
    return 0;
}

static int OptAddSol(cg_context_t *context)
{
    int error = 1;
    FILE *infile;
    int lnum;
    const char *filename = "model_data/addsol.txt";
    char line[200];
    double data[ADDSOL_DATA_PER_ROW];

    infile = fopen(filename, "rt");
    if (infile == NULL) goto fail;

    fprintf(context->outfile, "\n");

    lnum = 0;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;
        CHECK(ScanRealArray(context, filename, lnum, line, ADDSOL_DATA_PER_ROW, data));
        switch (context->language)
        {
        case CODEGEN_LANGUAGE_PYTHON:
            CHECK(OptAddSolPython(context, data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7]));
            break;

        case CODEGEN_LANGUAGE_C:
            CHECK(OptAddSolC(context, data));
            break;

        case CODEGEN_LANGUAGE_CSHARP:
            CHECK(OptAddSolCsharp(context, data));
            break;

        case CODEGEN_LANGUAGE_JS:
            CHECK(OptAddSolJS(context, data));
            break;

        default:
            error = LogError(context, "OptAddSol: Unsupported language %d\n", context->language);
            goto fail;
        }
    }

    error = 0;
fail:
    if (infile == NULL)
        error = LogError(context, "Cannot open input file: %s", filename);
    else
        fclose(infile);

    return error;
}

static double RoundAngle(double x)
{
    /*
        In the JavaScript code, it is more compact to store
        (quarter arcminutes) / 10 than degrees and sidereal hours.

        Before:
            -rw-r--r-- 1 don don 249078 May  4 10:06 astronomy.js
            -rw-r--r-- 1 don don  86423 May  4 10:06 astronomy.min.js

        After:
            -rw-r--r-- 1 don don 237503 May  4 12:17 astronomy.js
            -rw-r--r-- 1 don don  80904 May  4 12:17 astronomy.min.js
    */
    return round((4.0 * 60.0) * x) / 10.0;
}

#define NUM_CONSTELLATIONS 88

static int ConstellationData(cg_context_t *context)
{
    int error;
    const char *dictFileName   = "constellation/constellation.dictionary";
    const char *borderFileName = "constellation/constellation.borders";
    FILE *infile;
    int lnum, len, index;
    char line[100];
    char dict[NUM_CONSTELLATIONS][25];
    char *d;
    double dec, ra_lo, ra_hi;
    char symbol[4];

    /* Generate the table of constellation symbols and names. */

    infile = fopen(dictFileName, "rt");
    if (infile == NULL)
        CHECK(LogError(context, "Cannot open input file: %s", dictFileName));

    lnum = 0;
    while (fgets(line, sizeof(line), infile))
    {
        if (lnum == NUM_CONSTELLATIONS)
            CHECK(LogError(context, "Too many constellations in file: %s", dictFileName));

        d = dict[lnum++];
        len = (int)strlen(line);
        if (len < 26)
            CHECK(LogError(context, "File %s line %d is too short", dictFileName, lnum));

        if (!isupper(line[0]) || !isalpha(line[1]) || !isalpha(line[2]) || (line[3] != ' '))
            CHECK(LogError(context, "File %s line %d has invalid constellation symbol", dictFileName, lnum));

        /* Tokenize the symbol and name and store them in the 'dict' array. */
        /* We need 'dict' for symbol lookups when we process the boundary file below. */
        memcpy(d, line, 24);
        for (index=23; index >= 0 && d[index] == ' '; --index)
            d[index] = '\0';
        d[3] = '\0';

        len = (int)strlen(d+4);  /* capture variable length of constellation name for formatting columns nicely. */

        switch (context->language)
        {
        case CODEGEN_LANGUAGE_C:
            if (lnum == 1)
            {
                fprintf(context->outfile, "#define NUM_CONSTELLATIONS   %d\n\n", NUM_CONSTELLATIONS);
                fprintf(context->outfile, "static const constel_info_t ConstelInfo[] = {\n    ");
            }
            else
                fprintf(context->outfile, ",   ");
            fprintf(context->outfile, "/* %2d */ { \"%s\", \"%s\"%*s }\n", lnum-1, d, d+4, (20-len), "");
            break;

        case CODEGEN_LANGUAGE_CSHARP:
            if (lnum == 1)
            {
                fprintf(context->outfile, "        private static readonly constel_info_t[] ConstelNames = new constel_info_t[]\n");
                fprintf(context->outfile, "        {\n");
                fprintf(context->outfile, "            ");
            }
            else
                fprintf(context->outfile, "        ,   ");
            fprintf(context->outfile, "new constel_info_t(\"%s\", \"%s\"%*s)  // %2d\n", d, d+4, (20-len), "", lnum-1);
            break;

        case CODEGEN_LANGUAGE_JS:
            if (lnum == 1)
                fprintf(context->outfile, "const ConstelNames = [\n ");
            else
                fprintf(context->outfile, ",");
            fprintf(context->outfile, "   ['%s', '%s'%*s]  // %2d\n", d, d+4, (20-len), "", lnum-1);
            break;

        case CODEGEN_LANGUAGE_PYTHON:
            if (lnum == 1)
                fprintf(context->outfile, "_ConstelNames = (\n ");
            else
                fprintf(context->outfile, ",");
            fprintf(context->outfile, "   ('%s', '%s'%*s)  # %2d\n", d, d+4, (20-len), "", lnum-1);
            break;

        default:
             CHECK(LogError(context, "ConstellationData(1): Unsupported language %d", context->language));
        }
    }

    if (lnum != NUM_CONSTELLATIONS)
        CHECK(LogError(context, "Not enough constellations in file %s : expected %d", dictFileName, NUM_CONSTELLATIONS));

    fclose(infile);
    infile = fopen(borderFileName, "rt");
    if (infile == NULL)
        CHECK(LogError(context, "Cannot open input file: %s", borderFileName));

    /* End the name table. */
    /* Generate the table of constellation boundary lines. */

    switch (context->language)
    {
    case CODEGEN_LANGUAGE_C:
        fprintf(context->outfile, "};\n\n");
        fprintf(context->outfile, "static const constel_boundary_t ConstelBounds[] = {\n");
        break;

    case CODEGEN_LANGUAGE_CSHARP:
        fprintf(context->outfile, "        };\n\n");
        fprintf(context->outfile, "        private static readonly constel_boundary_t[] ConstelBounds = new constel_boundary_t[]\n");
        fprintf(context->outfile, "        {\n");
        break;

    case CODEGEN_LANGUAGE_JS:
        fprintf(context->outfile, "];\n\n");
        fprintf(context->outfile, "const ConstelBounds = [\n");
        break;

    case CODEGEN_LANGUAGE_PYTHON:
        fprintf(context->outfile, ")\n\n");
        fprintf(context->outfile, "_ConstelBounds = (\n");
        break;

    default:
        CHECK(LogError(context, "ConstellationData(2): Unsupported language %d", context->language));
    }

    /* Translate the constellation boundaries. */
    lnum = 0;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;
        if (4 != sscanf(line, "%lf %lf %lf %3s", &ra_lo, &ra_hi, &dec, symbol))
            CHECK(LogError(context, "Invalid constellation data at %s line %d", borderFileName, lnum));

        /* Search for the constellation symbol in the name table, to find its array index. */
        for (index=0; index < NUM_CONSTELLATIONS; ++index)
            if (0 == strcmp(symbol, dict[index]))
                break;

        if (index == NUM_CONSTELLATIONS)
            CHECK(LogError(context, "Invalid constellation symbol '%s' at %s line %d", symbol, borderFileName, lnum));

        /* For compact encoding, store each angle in tenths of quarter-arcminutes. */
        ra_lo = RoundAngle(ra_lo);
        ra_hi = RoundAngle(ra_hi);
        dec = RoundAngle(dec);

        switch (context->language)
        {
        case CODEGEN_LANGUAGE_C:
            fprintf(context->outfile, "%c   { %2d, %6lg, %6lg, %6lg }    /* %s */\n",
                ((lnum == 1) ? ' ' : ','),
                index, ra_lo, ra_hi, dec, symbol);
            break;

        case CODEGEN_LANGUAGE_CSHARP:
            fprintf(context->outfile, "        %c   new constel_boundary_t(%2d, %6lg, %6lg, %6lg)    // %s\n",
                ((lnum == 1) ? ' ' : ','),
                index, ra_lo, ra_hi, dec, symbol);
            break;

        case CODEGEN_LANGUAGE_JS:
            fprintf(context->outfile, "%c   [ %2d, %6lg, %6lg, %6lg ]    // %s\n",
                ((lnum == 1) ? ' ' : ','),
                index, ra_lo, ra_hi, dec, symbol);
            break;

        case CODEGEN_LANGUAGE_PYTHON:
            fprintf(context->outfile, "%c   ( %2d, %6lg, %6lg, %6lg )    # %s\n",
                ((lnum == 1) ? ' ' : ','),
                index, ra_lo, ra_hi, dec, symbol);
            break;

        default:
            CHECK(LogError(context, "ConstellationData(3): Unsupported language %d", context->language));
        }
    }

    switch (context->language)
    {
    case CODEGEN_LANGUAGE_C:
        fprintf(context->outfile, "};\n\n");
        fprintf(context->outfile, "#define NUM_CONSTEL_BOUNDARIES  %d\n\n", lnum);
        break;

    case CODEGEN_LANGUAGE_CSHARP:
        fprintf(context->outfile, "        };\n\n");
        break;

    case CODEGEN_LANGUAGE_JS:
        fprintf(context->outfile, "];\n\n");
        break;

    case CODEGEN_LANGUAGE_PYTHON:
        fprintf(context->outfile, ")\n\n");
        break;

    default:
        CHECK(LogError(context, "ConstellationData(4): Unsupported language %d", context->language));
    }

    error = 0;
fail:
    if (infile) fclose(infile);
    return error;
}


static int PlutoConstants_C(cg_context_t *context)
{
    fprintf(context->outfile, "#define PLUTO_NUM_STATES  %d\n", PLUTO_NUM_STATES);
    fprintf(context->outfile, "#define PLUTO_TIME_STEP   %d\n", PLUTO_TIME_STEP);
    fprintf(context->outfile, "#define PLUTO_DT          %d\n\n", PLUTO_DT);
    fprintf(context->outfile, "#define PLUTO_NSTEPS      %d\n\n", PLUTO_NSTEPS);
    return 0;
}


static int PlutoStateTable_C(cg_context_t *context, const top_model_t *model)
{
    int error = 1;
    double tt;
    top_rectangular_t equ;
    int i;

    fprintf(context->outfile, "static const body_state_t PlutoStateTable[] =\n");
    fprintf(context->outfile, "{\n");

    for (i=0; i < PLUTO_NUM_STATES; ++i)
    {
        tt = i*PLUTO_TIME_STEP + PLUTO_TT1;
        CHECK(TopPosition(model, tt, &equ));

        fprintf(context->outfile,
            "%c   { %10.1lf, {%17.13lf, %17.13lf, %17.13lf}, {%20.13le, %20.13le, %20.13le} }\n",
            (i==0 ? ' ' : ','),
            tt, equ.x, equ.y, equ.z, equ.vx, equ.vy, equ.vz);
    }

    fprintf(context->outfile, "}");

    error = 0;
fail:
    return error;
}


static int PlutoStateTable_CSharp(cg_context_t *context, const top_model_t *model)
{
    int error = 1;
    double tt;
    top_rectangular_t equ;
    int i;

    fprintf(context->outfile, "        private const int PLUTO_NUM_STATES = %d;\n", PLUTO_NUM_STATES);
    fprintf(context->outfile, "        private const int PLUTO_TIME_STEP  = %d;\n", PLUTO_TIME_STEP);
    fprintf(context->outfile, "        private const int PLUTO_DT         = %d;\n", PLUTO_DT);
    fprintf(context->outfile, "        private const int PLUTO_NSTEPS     = %d;\n", PLUTO_NSTEPS);
    fprintf(context->outfile, "\n");
    fprintf(context->outfile, "        private static readonly body_state_t[] PlutoStateTable = new body_state_t[]\n");
    fprintf(context->outfile, "        {\n");

    for (i=0; i < PLUTO_NUM_STATES; ++i)
    {
        tt = i*PLUTO_TIME_STEP + PLUTO_TT1;
        CHECK(TopPosition(model, tt, &equ));

        fprintf(context->outfile,
            "        %c   new body_state_t(%10.1lf, new TerseVector(%17.13lf, %17.13lf, %17.13lf), new TerseVector(%20.13le, %20.13le, %20.13le))\n",
            (i==0 ? ' ' : ','),
            tt, equ.x, equ.y, equ.z, equ.vx, equ.vy, equ.vz);
    }

    fprintf(context->outfile, "        }");

    error = 0;
fail:
    return error;
}


static int PlutoStateTable_JS(cg_context_t *context, const top_model_t *model)
{
    int error = 1;
    double tt;
    top_rectangular_t equ;
    int i;

    fprintf(context->outfile, "const PLUTO_NUM_STATES = %d;\n", PLUTO_NUM_STATES);
    fprintf(context->outfile, "const PLUTO_TIME_STEP  = %d;\n", PLUTO_TIME_STEP);
    fprintf(context->outfile, "const PLUTO_DT         = %d;\n", PLUTO_DT);
    fprintf(context->outfile, "const PLUTO_NSTEPS     = %d;\n", PLUTO_NSTEPS);
    fprintf(context->outfile, "\n");
    fprintf(context->outfile, "const PlutoStateTable: BodyStateTableEntry[] = [\n");

    for (i=0; i < PLUTO_NUM_STATES; ++i)
    {
        tt = i*PLUTO_TIME_STEP + PLUTO_TT1;
        CHECK(TopPosition(model, tt, &equ));

        fprintf(context->outfile,
            "%c   [%10.1lf, [%17.13lf, %17.13lf, %17.13lf], [%20.13le, %20.13le, %20.13le]]\n",
            (i==0 ? ' ' : ','),
            tt, equ.x, equ.y, equ.z, equ.vx, equ.vy, equ.vz);
    }

    fprintf(context->outfile, "];\n");

    error = 0;
fail:
    return error;
}


static int PlutoStateTable_Python(cg_context_t *context, const top_model_t *model)
{
    int error = 1;
    double tt;
    top_rectangular_t equ;
    int i;

    fprintf(context->outfile, "_PLUTO_NUM_STATES = %d\n", PLUTO_NUM_STATES);
    fprintf(context->outfile, "_PLUTO_TIME_STEP  = %d\n", PLUTO_TIME_STEP);
    fprintf(context->outfile, "_PLUTO_DT         = %d\n", PLUTO_DT);
    fprintf(context->outfile, "_PLUTO_NSTEPS     = %d\n", PLUTO_NSTEPS);
    fprintf(context->outfile, "\n");
    fprintf(context->outfile, "_PlutoStateTable = [\n");

    for (i=0; i < PLUTO_NUM_STATES; ++i)
    {
        tt = i*PLUTO_TIME_STEP + PLUTO_TT1;
        CHECK(TopPosition(model, tt, &equ));

        fprintf(context->outfile,
            "%c   [%10.1lf, [%17.13lf, %17.13lf, %17.13lf], [%20.13le, %20.13le, %20.13le]]\n",
            (i==0 ? ' ' : ','),
            tt, equ.x, equ.y, equ.z, equ.vx, equ.vy, equ.vz);
    }

    fprintf(context->outfile, "]\n");

    error = 0;
fail:
    return error;
}


static int PlutoStateTable(cg_context_t *context)
{
    int error = 1;
    top_model_t model;
    TopInitModel(&model);

    CHECK(TopLoadModel(&model, "TOP2013.dat", 9));

    switch (context->language)
    {
    case CODEGEN_LANGUAGE_C:
        CHECK(PlutoStateTable_C(context, &model));
        break;

    case CODEGEN_LANGUAGE_CSHARP:
        CHECK(PlutoStateTable_CSharp(context, &model));
        break;

    case CODEGEN_LANGUAGE_JS:
        CHECK(PlutoStateTable_JS(context, &model));
        break;

    case CODEGEN_LANGUAGE_PYTHON:
        CHECK(PlutoStateTable_Python(context, &model));
        break;

    default:
        CHECK(LogError(context, "PlutoStateTable: Unsupported target language %d", context->language));
    }

    error = 0;
fail:
    TopFreeModel(&model);
    return error;
}


/*-------------------------- begin Jupiter moons -----------------------------*/

static int ReadFixLine(char *line, size_t size, FILE *infile)
{
    if (!fgets(line, (int)size, infile))
        return 0;

    // Convert FORTRAN double-precision exponential notation like
    // 0.333688627964400D+06
    // to the format that sscanf() can handle:
    // 0.333688627964400e+06
    for (int i = 0; line[i]; ++i)
        if (isdigit(line[i]) && line[i+1] == 'D' && (line[i+2] == '+' || line[i+2] == '-') && isdigit(line[i+3]))
            line[i+1] = 'e';

    return 1;
}


#define REQUIRE_NSCAN(required) \
    do{if ((nscanned) != (required)) FAIL("LoadJupiterMoonModel(%s line %d): expected %d tokens, found %d\n", filename, lnum, required, nscanned);}while(0)


int LoadJupiterMoonModel(const char *filename, jupiter_moon_model_t *model)
{
    int error;
    FILE *infile;
    int lnum;
    char line[200];
    int nscanned, check, mindex, nterms, tindex, var, next_term_index, read_mean_longitude;
    vsop_series_t *series = NULL;

    infile = fopen(filename, "rt");
    if (infile == NULL)
        FAIL("LoadJupiterMoonModel: Cannot open file for read: %s\n", filename);

    lnum = 0;
    mindex = -1;
    nterms = 0;
    tindex = 0;
    var = 3;
    next_term_index = 0;
    read_mean_longitude = 0;
    while (ReadFixLine(line, sizeof(line), infile))
    {
        ++lnum;
        if (lnum == 19)
        {
            /* G(M+m), where M = Jupiter mass, m = moon mass. */
            nscanned = sscanf(line, "%lf %lf %lf %lf", &model->moon[0].mu, &model->moon[1].mu, &model->moon[2].mu, &model->moon[3].mu);
            REQUIRE_NSCAN(4);
        }
        else if (lnum == 20)
        {
            double cp, sp, ci, si;

            /* read angles that convert Jupiter equatorial coordinates into Earth J2000 equatorial coordinates. */
            nscanned = sscanf(line, "%lf %lf", &model->psi, &model->incl);
            REQUIRE_NSCAN(2);

            /* calculate the corresponding rotation matrix. */
            ci = cos(model->incl);
            si = sin(model->incl);
            cp = cos(model->psi);
            sp = sin(model->psi);

            model->rot[0][0] = +cp;
            model->rot[1][0] = -sp*ci;
            model->rot[2][0] = +sp*si;
            model->rot[0][1] = +sp;
            model->rot[1][1] = +cp*ci;
            model->rot[2][1] = -cp*si;
            model->rot[0][2] = 0.0;
            model->rot[1][2] = +si;
            model->rot[2][2] = +ci;
        }
        else if (lnum >= 21)
        {
            if (tindex < nterms)
            {
                if (tindex == 0 && var == 1 && !read_mean_longitude)
                {
                    /* Special case: there is an extra row of mean longitude values here. */
                    nscanned = sscanf(line, "%lf %lf",
                        &model->moon[mindex].al[0],
                        &model->moon[mindex].al[1]);

                    REQUIRE_NSCAN(2);
                    read_mean_longitude = 1;
                }
                else
                {
                    /* Load the next term for the current moon and variable. */
                    nscanned = sscanf(line, "%d %lf %lf %lf",
                        &check,
                        &series->term[tindex].amplitude,
                        &series->term[tindex].phase,
                        &series->term[tindex].frequency);

                    REQUIRE_NSCAN(4);

                    if (check != tindex + 1)
                        FAIL("LoadJupiterMoonModel(%s line %d): invalid check value %d (expected %d)\n", filename, lnum, check, tindex + 1);

                    ++tindex;
                }
            }
            else
            {
                read_mean_longitude = 0;
                tindex = 0;
                if (var == 3)
                {
                    /* Start the next moon. */
                    var = 0;
                    ++mindex;
                    if (mindex == NUM_JUPITER_MOONS)
                        break;  /* We are done loading data! (We don't care about the Tchebycheff polynomials.) */
                }
                else
                {
                    /* Start the next variable in the current moon. */
                    ++var;
                }

                /* Select which variable corresponds to the value of 'var'. */
                switch (var)
                {
                    case 0:  series = &model->moon[mindex].a;     break;
                    case 1:  series = &model->moon[mindex].l;     break;
                    case 2:  series = &model->moon[mindex].z;     break;
                    case 3:  series = &model->moon[mindex].zeta;  break;
                    default:
                        FAIL("LoadJupiterMoonModel: Invalid var = %d\n", var);
                }

                /* Read the number of terms in this variable. */
                nscanned = sscanf(line, "%d", &nterms);
                REQUIRE_NSCAN(1);
                if (nterms < 0 || nterms > MAX_JM_SERIES)
                    FAIL("LoadJupiterMoonModel(%s line %d): Invalid nterms = %d\n", filename, lnum, nterms);

                /* Allocate terms for this series from the term buffer. */
                series->nterms_calc = series->nterms_total = nterms;
                series->term = &model->buffer[next_term_index];
                next_term_index += nterms;
                if (next_term_index > MAX_JUPITER_TERMS)
                    FAIL("LoadJupiterMoonModel(%s line %d): Exhausted the term buffer!\n", filename, lnum);
            }
        }
    }

    if (mindex != NUM_JUPITER_MOONS)
        FAIL("LoadJupiterMoonModel(%s): only loaded %d moons.\n", filename, mindex);

    error = 0;
fail:
    if (infile != NULL) fclose(infile);
    return error;
}


static int JupiterMoonWriteVar(FILE *outfile, const vsop_series_t *series, const double *al)
{
    int tindex;

    fprintf(outfile, "%d\n", series->nterms_calc);

    if (al != NULL)
        fprintf(outfile, "%23.16le %23.16le\n", al[0], al[1]);

    for (tindex = 0; tindex < series->nterms_calc; ++tindex)
    {
        fprintf(outfile, "%3d %23.16lf %20.13le %20.13le\n",
            1 + tindex,
            series->term[tindex].amplitude,
            series->term[tindex].phase,
            series->term[tindex].frequency);
    }

    return 0;
}



int SaveJupiterMoonModel(const char *filename, const jupiter_moon_model_t *model)
{
    int error, mindex;
    FILE *outfile;

    outfile = fopen(filename, "wt");
    if (outfile == NULL)
        FAIL("SaveJupiterMoonModel: cannot open output file: %s\n", filename);

    /* As a hack, we emit 18 blank lines, because we ignore those from the original L1.2 file. */
    for (mindex = 0; mindex < 18; ++mindex)
        fprintf(outfile, "\n");

    /* Write G(M+m) constants */
    fprintf(outfile, "%23.16le %23.16le %23.16le %23.16le\n", model->moon[0].mu, model->moon[1].mu, model->moon[2].mu, model->moon[3].mu);

    /* Write Jupiter equatorial orientation angles. */
    fprintf(outfile, "%23.16le %23.16le\n", model->psi, model->incl);

    for (mindex = 0; mindex < NUM_JUPITER_MOONS; ++mindex)
    {
        CHECK(JupiterMoonWriteVar(outfile, &model->moon[mindex].a,    NULL));
        CHECK(JupiterMoonWriteVar(outfile, &model->moon[mindex].l,    model->moon[mindex].al));
        CHECK(JupiterMoonWriteVar(outfile, &model->moon[mindex].z,    NULL));
        CHECK(JupiterMoonWriteVar(outfile, &model->moon[mindex].zeta, NULL));
    }

    fprintf(outfile, "\n");     /* final blank line needed for loader to know it has all 4 moons. */

    error = 0;
fail:
    if (outfile != NULL) fclose(outfile);
    return error;
}


static int JupiterMoons_C(cg_context_t *context, const jupiter_moon_model_t *model)
{
    int mindex, var, i;
    const char *moon_name[] = { "Io", "Europa", "Ganymede", "Callisto" };
    const char *var_name[] = { "a", "l", "z", "zeta" };
    vsop_series_t series[4];

    fprintf(context->outfile, "static const astro_rotation_t Rotation_JUP_EQJ =\n");
    fprintf(context->outfile, "{\n");
    fprintf(context->outfile, "    ASTRO_SUCCESS,\n");
    fprintf(context->outfile, "    {\n");
    fprintf(context->outfile, "        { %23.16le, %23.16le, %23.16le },\n", model->rot[0][0], model->rot[0][1], model->rot[0][2]);
    fprintf(context->outfile, "        { %23.16le, %23.16le, %23.16le },\n", model->rot[1][0], model->rot[1][1], model->rot[1][2]);
    fprintf(context->outfile, "        { %23.16le, %23.16le, %23.16le }\n",  model->rot[2][0], model->rot[2][1], model->rot[2][2]);
    fprintf(context->outfile, "    }\n");
    fprintf(context->outfile, "};\n\n");

    for (mindex = 0; mindex < NUM_JUPITER_MOONS; ++mindex)
    {
        series[0] = model->moon[mindex].a;
        series[1] = model->moon[mindex].l;
        series[2] = model->moon[mindex].z;
        series[3] = model->moon[mindex].zeta;
        for (var = 0; var < NUM_JM_VARS; ++var)
        {
            int nterms = series[var].nterms_total;
            fprintf(context->outfile, "static const vsop_term_t jm_%s_%s[] =\n", moon_name[mindex], var_name[var]);
            fprintf(context->outfile, "{\n");
            for (i = 0; i < nterms; ++i)
            {
                const vsop_term_t *term = &series[var].term[i];
                fprintf(context->outfile, "    { %19.16lf, %23.16le, %23.16le }%s\n",
                    term->amplitude,
                    term->phase,
                    term->frequency,
                    (i+1 < nterms) ? "," : "");
            }
            fprintf(context->outfile, "};\n\n");
        }
    }

    fprintf(context->outfile, "static const jupiter_moon_t JupiterMoonModel[] =\n");
    fprintf(context->outfile, "{\n");
    for (mindex = 0; mindex < NUM_JUPITER_MOONS; ++mindex)
    {
        series[0] = model->moon[mindex].a;
        series[1] = model->moon[mindex].l;
        series[2] = model->moon[mindex].z;
        series[3] = model->moon[mindex].zeta;
        fprintf(context->outfile, "    { %23.16le, { %23.16le, %23.16le }", model->moon[mindex].mu, model->moon[mindex].al[0], model->moon[mindex].al[1]);
        for (var = 0; var < NUM_JM_VARS; ++var)
            fprintf(context->outfile, "\n        , { %3d, jm_%s_%-4s }", series[var].nterms_total, moon_name[mindex], var_name[var]);
        fprintf(context->outfile, " }%s\n", ((mindex + 1 < NUM_JUPITER_MOONS) ? ",\n" : ""));
    }
    fprintf(context->outfile, "}");
    return 0;
}


static int JupiterMoons_CSharp(cg_context_t *context, const jupiter_moon_model_t *model)
{
    int mindex, var, i;
    vsop_series_t series[4];
    const char *moon_name[] = { "Io", "Europa", "Ganymede", "Callisto" };
    const char *var_name[] = { "a", "l", "z", "zeta" };

    fprintf(context->outfile, "        private static readonly RotationMatrix Rotation_JUP_EQJ = new RotationMatrix(\n");
    fprintf(context->outfile, "            new double[3,3]\n");
    fprintf(context->outfile, "            {\n");
    fprintf(context->outfile, "                { %23.16le, %23.16le, %23.16le },\n", model->rot[0][0], model->rot[0][1], model->rot[0][2]);
    fprintf(context->outfile, "                { %23.16le, %23.16le, %23.16le },\n", model->rot[1][0], model->rot[1][1], model->rot[1][2]);
    fprintf(context->outfile, "                { %23.16le, %23.16le, %23.16le }\n",  model->rot[2][0], model->rot[2][1], model->rot[2][2]);
    fprintf(context->outfile, "            }\n");
    fprintf(context->outfile, "        );\n\n");

    fprintf(context->outfile, "        private static readonly jupiter_moon_t[] JupiterMoonModel = new jupiter_moon_t[] {\n");
    for (mindex = 0; mindex < NUM_JUPITER_MOONS; ++mindex)
    {
        series[0] = model->moon[mindex].a;
        series[1] = model->moon[mindex].l;
        series[2] = model->moon[mindex].z;
        series[3] = model->moon[mindex].zeta;
        fprintf(context->outfile, "            // [%d] %s\n", mindex, moon_name[mindex]);
        fprintf(context->outfile, "            new jupiter_moon_t {\n");
        fprintf(context->outfile, "                mu = %23.16le,\n", model->moon[mindex].mu);
        fprintf(context->outfile, "                al0 = %23.16le,\n", model->moon[mindex].al[0]);
        fprintf(context->outfile, "                al1 = %23.16le,\n", model->moon[mindex].al[1]);
        for (var = 0; var < NUM_JM_VARS; ++var)
        {
            int nterms = series[var].nterms_total;
            fprintf(context->outfile, "                %s = new vsop_term_t[] {\n", var_name[var]);
            for (i = 0; i < nterms; ++i)
            {
                const vsop_term_t *term = &series[var].term[i];
                fprintf(context->outfile, "                    new vsop_term_t(%19.16lf, %23.16le, %23.16le)%s\n",
                    term->amplitude,
                    term->phase,
                    term->frequency,
                    (i+1 < nterms) ? "," : "");
            }
            fprintf(context->outfile, "                }%s\n", ((var+1 < NUM_JM_VARS) ? "," : ""));
        }
        fprintf(context->outfile, "            }%s\n", ((mindex+1 < NUM_JUPITER_MOONS) ? ",\n" : ""));
    }
    fprintf(context->outfile, "        }");
    return 0;
}


static int JupiterMoons_JS(cg_context_t *context, const jupiter_moon_model_t *model)
{
    int mindex, var, i;
    vsop_series_t series[4];
    const char *moon_name[] = { "Io", "Europa", "Ganymede", "Callisto" };
    const char *var_name[] = { "a", "l", "z", "zeta" };

    fprintf(context->outfile, "const Rotation_JUP_EQJ = new RotationMatrix([\n");
    fprintf(context->outfile, "    [ %23.16le, %23.16le, %23.16le ],\n", model->rot[0][0], model->rot[0][1], model->rot[0][2]);
    fprintf(context->outfile, "    [ %23.16le, %23.16le, %23.16le ],\n", model->rot[1][0], model->rot[1][1], model->rot[1][2]);
    fprintf(context->outfile, "    [ %23.16le, %23.16le, %23.16le ]\n",  model->rot[2][0], model->rot[2][1], model->rot[2][2]);
    fprintf(context->outfile, "]);\n\n");

    fprintf(context->outfile, "const JupiterMoonModel: jupiter_moon_t[] = [\n");
    for (mindex = 0; mindex < NUM_JUPITER_MOONS; ++mindex)
    {
        series[0] = model->moon[mindex].a;
        series[1] = model->moon[mindex].l;
        series[2] = model->moon[mindex].z;
        series[3] = model->moon[mindex].zeta;
        fprintf(context->outfile, "    // [%d] %s\n", mindex, moon_name[mindex]);
        fprintf(context->outfile, "    {\n");
        fprintf(context->outfile, "        mu: %23.16le,\n", model->moon[mindex].mu);
        fprintf(context->outfile, "        al: [%23.16le, %23.16le],\n", model->moon[mindex].al[0], model->moon[mindex].al[1]);
        for (var = 0; var < NUM_JM_VARS; ++var)
        {
            int nterms = series[var].nterms_total;
            fprintf(context->outfile, "        %s: [\n", var_name[var]);
            for (i = 0; i < nterms; ++i)
            {
                const vsop_term_t *term = &series[var].term[i];
                fprintf(context->outfile, "            [ %19.16lf, %23.16le, %23.16le ]%s\n",
                    term->amplitude,
                    term->phase,
                    term->frequency,
                    (i+1 < nterms) ? "," : "");
            }
            fprintf(context->outfile, "        ]%s\n", ((var+1 < NUM_JM_VARS) ? "," : ""));
        }
        fprintf(context->outfile, "    }%s\n", ((mindex+1 < NUM_JUPITER_MOONS) ? ",\n" : ""));
    }
    fprintf(context->outfile, "]");
    return 0;
}


static int JupiterMoons_Python(cg_context_t *context, const jupiter_moon_model_t *model)
{
    int mindex, var, i;
    vsop_series_t series[4];
    const char *moon_name[] = { "Io", "Europa", "Ganymede", "Callisto" };
    const char *var_name[] = { "a", "l", "z", "zeta" };

    fprintf(context->outfile, "_Rotation_JUP_EQJ = RotationMatrix([\n");
    fprintf(context->outfile, "    [ %23.16le, %23.16le, %23.16le ],\n", model->rot[0][0], model->rot[0][1], model->rot[0][2]);
    fprintf(context->outfile, "    [ %23.16le, %23.16le, %23.16le ],\n", model->rot[1][0], model->rot[1][1], model->rot[1][2]);
    fprintf(context->outfile, "    [ %23.16le, %23.16le, %23.16le ]\n",  model->rot[2][0], model->rot[2][1], model->rot[2][2]);
    fprintf(context->outfile, "])\n\n");

    fprintf(context->outfile, "_JupiterMoonModel = [\n");
    for (mindex = 0; mindex < NUM_JUPITER_MOONS; ++mindex)
    {
        series[0] = model->moon[mindex].a;
        series[1] = model->moon[mindex].l;
        series[2] = model->moon[mindex].z;
        series[3] = model->moon[mindex].zeta;
        fprintf(context->outfile, "    # [%d] %s\n", mindex, moon_name[mindex]);
        fprintf(context->outfile, "    [\n");
        fprintf(context->outfile, "        %23.16le, %23.16le, %23.16le, # mu, al0, al1\n", model->moon[mindex].mu, model->moon[mindex].al[0], model->moon[mindex].al[1]);
        for (var = 0; var < NUM_JM_VARS; ++var)
        {
            int nterms = series[var].nterms_total;
            fprintf(context->outfile, "        [   # %s\n", var_name[var]);
            for (i = 0; i < nterms; ++i)
            {
                const vsop_term_t *term = &series[var].term[i];
                fprintf(context->outfile, "            [ %19.16lf, %23.16le, %23.16le ]%s\n",
                    term->amplitude,
                    term->phase,
                    term->frequency,
                    (i+1 < nterms) ? "," : "");
            }
            fprintf(context->outfile, "        ]%s\n", ((var+1 < NUM_JM_VARS) ? "," : ""));
        }
        fprintf(context->outfile, "    ]%s\n", ((mindex+1 < NUM_JUPITER_MOONS) ? ",\n" : ""));
    }
    fprintf(context->outfile, "]");
    return 0;
}


static int JupiterMoons(cg_context_t *context)
{
    int error;
    jupiter_moon_model_t *model;

    model = calloc(1, sizeof(jupiter_moon_model_t));
    if (model == NULL)
        FAIL("JupiterMoons: memory allocation failure!\n");

    CHECK(LoadJupiterMoonModel("output/jupiter_moons.txt", model));

    switch (context->language)
    {
    case CODEGEN_LANGUAGE_C:
        CHECK(JupiterMoons_C(context, model));
        break;

    case CODEGEN_LANGUAGE_CSHARP:
        CHECK(JupiterMoons_CSharp(context, model));
        break;

    case CODEGEN_LANGUAGE_JS:
        CHECK(JupiterMoons_JS(context, model));
        break;

    case CODEGEN_LANGUAGE_PYTHON:
        CHECK(JupiterMoons_Python(context, model));
        break;

    default:
        CHECK(LogError(context, "JupiterMoons: Unsupported target language %d", context->language));
    }

    error = 0;
fail:
    free(model);
    return error;
}


/*-------------------------- end Jupiter moons -----------------------------*/


static int LogError(const cg_context_t *context, const char *format, ...)
{
    va_list v;
    va_start(v, format);
    fprintf(stderr, "ERROR(%s %d): ", context->inFileName, context->lnum);
    vfprintf(stderr, format, v);
    fprintf(stderr, "\n");
    va_end(v);
    return 1;
}

static const cg_directive_entry DirectiveTable[] =
{
    { "C_VSOP",             CVsop               },
    { "CSHARP_VSOP",        CsharpVsop          },
    { "LIST_VSOP",          ListVsop            },
    { "IAU_DATA",           OptIauData          },
    { "ADDSOL",             OptAddSol           },
    { "CONSTEL",            ConstellationData   },
    { "C_PLUTO_CONST",      PlutoConstants_C    },
    { "PLUTO_TABLE",        PlutoStateTable     },
    { "JUPITER_MOONS",      JupiterMoons        },
    { NULL, NULL }  /* Marks end of list */
};

static int ProcessDirective(cg_context_t *context)
{
    int i;
    for (i = 0; DirectiveTable[i].verb != NULL; ++i)
        if (!strcmp(context->verb, DirectiveTable[i].verb))
            return DirectiveTable[i].func(context);

    return LogError(context, "Unknown verb '%s'", context->verb);
}
