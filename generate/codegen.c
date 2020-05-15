/*
    MIT License

    Copyright (c) 2019-2020 Don Cross <cosinekitty@gmail.com>

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

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "novas.h"
#include "codegen.h"
#include "vsop.h"
#include "ephfile.h"

#define EXTRAPOLATE_DT

#define CG_MAX_LINE_LENGTH  200
#define MAX_DATA_PER_LINE    20
#define IAU_DATA_PER_ROW     11
#define IAU_DATA_NUM_ROWS    77
#define ADDSOL_DATA_PER_ROW   8

static const double MJD_BASIS = 2400000.5;

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

static int ListChebyshev(cg_context_t *context)
{
    int error = 1;
    int body, i, record_index;
    char filename[100];
    eph_file_reader_t reader;
    eph_record_t record;

    if (1 != sscanf(context->args, "%d", &body) || body < 0 || body > 8)
    {
        error = LogError(context, "Chebyshev body name is invalid.");
        goto fail;
    }

    snprintf(filename, sizeof(filename), "output/%02d.eph", body);
    error = EphFileOpen(&reader, filename);
    if (error)
    {
        LogError(context, "EphFileOpen returned error %d for file: %s", error, filename);
        goto fail;
    }

    fprintf(context->outfile, "[\n");
    for (record_index=0; EphReadRecord(&reader, &record); ++record_index)
    {
        if (record_index > 0)
            fprintf(context->outfile, ",\n");

        fprintf(context->outfile, "{ 'tt':%lf, 'ndays':%lf, 'coeff':[\n", record.jdStart - T0, record.jdDelta);
        for (i=0; i < record.numpoly; ++i)
        {
            fprintf(context->outfile, "    [%0.12lf, %0.12lf, %0.12lf]%s\n",
                record.coeff[0][i],
                record.coeff[1][i],
                record.coeff[2][i],
                (i+1 < record.numpoly) ? "," : "]");
        }
        fprintf(context->outfile, "}");
    }
    fprintf(context->outfile, "]");

    if (record.error)
    {
        LogError(context, "Error %d in EphReadRecord for line %d in file %s", record.error, reader.lnum, filename);
        error = record.error;
        goto fail;
    }

fail:
    EphFileClose(&reader);
    return error;
}

static int CChebyshev(cg_context_t *context)
{
    int error = 1;
    int body, i, record_index, record_count;
    char filename[100];
    eph_file_reader_t reader;
    eph_record_t record;

    if (1 != sscanf(context->args, "%d", &body) || body < 0 || body > 8)
    {
        error = LogError(context, "Chebyshev body name is invalid.");
        goto fail;
    }

    snprintf(filename, sizeof(filename), "output/%02d.eph", body);
    error = EphFileOpen(&reader, filename);
    if (error)
    {
        LogError(context, "EphFileOpen #1 returned error %d for file: %s", error, filename);
        goto fail;
    }

    for (record_index=0; EphReadRecord(&reader, &record); ++record_index)
    {
        fprintf(context->outfile, "static const astro_cheb_coeff_t cheb_%d_%d[] =\n", body, record_index);
        fprintf(context->outfile, "{\n");
        for (i=0; i < record.numpoly; ++i)
        {
            fprintf(context->outfile, "    { { %16.12lf, %16.12lf, %16.12lf } }%s\n",
                record.coeff[0][i],
                record.coeff[1][i],
                record.coeff[2][i],
                (i+1 < record.numpoly) ? "," : "");
        }
        fprintf(context->outfile, "};\n\n");
    }
    record_count = record_index;

    if (record.error)
    {
        LogError(context, "Error %d in EphReadRecord#1 for line %d in file %s", record.error, reader.lnum, filename);
        error = record.error;
        goto fail;
    }

    EphFileClose(&reader);
    error = EphFileOpen(&reader, filename);
    if (error)
    {
        LogError(context, "EphFileOpen #2 returned error %d for file: %s", error, filename);
        goto fail;
    }

    fprintf(context->outfile, "static const astro_cheb_record_t cheb_%d[] =\n{\n", body);
    for (record_index=0; EphReadRecord(&reader, &record); ++record_index)
    {
        fprintf(context->outfile, "    { %10.1lf, %7.1lf, ARRAYSIZE(cheb_%d_%d), cheb_%d_%d }%s\n",
            record.jdStart - T0,
            record.jdDelta,
            body,
            record_index,
            body,
            record_index,
            (record_index+1 < record_count) ? "," : "");
    }
    fprintf(context->outfile, "}");

    if (record.error)
    {
        LogError(context, "Error %d in EphReadRecord#2 for line %d in file %s", record.error, reader.lnum, filename);
        error = record.error;
        goto fail;
    }

fail:
    EphFileClose(&reader);
    return error;
}

static int CsharpChebyshev(cg_context_t *context)
{
    int error = 1;
    int body, i, record_index, record_count;
    char filename[100];
    eph_file_reader_t reader;
    eph_record_t record;

    if (1 != sscanf(context->args, "%d", &body) || body < 0 || body > 8)
    {
        error = LogError(context, "Chebyshev body name is invalid.");
        goto fail;
    }

    snprintf(filename, sizeof(filename), "output/%02d.eph", body);
    error = EphFileOpen(&reader, filename);
    if (error)
    {
        LogError(context, "EphFileOpen #1 returned error %d for file: %s", error, filename);
        goto fail;
    }

    for (record_index=0; EphReadRecord(&reader, &record); ++record_index)
    {
        fprintf(context->outfile, "        private static readonly astro_cheb_coeff_t[] cheb_%d_%d =\n", body, record_index);
        fprintf(context->outfile, "        {\n");
        for (i=0; i < record.numpoly; ++i)
        {
            fprintf(context->outfile, "            new astro_cheb_coeff_t(%16.12lf, %16.12lf, %16.12lf)%s\n",
                record.coeff[0][i],
                record.coeff[1][i],
                record.coeff[2][i],
                (i+1 < record.numpoly) ? "," : "");
        }
        fprintf(context->outfile, "        };\n\n");
    }
    record_count = record_index;

    if (record.error)
    {
        LogError(context, "Error %d in EphReadRecord#1 for line %d in file %s", record.error, reader.lnum, filename);
        error = record.error;
        goto fail;
    }

    EphFileClose(&reader);
    error = EphFileOpen(&reader, filename);
    if (error)
    {
        LogError(context, "EphFileOpen #2 returned error %d for file: %s", error, filename);
        goto fail;
    }

    fprintf(context->outfile, "        private static readonly astro_cheb_record_t[] cheb_%d =\n        {\n", body);
    for (record_index=0; EphReadRecord(&reader, &record); ++record_index)
    {
        fprintf(context->outfile, "            new astro_cheb_record_t(%10.1lf, %7.1lf, cheb_%d_%d)%s\n",
            record.jdStart - T0,
            record.jdDelta,
            body,
            record_index,
            (record_index+1 < record_count) ? "," : "");
    }
    fprintf(context->outfile, "        }");

    if (record.error)
    {
        LogError(context, "Error %d in EphReadRecord#2 for line %d in file %s", record.error, reader.lnum, filename);
        error = record.error;
        goto fail;
    }

fail:
    EphFileClose(&reader);
    return error;
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
    char sname[100];

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
    const char *coord_name[3] = { "lat", "lon", "rad" };

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
    char sname[100];

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
    const char *coord_name[3] = { "lat", "lon", "rad" };

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

static int GenArrayEnd(cg_context_t *context)
{
    switch (context->language)
    {
    case CODEGEN_LANGUAGE_JS:
        fprintf(context->outfile, "\n]");
        return 0;

    case CODEGEN_LANGUAGE_C:
        fprintf(context->outfile, "\n}");
        return 0;

    case CODEGEN_LANGUAGE_CSHARP:
        fprintf(context->outfile, "\n        }");
        return 0;

    case CODEGEN_LANGUAGE_PYTHON:
        fprintf(context->outfile, "\n]");
        return 0;

    default:
        return LogError(context, "GenArrayEnd: Unknown language type %d", context->language);
    }
}

static int GenDeltaTArrayEntry(cg_context_t *context, int count, double mjd, const char *dt_text)
{
    switch (context->language)
    {
    case CODEGEN_LANGUAGE_C:
        fprintf(context->outfile, "%s\n", (count==1) ? "{" : ",");
        fprintf(context->outfile, "{ %0.1lf, %s }", mjd, dt_text);
        return 0;

    case CODEGEN_LANGUAGE_CSHARP:
        fprintf(context->outfile, "%s\n", (count==1) ? "new deltat_entry_t[]\n        {" : ",");
        fprintf(context->outfile, "            new deltat_entry_t { mjd=%0.1lf, dt=%s }", mjd, dt_text);
        return 0;

    case CODEGEN_LANGUAGE_JS:
        fprintf(context->outfile, "%s\n", (count==1) ? "[" : ",");
        fprintf(context->outfile, "[%7lg, %s]", mjd, dt_text);
        return 0;

    case CODEGEN_LANGUAGE_PYTHON:
        fprintf(context->outfile, "%s\n", (count==1) ? "[" : ",");
        fprintf(context->outfile, "_delta_t_entry_t(%0.1lf, %s)", mjd, dt_text);
        return 0;

    default:
        return LogError(context, "GenDeltaTArrayEntry: Unknown language type %d", context->language);
    }
}

#define MJD_FROM_DATE(year, month, day)     (julian_date((short)(year), (short)(month), (short)(day), 0.0) - MJD_BASIS)
#define MJD_FROM_YEAR(year)                 MJD_FROM_DATE(year, 1, 1)

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

#ifdef EXTRAPOLATE_DT

static int GenExtrapolatedDeltaT(cg_context_t *context, int year, int count)
{
    double  mjd, dt;
    char    dt_text[20];

    mjd = MJD_FROM_YEAR(year);
    dt = ExtrapolatedDeltaT(year);
    snprintf(dt_text, sizeof(dt_text), "%0.2lf", dt);
    return GenDeltaTArrayEntry(context, count, mjd, dt_text);
}

#endif /* EXTRAPOLATE_DT */

static int GenDeltaT(cg_context_t *context)
{
    FILE *infile;
    int error=1, lnum, count=0;
    const char *filename;
    char line[100];
    char dt_text[20];
    int year, frac_year, month, day;
    double dt, float_year;
    double mjd = 0.0;
    double last_mjd;

    filename = "delta_t/historic.txt";
    infile = fopen(filename, "rt");
    if (infile == NULL) goto fail;
    lnum = 0;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;
        if (lnum < 3) continue;     /* first 2 lines are headers */
        if (3 != sscanf(line, "%d.%d %20s", &year, &frac_year, dt_text) || 1 != sscanf(dt_text, "%lf", &dt))
        {
            error = LogError(context, "Line %d of file %s has invalid format.\n", lnum, filename);
            goto fail;
        }

        if (frac_year == 0)
        {
            /* reduce the data size */
            if (year < 1750 && year % 20 != 0) continue;
            if (year < 1850 && year % 10 != 0) continue;
            if (year % 5 != 0) continue;

            mjd = MJD_FROM_YEAR(year);
            ++count;
            CHECK(GenDeltaTArrayEntry(context, count, mjd, dt_text));
        }
    }
    fclose(infile);

    filename = "delta_t/recent.txt";
    infile = fopen(filename, "rt");
    if (infile == NULL) goto fail;
    lnum = 0;
    last_mjd = mjd;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;
        if (4 != sscanf(line, "%d %d %d %20s", &year, &month, &day, dt_text) || 1 != sscanf(dt_text, "%lf", &dt))
        {
            error = LogError(context, "Line %d of file %s has invalid format.", lnum, filename);
            goto fail;
        }

        /* reduce the data size by keeping only 1 sample per year */
        if (month != 1) continue;

        mjd = MJD_FROM_DATE(year, month, day);
        if (mjd > last_mjd)
        {
            ++count;
            CHECK(GenDeltaTArrayEntry(context, count, mjd, dt_text));
        }
    }
    fclose(infile);

    filename = "delta_t/predicted.txt";
    infile = fopen(filename, "rt");
    if (infile == NULL) goto fail;
    last_mjd = mjd;
    lnum = 0;
    float_year = 0.0;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;
        if (lnum < 2) continue;     /* skip header line */
        if (3 != sscanf(line, "%lf %lf %20s", &mjd, &float_year, dt_text) || 1 != sscanf(dt_text, "%lf", &dt))
        {
            error = LogError(context, "Line %d of file %s has invalid format.", lnum, filename);
            goto fail;
        }

        /* reduce the data size by keeping only 1 sample per year */
        if (float_year != floor(float_year)) continue;

        if (mjd > last_mjd)
        {
            ++count;
            CHECK(GenDeltaTArrayEntry(context, count, mjd, dt_text));
        }
    }

    /* Keep the final data point to maximize the extent of predictions into the future. */
    if (mjd > last_mjd && float_year != floor(float_year))
    {
        ++count;
        CHECK(GenDeltaTArrayEntry(context, count, mjd, dt_text));
    }

#ifdef EXTRAPOLATE_DT
    /* Add an extrapolated delta_t value for every decade 2030..2200. */
    for (year = 2030; year <= 2200; year += 10)
    {
        ++count;
        CHECK(GenExtrapolatedDeltaT(context, year, count));
    }
#endif

    CHECK(GenArrayEnd(context));

    if (count < 2)
    {
        error = LogError(context, "There must be at least 2 delta_t data!");
        goto fail;
    }

    error = 0;

fail:
    if (infile == NULL)
        error = LogError(context, "Cannot open input file: %s", filename);
    else
        fclose(infile);
    return error;
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

    if (nprinted >= size)
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

    fprintf(context->outfile, "    { nals:[ %2.0lf", data[0]);

    for (i=1; i < 5; ++i)
        fprintf(context->outfile, ", %2.0lf", data[i]);

    fprintf(context->outfile, " ], cls:[ %12.0lf", data[i]);
    for (++i; i < 11; ++i)
        fprintf(context->outfile, ", %12.0lf", data[i]);

    fprintf(context->outfile, " ] }");
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

        ra_lo /= 15.0;      /* convert RA from degrees to sidereal hours */
        ra_hi /= 15.0;      /* convert RA from degrees to sidereal hours */

        /* Search for the constellation symbol in the name table, to find its array index. */
        for (index=0; index < NUM_CONSTELLATIONS; ++index)
            if (0 == strcmp(symbol, dict[index]))
                break;

        if (index == NUM_CONSTELLATIONS)
            CHECK(LogError(context, "Invalid constellation symbol '%s' at %s line %d", symbol, borderFileName, lnum));

        switch (context->language)
        {
        case CODEGEN_LANGUAGE_C:
            fprintf(context->outfile, "%c   { %2d, %17.14lf, %17.14lf, %17.14lf }    /* %s */\n",
                ((lnum == 1) ? ' ' : ','),
                index, ra_lo, ra_hi, dec, symbol);
            break;

        case CODEGEN_LANGUAGE_CSHARP:
            fprintf(context->outfile, "        %c   new constel_boundary_t(%2d, %17.14lf, %17.14lf, %17.14lf)    // %s\n",
                ((lnum == 1) ? ' ' : ','),
                index, ra_lo, ra_hi, dec, symbol);
            break;

        case CODEGEN_LANGUAGE_JS:
            fprintf(context->outfile, "%c   [ %2d, %6lg, %6lg, %6lg ]    // %s\n",
                ((lnum == 1) ? ' ' : ','),
                index,
                RoundAngle(15.0 * ra_lo),
                RoundAngle(15.0 * ra_hi),
                RoundAngle(dec),
                symbol);
            break;

        case CODEGEN_LANGUAGE_PYTHON:
            fprintf(context->outfile, "%c   ( %2d, %17.14lf, %17.14lf, %17.14lf )    # %s\n",
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
    { "LIST_CHEBYSHEV",     ListChebyshev       },
    { "C_CHEBYSHEV",        CChebyshev          },
    { "CSHARP_CHEBYSHEV",   CsharpChebyshev     },
    { "DELTA_T",            GenDeltaT           },
    { "IAU_DATA",           OptIauData          },
    { "ADDSOL",             OptAddSol           },
    { "CONSTEL",            ConstellationData   },
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
