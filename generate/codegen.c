/*
    MIT License

    Copyright (c) 2019 Don Cross <cosinekitty@gmail.com>

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

#define CG_MAX_LINE_LENGTH  200

static const double MJD_BASIS = 2400000.5;

typedef struct
{
    FILE *outfile;
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

static int JsChebyshev(cg_context_t *context)
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
            
        fprintf(context->outfile, "{ tt:%lf, ndays:%lf, coeff:[\n", record.jdStart - T0, record.jdDelta);
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

static int JsVsop(cg_context_t *context)
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

static int JsDeltaT(cg_context_t *context)
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

            mjd = julian_date((short)year, 1, 1, 0.0) - MJD_BASIS;
            ++count;            
            fprintf(context->outfile, "%s\n", (count==1) ? "[" : ",");
            fprintf(context->outfile, "{ mjd:%0.1lf, dt:%s }", mjd, dt_text);
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

        mjd = julian_date((short)year, (short)month, (short)day, 0.0) - MJD_BASIS;
        if (mjd > last_mjd)
        {
            ++count;
            fprintf(context->outfile, "%s\n", (count==1) ? "[" : ",");
            fprintf(context->outfile, "{ mjd:%0.1lf, dt:%s }", mjd, dt_text);
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
            fprintf(context->outfile, "%s\n", (count==1) ? "[" : ",");
            fprintf(context->outfile, "{ mjd:%0.1lf, dt:%s }", mjd, dt_text);
        }
    }

    /* Keep the final data point to maximize the extent of predictions into the future. */
    if (mjd > last_mjd && float_year != floor(float_year))
    {
        ++count;
        fprintf(context->outfile, "%s\n", (count==1) ? "[" : ",");
        fprintf(context->outfile, "{ mjd:%0.1lf, dt:%s }", mjd, dt_text);
    }
    fprintf(context->outfile, "\n]");

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
    { "JS_VSOP", JsVsop },
    { "JS_CHEBYSHEV", JsChebyshev },
    { "JS_DELTA_T", JsDeltaT },
    { NULL, NULL }
};

static int ProcessDirective(cg_context_t *context)
{
    int i;
    for (i = 0; DirectiveTable[i].verb != NULL; ++i)
        if (!strcmp(context->verb, DirectiveTable[i].verb))
            return DirectiveTable[i].func(context);

    return LogError(context, "Unknown verb '%s'", context->verb);
}
