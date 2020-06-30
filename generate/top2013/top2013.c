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

    TopInitModel(model);
    infile = fopen(filename, "rt");
    if (infile == NULL)
        FAIL("TopLoadModel: cannot open file: %s\n", filename);

    error = 0;
fail:
    if (infile != NULL) fclose(infile);
    return error;
}
