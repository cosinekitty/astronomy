/*
    MIT License

    Copyright (c) 2019-2024 Don Cross <cosinekitty@gmail.com>

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
#ifndef __DDC_ASTRO_CODEGEN_H
#define __DDC_ASTRO_CODEGEN_H

#include "vsop.h"

#define CHECK(x)    do{if(0 != (error = (x))) goto fail;}while(0)
#define FAIL(...)   do{fprintf(stderr, __VA_ARGS__); error = 1; goto fail;}while(0)

typedef enum
{
    CODEGEN_LANGUAGE_UNKNOWN,
    CODEGEN_LANGUAGE_JS,
    CODEGEN_LANGUAGE_C,
    CODEGEN_LANGUAGE_CSHARP,
    CODEGEN_LANGUAGE_PYTHON,
    CODEGEN_LANGUAGE_KOTLIN
}
cg_language_t;

int GenerateCode(
    cg_language_t language,
    const char *outCodeFileName,
    const char *inTemplateFileName,
    const char *dataPath);

double ExtrapolatedDeltaT(int year);

/*-------------------------- Jupiter moons -----------------------------*/

#define NUM_JUPITER_MOONS     4
#define MAX_JM_SERIES        50
#define NUM_JM_VARS           4
#define MAX_JUPITER_TERMS   (NUM_JUPITER_MOONS * NUM_JM_VARS * MAX_JM_SERIES)

typedef struct
{
    double          mu;          /* mu = G(M+m), where M = Jupiter mass, m = moon mass. */
    double          al[2];       /* mean longitude coefficients */
    vsop_series_t   a;
    vsop_series_t   l;
    vsop_series_t   z;
    vsop_series_t   zeta;
}
jupiter_moon_t;


typedef struct
{
    /* angles that rotate Jupiter equatorial system to EQJ */
    double          incl;
    double          psi;

    /* a rotation matrix calculated from incl and psi, for the convenience of the code generators. */
    double          rot[3][3];

    jupiter_moon_t  moon[NUM_JUPITER_MOONS];
    vsop_term_t     buffer[MAX_JUPITER_TERMS];
}
jupiter_moon_model_t;


int LoadJupiterMoonModel(const char *filename, jupiter_moon_model_t *model);
int SaveJupiterMoonModel(const char *filename, const jupiter_moon_model_t *model);

#endif /* __DDC_ASTRO_CODEGEN_H */
