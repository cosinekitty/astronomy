/*
    MIT License

    Copyright (c) 2019-2023 Don Cross <cosinekitty@gmail.com>

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "astro_vector.h"

int VectorListAppend(VectorListType *list, double t, const double v[3])
{
    VectorType *nv;

    /* Make the array bigger as needed. */
    if (list->length == list->size)
    {
        int newsize = (list->size > 0) ? (2 * list->size) : 1024;
        VectorType *bigger = calloc(newsize, sizeof(VectorType));
        if (bigger == NULL)
        {
            fprintf(stderr, "VectorListAppend: out of memory!\n");
            return 1;
        }
        memcpy(bigger, list->array, list->length * sizeof(VectorType));
        free(list->array);
        list->array = bigger;
        list->size = newsize;
    }

    nv = &list->array[(list->length)++];
    nv->t = t;
    nv->v[0] = v[0];
    nv->v[1] = v[1];
    nv->v[2] = v[2];
    return 0;
}

void ClearVectorList(VectorListType *list)
{
    list->length = 0;
}

void FreeVectorList(VectorListType *list)
{
    free(list->array);
    memset(list, 0, sizeof(VectorListType));
}

typedef struct
{
    double sx;          /* sum(x) */
    double sxt;         /* sum(x*t) */
    double sxt2;        /* sum(x*t*t) */
}
ComponentSolver;      /* holds stuff needed to regress one of the vector columns (x, y, or z) */

typedef struct
{
    int N;              /* sum(1) = number of data points */
    double st1;         /* sum(t) */
    double st2;         /* sum(t^2) */
    double st3;         /* sum(t^3) */
    double st4;         /* sum(t^4) */
    ComponentSolver comp[3];        /* one component solver per (x, y, z) */
    double M[3][3];
}
VectorSolver;

static void Splice(double D[3][3], const double M[3][3], int col, const ComponentSolver *comp)
{
    memcpy(D, M, 9*sizeof(double));
    D[0][col] = comp->sxt2;
    D[1][col] = comp->sxt;
    D[2][col] = comp->sx;
}

int QuadraticRegression(const VectorListType *list, ParabolaType parab[3])
{
    int i, k;
    VectorSolver s;
    double A[3][3], B[3][3], C[3][3];
    double detM;

    memset(&s, 0, sizeof(s));
    s.N = list->length;
    for (i=0; i < list->length; ++i)
    {
        double t = list->array[i].t;
        double t2 = t*t;
        s.st1 += t;
        s.st2 += t2;
        s.st3 += t2*t;
        s.st4 += t2*t2;

        for (k=0; k<3; ++k)
        {
            double x = list->array[i].v[k];
            s.comp[k].sx += x;
            s.comp[k].sxt += x*t;
            s.comp[k].sxt2 += x*t2;
        }
    }

    s.M[0][0] = s.st4;  s.M[0][1] = s.st3;  s.M[0][2] = s.st2;
    s.M[1][0] = s.st3;  s.M[1][1] = s.st2;  s.M[1][2] = s.st1;
    s.M[2][0] = s.st2;  s.M[2][1] = s.st1;  s.M[2][2] = s.N;

    detM = Det(s.M);
    if (fabs(detM) < 1.0e-12)
    {
        fprintf(stderr, "QuadraticRegression: Cannot solve system of equations. det(M) = %lg, N = %d\n", detM, s.N);
        return 1;
    }

    for (k=0; k<3; ++k)
    {
        Splice(A, s.M, 0, &s.comp[k]);
        Splice(B, s.M, 1, &s.comp[k]);
        Splice(C, s.M, 2, &s.comp[k]);
        parab[k].a = Det(A) / detM;
        parab[k].b = Det(B) / detM;
        parab[k].c = Det(C) / detM;
    }

    return 0;
}

void QuadraticEval(const ParabolaType parab[3], double t, VectorType *approx)
{
    int k;
    approx->t = t;
    for (k=0; k<3; ++k)
        approx->v[k] = parab[k].a*t*t + parab[k].b*t + parab[k].c;
}

double Det(const double a[3][3])     /* determinant of a 3x3 matrix */
{
    return a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] 
         - a[0][2]*a[1][1]*a[2][0] - a[0][0]*a[1][2]*a[2][1] - a[0][1]*a[1][0]*a[2][2];
}

int TestRegression(const VectorListType *list, const ParabolaType parab[3], char tag)
{
    int i, k;
    const VectorType *exact;
    double sum, rms;
    VectorType approx;

    if (list->length <= 0)
    {
        fprintf(stderr, "TestRegression: Invalid list length %d\n", list->length);
        return 1;
    }

    sum = 0.0;
    for (i=0; i < list->length; ++i)
    {
        exact = &list->array[i];
        QuadraticEval(parab, exact->t, &approx);
        for (k=0; k < 3; ++k)
        {
            double diff = approx.v[k] - exact->v[k];
            sum += diff * diff;
        }
    }

    rms = sqrt(sum) / list->length;
    fprintf(stderr, "TestRegression(%c): N=%d, RMS=%lg\n", tag, list->length, rms);
    printf("# %c N=%d RMS=%lg\n", tag, list->length, rms);
    return 0;
}

VectorType Add(VectorType a, VectorType b)
{
    VectorType c;
    c.t = 0.0;      /* we lose time information */
    c.v[0] = a.v[0] + b.v[0];
    c.v[1] = a.v[1] + b.v[1];
    c.v[2] = a.v[2] + b.v[2];
    return c;
}

VectorType Subtract(VectorType a, VectorType b)
{
    VectorType c;
    c.t = 0.0;      /* we lose time information */
    c.v[0] = a.v[0] - b.v[0];
    c.v[1] = a.v[1] - b.v[1];
    c.v[2] = a.v[2] - b.v[2];
    return c;
}

double Length(VectorType a)
{
    return sqrt(a.v[0]*a.v[0] + a.v[1]*a.v[1] + a.v[2]*a.v[2]);
}

double Dot(VectorType a, VectorType b)
{
    return a.v[0]*b.v[0] + a.v[1]*b.v[1] + a.v[2]*b.v[2];
}

double SafeAcos(double x)
{
    double a;

    if (x < -1.0)
    {
        if (x < -1.000000001)
        {
            fprintf(stderr, "SafeAcos(FATAL): Argument too small (%lg)\n", x);
            exit(1);
        }
        return 3.14159265358979323846;
    }

    if (x > +1.0)
    {
        if (x > +1.000000001)
        {
            fprintf(stderr, "SafeAcos(FATAL): Argument too large (%lg)\n", x);
            exit(1);
        }
        return 0.0;
    }

    a = acos(x);
    if (isnan(a))
    {
        fprintf(stderr, "SafeAcos(FATAL): acos(%lg) = NAN\n", x);
        exit(1);
    }
    return a;
}

double AngleBetween(VectorType a, VectorType b)
{
    return SafeAcos(Dot(UnitVector(a), UnitVector(b)));
}

VectorType UnitVector(VectorType a)
{
    double mag = Length(a);
    if (fabs(mag) < 1.0e-12)
    {
        fprintf(stderr, "UnitVector: vector is too short!\n");
        exit(1);    /* not good -- leaks memory, leaves files unflushed, etc. */
    }
    a.v[0] /= mag;
    a.v[1] /= mag;
    a.v[2] /= mag;
    return a;
}
