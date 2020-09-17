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

#include <math.h>
#include <string.h>
#include "chebyshev.h"

/*
    Initialize a Chebyshev encoder.
    Prepares it for a series of calls to ChebApprox.
*/
int ChebInit(
    ChebEncoder *encoder,           /* the encoder to be initialized */
    int numpoly,                    /* the number of polynomials to be calculated: 1..CHEB_MAX_POLYS */
    int dimension)                  /* how many dimensions are in the vector returned by the function to be approximated: 1..CHEB_MAX_DIM */
{
    int j, k;
    const double pi = 3.1415926535897932384626433832795;

    memset(encoder, 0, sizeof(ChebEncoder));

    if (numpoly < 1 || numpoly > CHEB_MAX_POLYS)
        return 1;   /* invalid number of polynomials */

    if (dimension < 1 || dimension > CHEB_MAX_DIM)
        return 2;   /* invalid number of dimensions */

    encoder->NumPoly = numpoly;
    encoder->Dimension = dimension;

    for (j=0; j < numpoly; ++j)
        for (k=0; k < numpoly; ++k)
            encoder->Alpha[j][k] = cos((pi * j * (k + 0.5)) / numpoly);

    return 0;
}

/*
    Calculate Chebyshev polynomial coefficients that approximate
    the supplied function 'func'.
*/
int ChebGenerate(
    const ChebEncoder *encoder,     /* encoder that has already been initialized by a call to ChebInit() */
    ChebFuncType func,              /* the function to be sampled and later approximated */
    const void *context,            /* any function-specific data needed to evaluate the function */
    double t_min,                   /* the minimum value of the independent variable to be fed to func */
    double t_max,                   /* the maximum value of the independent variable to be fed to func */
    double coeff[CHEB_MAX_DIM][CHEB_MAX_POLYS])    /* the array[encoder->Dimension][encoder->NumPoly] of output coefficients */
{
    int j, k, d, error;
    double f[CHEB_MAX_POLYS][CHEB_MAX_DIM];
    const int n = encoder->NumPoly;
    const double t_center = (t_max + t_min) / 2.0;
    const double t_diff = (t_max - t_min) / 2.0;

    /* Pre-evaluate all the function values we will need. */
    for (k=0; k < n; ++k)
    {
        /* Cache the function value in f[k] so we can re-use it below, j times. */
        error = func(context, (t_diff * encoder->Alpha[1][k]) + t_center, f[k]);
        if (error)
            return error;
    }

    for (d=0; d < encoder->Dimension; ++d)
    {
        for (j=0; j < n; ++j)
        {
            double sum = 0.0;
            for (k=0; k < n; ++k)
                sum += encoder->Alpha[j][k] * f[k][d];

            coeff[d][j] = (2.0 / n) * sum;
        }
    }

    return 0;
}

/*
    Use Chebyshev polynomials to approximate the value of a function.
*/
void ChebApprox(
    int numpoly,                    /* the number of polynomial coefficients: 1..CHEB_MAX_POLYS */
    int dimension,                  /* the number of dimensions: 1..CHEB_MAX_DIM */
    const double coeff[CHEB_MAX_DIM][CHEB_MAX_POLYS],        /* array of coefficients */
    double x,                       /* dimensionless independent variable in the range [-1, +1] */
    double result[CHEB_MAX_DIM])    /* output approximation vector */
{
    int k, d;
    double sum, p0, p1, p2;

    for (d=0; d < dimension; ++d)
    {
        sum = 0.0;
        if (0 < numpoly)
        {
            p0 = 1.0;
            sum += coeff[d][0];
            if (1 < numpoly)
            {
                p1 = x;
                sum += coeff[d][1] * p1;
                for (k=2; k < numpoly; ++k)
                {
                    p2 = (2.0 * x * p1) - p0;
                    sum += coeff[d][k] * p2;
                    p0 = p1;
                    p1 = p2;
                }
            }
            sum -= coeff[d][0] / 2.0;
        }
        result[d] = sum;
    }
}

/*
    A utility function to scale a dimensional variable 't' in the problem domain
    (for example, a time offset expressed in seconds)
    to a dimensionless variable 'x' in the range [-1, +1] as required by ChebApprox().
*/
double ChebScale(
    double t_min,   /* the lower bound of the domain-specific range */
    double t_max,   /* the upper bound of the domain-specific range */
    double t)       /* a dimensional variable in the range [t_min, t_max]. */
{
    return (2*t - (t_max + t_min)) / (t_max - t_min);
}
