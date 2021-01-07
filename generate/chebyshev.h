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
#ifndef __DDC_CHEBYSHEV_H
#define __DDC_CHEBYSHEV_H

#define CHEB_MAX_POLYS  50
#define CHEB_MAX_DIM     3

/*
    In order to sample an arbitrary mathematical function of one real variable,
    we function pointer type that accepts an arbitrary context along with
    the independent variable.
    The context is used to hold whatever parameters are needed to evaluate
    the function.  In the case of ephemeris resampling, the context will
    include the solar system body identifier.
    The function stores the N-dimensional vector result of its calculation
    at the array indicated by the third parameter 'double[CHEB_MAX_DIM]'.
    The function returns 0 on success, nonzero on error.
    This is needed because ephemeris functions can return errors.
*/
typedef int (* ChebFuncType) (const void *, double, double[CHEB_MAX_DIM]);

typedef struct sChebEncoder
{
    int NumPoly;    /* number of Chebyshev polynomials to be calculated; 1..CHEB_MAX_POLYS */
    int Dimension;  /* how many dimensions are in the vector returned by the function to be approximated: 1..CHEB_MAX_DIM */
    double Alpha [CHEB_MAX_POLYS] [CHEB_MAX_POLYS];       /* Alpha[j][k] = cos((pi*j*(k + 1/2))/N) */
}
ChebEncoder;

/*
    Initialize a Chebyshev encoder.
    Prepares it for a series of calls to ChebApprox.
*/
int ChebInit(
    ChebEncoder *encoder,           /* the encoder to be initialized */
    int numpoly,                    /* the number of polynomials to be calculated: 1..CHEB_MAX_POLYS */
    int dimension);                 /* how many dimensions are in the vector returned by the function to be approximated: 1..CHEB_MAX_DIM */

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
    double coeff[CHEB_MAX_DIM][CHEB_MAX_POLYS]);                /* the array[encoder->Dimension][encoder->NumPoly] of output coefficients */

/*
    Use Chebyshev polynomials to approximate the value of a function.
*/
void ChebApprox(
    int numpoly,                    /* the number of polynomial coefficients: 1..CHEB_MAX_POLYS */
    int dimension,                  /* the number of dimensions: 1..CHEB_MAX_DIM */
    const double coeff[CHEB_MAX_DIM][CHEB_MAX_POLYS],        /* array of coefficients */
    double x,                       /* dimensionless independent variable in the range [-1, +1] */
    double result[CHEB_MAX_DIM]);   /* output approximation vector */

/*
    A utility function to scale a dimensional variable 't' in the problem domain
    (for example, a time offset expressed in seconds)
    to a dimensionless variable 'x' in the range [-1, +1] as required by ChebApprox().
*/
double ChebScale(
    double t_min,   /* the lower bound of the domain-specific range */
    double t_max,   /* the upper bound of the domain-specific range */
    double t);      /* a dimensional variable in the range [t_min, t_max]. */

#endif /* __DDC_CHEBYSHEV_H */
