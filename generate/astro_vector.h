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

#ifndef __DDC_ASTRO_VECTOR_H
#define __DDC_ASTRO_VECTOR_H

typedef struct
{
    double t;
    double v[3];
}
VectorType;

typedef struct
{
    int size;
    int length;
    VectorType *array;
}
VectorListType;

typedef struct
{
    double a;
    double b;
    double c;
}
ParabolaType;

VectorType Subtract(VectorType a, VectorType b);
VectorType Add(VectorType a, VectorType b);
double Length(VectorType a);
double Dot(VectorType a, VectorType b);
double AngleBetween(VectorType a, VectorType b);
double SafeAcos(double x);
VectorType UnitVector(VectorType a);

void InitVectorList(VectorListType *list);
int VectorListAppend(VectorListType *list, double t, const double v[3]);
void ClearVectorList(VectorListType *list);
void FreeVectorList(VectorListType *list);

int QuadraticRegression(const VectorListType *list, ParabolaType parab[3]);
void QuadraticEval(const ParabolaType parab[3], double t, VectorType *approx);

int TestRegression(const VectorListType *list, const ParabolaType parab[3], char tag);

double Det(const double a[3][3]);     /* determinant of a 3x3 matrix */

#endif /* __DDC_ASTRO_VECTOR_H */
