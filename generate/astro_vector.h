/*
    astro_vector.h  -  Don Cross  -  2019-03-14
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
