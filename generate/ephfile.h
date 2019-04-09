/*
    ephfile.h  -  Don Cross  -  2019-03-04.
*/
#ifndef __DDC_EPH_READER
#define __DDC_EPH_READER

#include <stdio.h>
#include "chebyshev.h"

typedef struct
{
    FILE *infile;
    int lnum;
    int body;
}
eph_file_reader_t;

typedef struct
{
    int error;
    double jdStart;
    double jdDelta;
    int numpoly;
    double coeff[3][CHEB_MAX_POLYS];
}
eph_record_t;

int EphFileOpen(eph_file_reader_t *reader, const char *filename);
int EphReadRecord(eph_file_reader_t *reader, eph_record_t *record);
void EphFileClose(eph_file_reader_t *reader);

#endif /* __DDC_EPH_READER */
