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
