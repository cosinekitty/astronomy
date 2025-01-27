/*
    MIT License

    Copyright (c) 2019-2025 Don Cross <cosinekitty@gmail.com>

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
#include <stdlib.h>
#include <ctype.h>
#include "ephfile.h"

int EphFileOpen(eph_file_reader_t *reader, const char *filename)
{
    char line[128];
    char *key;
    char *value;

    memset(reader, 0, sizeof(eph_file_reader_t));
    reader->infile = fopen(filename, "rt");
    if (reader->infile == NULL)
        return 1;   /* cannot open file */

    while (fgets(line, sizeof(line), reader->infile))
    {
        ++(reader->lnum);

        /* parse "key=value" from line */
        for (key=line; *key != '\0' && isspace(*key); ++key);
        if (*key == '\0')
            return 0;   /* found end of header */

        /* search for '=' separating key and value */
        for (value=key+1; *value != '\0' && *value != '='; ++value);
        if (*value == '\0')
            break;

        *value++ = '\0';      /* terminate the key string by stomping on the '=' */
        if (!strcmp(key, "body"))
            reader->body = atoi(value);
    }

    fclose(reader->infile);
    reader->infile = NULL;
    return 2;   /* syntax error */
}

int EphReadRecord(eph_file_reader_t *reader, eph_record_t *record)
{
    int nscanned;
    int k;
    char line[128];

    memset(record, 0, sizeof(eph_record_t));

    if (!fgets(line, sizeof(line), reader->infile))
        return 0;       /* normal end of file */

    ++(reader->lnum);
    nscanned = sscanf(line, "%lf %lf %d", &record->jdStart, &record->jdDelta, &record->numpoly);
    if (nscanned != 3)
    {
        record->error = 1;      /* syntax error */
        return 0;
    }

    if (record->numpoly < 1 || record->numpoly > CHEB_MAX_POLYS)
    {
        record->error = 2;      /* invalid number of polynomials */
        return 0;
    }

    for (k=0; k < record->numpoly; ++k)
    {
        if (!fgets(line, sizeof(line), reader->infile))
        {
            record->error = 3;      /* unexpected end of data in the middle of a record */
            return 0;
        }
        ++(reader->lnum);
        nscanned = sscanf(line, "%lf %lf %lf", &record->coeff[0][k], &record->coeff[1][k], &record->coeff[2][k]);
        if (nscanned != 3)
        {
            record->error = 4;      /* unexpected failure to read coefficients */
            return 0;
        }
    }

    return 1;   /* success */
}

void EphFileClose(eph_file_reader_t *reader)
{
    if (reader->infile)
    {
        fclose(reader->infile);
        reader->infile = NULL;
    }
}
