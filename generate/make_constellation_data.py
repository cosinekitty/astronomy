#!/usr/bin/env python3
#
#   Use star database from https://github.com/astronexus/HYG-Database
#   to generate constellation test data of the form:
#   ra dec symbol
#

import sys
import re

def Unquote(s):
    return s.replace('"', '')

def Translate(inFileName, outFileName):

    # Ignore disagreements with HYG v3.6.1.
    # See: https://github.com/astronexus/HYG-Database/issues/21
    ignore = set([78313, 117789])

    count = 0
    with open(inFileName, 'rt') as infile:
        with open(outFileName, 'wt') as outfile:
            lnum = 0
            for line in infile:
                lnum += 1
                token = line.strip().split(',')
                if lnum == 1:
                    columns = dict((Unquote(key), col) for (col, key) in enumerate(token))
                else:
                    id = int(token[columns['id']])
                    ra = float(token[columns['ra']])
                    dec = float(token[columns['dec']])
                    sym = Unquote(token[columns['con']])
                    if sym != '':
                        if not re.match(r'^[A-Z][a-zA-Z]{2}$', sym):
                            print('make_constellation_data: bad symbol "{}" in {} line {}'.format(sym, inFileName, lnum))
                            return 1

                        if id not in ignore:
                            count += 1
                            outfile.write('{:6d} {:10.6f} {:10.6f} {:3s}\n'.format(id, ra, dec, sym))
    print('make_constellation_data: Wrote {:d} test cases.'.format(count))
    return 0

if __name__ == '__main__':
    sys.exit(Translate('hyg_v36_1.csv', 'constellation/test_input.txt'))
