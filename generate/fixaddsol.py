#!/usr/bin/env python3
import sys
import re

def Translate(line):
    m = re.match(r'^    AddSol\(\s*([0-9\+\-\.]+)\s*,\s*([0-9\+\-\.]+)\s*,\s*([0-9\+\-\.]+)\s*,\s*([0-9\+\-\.]+)\s*,\s*([\+\-]?\d+)\s*,\s*([\+\-]?\d+)\s*,\s*([\+\-]?\d+)\s*,\s*([\+\-]?\d+)\s*\)\s*$', line)
    if m:
        cl = float(m.group(1))
        cs = float(m.group(2))
        cg = float(m.group(3))
        cp = float(m.group(4))
        p = int(m.group(5))
        q = int(m.group(6))
        r = int(m.group(7))
        s = int(m.group(8))

        text = '\n'
        text += '    # AddSol({}, {}, {}, {}, {}, {}, {}, {})\n'.format(cl, cs, cg, cp, p, q, r, s)
        text += '    (a, b) = (1, 0)\n'
        if p != 0:
            text += '    (a, b) = AddThe(a, b, co[{0:d}][1], si[{0:d}][1])\n'.format(p)
        if q != 0:
            text += '    (a, b) = AddThe(a, b, co[{0:d}][2], si[{0:d}][2])\n'.format(q)
        if r != 0:
            text += '    (a, b) = AddThe(a, b, co[{0:d}][3], si[{0:d}][3])\n'.format(r)
        if s != 0:
            text += '    (a, b) = AddThe(a, b, co[{0:d}][4], si[{0:d}][4])\n'.format(s)

        text += '    DLAM  += {} * b\n'.format(cl)
        text += '    DS    += {} * b\n'.format(cs)
        text += '    GAM1C += {} * a\n'.format(cg)
        text += '    SINPI += {} * a\n'.format(cp)
        
        return text
    return line

def FixAddSol(inFileName, outFileName):
    with open(inFileName, 'rt') as infile:
        with open(outFileName, 'wt') as outfile:
            lnum = 0
            for line in infile:
                lnum += 1
                xlat = Translate(line)
                outfile.write(xlat)

if __name__ == '__main__':
    FixAddSol('template/old.py', 'template/astronomy.py')
    sys.exit(0)
