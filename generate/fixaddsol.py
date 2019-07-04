#!/usr/bin/env python3
import sys
import re

def AddThe(first, text, p, i):
    if p != 0:
        cosi = 'co[{0}][{1}], si[{0}][{1}]'.format(p, i)
        if first:
            first = False
            text += '    (a, b) = ({0})\n'.format(cosi)
        else:
            text += '    (a, b) = AddThe(a, b, {0})\n'.format(cosi)
    return first, text

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

        first = True
        first, text = AddThe(first, text, p, 1)
        first, text = AddThe(first, text, q, 2)
        first, text = AddThe(first, text, r, 3)
        first, text = AddThe(first, text, s, 4)

        if cl != 0:
            text += '    DLAM  += {} * b\n'.format(cl)

        if cs != 0:
            text += '    DS    += {} * b\n'.format(cs)

        if cg != 0:
            text += '    GAM1C += {} * a\n'.format(cg)

        if cp != 0:
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
