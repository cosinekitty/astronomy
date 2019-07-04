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

        op = ''
        text += '    z = '
        if p != 0:
            text += op + 'ex[{}][1]'.format(p)
            op = ' * '
        if q != 0:
            text += op + 'ex[{}][2]'.format(q)
            op = ' * '
        if r != 0:
            text += op + 'ex[{}][3]'.format(r)
            op = ' * '
        if s != 0:
            text += op + 'ex[{}][4]'.format(s)
        text += '\n'

        if cl != 0:
            text += '    DLAM  += {} * z.imag\n'.format(cl)
        if cs != 0:
            text += '    DS    += {} * z.imag\n'.format(cs)
        if cg != 0:
            text += '    GAM1C += {} * z.real\n'.format(cg)
        if cp != 0:
            text += '    SINPI += {} * z.real\n'.format(cp)
        
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
