#!/usr/bin/env python3
#
#   Obtains a list of public functions in each langauge
#   version of Astronomy Engine. Compares them and verifies
#   that the same functions are available among all languages.
#
import sys
import re

def FuncList_Python(filename):
    funcset = set()
    with open(filename, 'rt') as infile:
        for line in infile:
            m = re.match(r'^def\s+([A-Z][A-Za-z_]+)\s*\(', line)
            if m:
                funcset.add(m.group(1))
    return funcset

def FuncList_C(filename):
    funcset = set()
    with open(filename, 'rt') as infile:
        for line in infile:
            m = re.match(r'^([a-z_]+).*?Astronomy_([A-Za-z_]+)\s*\(', line)
            if m:
                if m.group(1) != 'static':
                    funcset.add(m.group(2))
    return funcset

def FuncList_Csharp(filename):
    funcset = set()
    inAstronomyClass = False
    with open(filename, 'rt') as infile:
        for line in infile:
            if not inAstronomyClass:
                if line == '    public static class Astronomy\n':
                    inAstronomyClass = True
            else:
                if line == '    }\n':
                    break
                m = re.match(r'^        public static .*?([A-Za-z_]+)\(', line)
                if m:
                    funcset.add(m.group(1))
    return funcset

def Funclist_JavaScript(filename):
    funcset = set()
    with open(filename, 'rt') as infile:
        for line in infile:
            m = re.match(r'^Astronomy\.([A-Za-z_]+)\s*=\s*function\s*\(', line)
            if m:
                funcset.add(m.group(1))
    return funcset

if __name__ == '__main__':
    c = FuncList_C('../source/c/astronomy.c')
    cs = FuncList_Csharp('../source/csharp/astronomy.cs')
    js = Funclist_JavaScript('../source/js/astronomy.js')
    py = FuncList_Python('../source/python/astronomy.py')
    all = sorted(c | cs | js | py)
    for f in all:
        line = '{0:<25s}'.format(f)
        line += '   ' + ('c ' if f in c else '  ')
        line += '   ' + ('cs' if f in cs else '  ')
        line += '   ' + ('js' if f in js else '  ')
        line += '   ' + ('py' if f in py else '  ')
        print(line)

    sys.exit(0)
