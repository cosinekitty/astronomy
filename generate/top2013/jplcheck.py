#!/usr/bin/env python3
import sys
import re
from math import sqrt, degrees

class Record:
    def __init__(self, planet, jd, x, y, z, vx, vy, vz):
        self.planet = planet
        self.jd = jd
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz

def JplRead(fn, planet):
    with open(fn) as jplfile:
        seek = True
        jlist = []
        for line in jplfile:
            line = line.strip()
            if seek:
                if line == '$$SOE':
                    seek = False
                    row = 0
            else:
                if line == '$$EOE':
                    break
                if row == 0:
                    # 2411545.000000000 = A.D. 1890-Jun-26 12:00:00.0000 TDB
                    jd = float(line.split()[0])
                elif row == 1:
                    #  X = 1.213238742051184E+01 Y = 2.532261522896372E+01 Z = 1.006281191830243E+01
                    m = re.match(r'^\s*X\s*=\s*(\S+)\s+Y\s*=\s*(\S+)\s+Z\s*=\s*(\S+)', line)
                    if not m:
                        raise Exception('Bad data format in [{}]'.format(line))
                    x = float(m.group(1))
                    y = float(m.group(2))
                    z = float(m.group(3))
                elif row == 2:
                    #  VX=-2.880310222354256E-03 VY= 1.177013903016748E-03 VZ= 5.534782388307750E-04
                    m = re.match(r'^\s*VX\s*=\s*(\S+)\s+VY\s*=\s*(\S+)\s+VZ\s*=\s*(\S+)', line)
                    if not m:
                        raise Exception('Bad data format in [{}]'.format(line))
                    vx = float(m.group(1))
                    vy = float(m.group(2))
                    vz = float(m.group(3))
                else:
                    jlist.append(Record(planet, jd, x, y, z, vx, vy, vz))
                row = (row + 1) % 4
        return jlist


def TopRead(fn):
    with open(fn) as topfile:
        row = 0
        tlist = []
        for line in topfile:
            line = line.strip()
            token = line.split()
            if row == 0:
                planet = int(token[0])
                jd = float(token[1])
            elif row == 3:
                x, y, z, vx, vy, vz = [float(t) for t in token]
                tlist.append(Record(planet, jd, x, y, z, vx, vy,vz))
            row = (row + 1) % 4
        return tlist


def Compare(t, j):
    if t.jd != j.jd:
        raise Exception('Mismatching times')
    if t.planet != j.planet:
        raise Exception('Mismatching planets')
    range = sqrt(j.x**2 + j.y**2 + j.z**2)
    diff  = sqrt((t.x-j.x)**2 + (t.y-j.y)**2 + (t.z-j.z)**2)
    arcmin = 60 * degrees(diff / range)
    print('planet {:1d}  jd {:10.1f}  arcmin {:10.6f}'.format(t.planet, t.jd, arcmin))


def main():
    tlist = TopRead('correct.txt')
    print('TOP2013 record count =', len(tlist))
    jlist = []
    for planet in range(5, 10):
        jlist += JplRead('jplhor_{:d}.txt'.format(planet), planet)
    print('JPLHOR record count =', len(jlist))
    for (t, j) in zip(tlist, jlist):
        Compare(t, j)
    return 0

if __name__ == '__main__':
    sys.exit(main())
