#!/usr/bin/env python3
import sys
import re
import math


class Sample:
    def __init__(self, tt, pos, vel):
        self.tt = tt
        self.pos = pos
        self.vel = vel


def LoadRef(filename):
    data = [[], [], [], []]
    with open(filename, 'rt') as infile:
        lnum = 0
        for line in infile:
            lnum += 1
            token = line.split()
            if lnum == 1:
                if token != ['tt', 'm', 'px', 'py', 'pz', 'vx', 'vy', 'vz']:
                    print('ERROR: invalid reference file:', filename)
                    sys.exit(1)
            else:
                tt = float(token[0])
                moon = int(token[1])
                pos = [float(t) for t in token[2:5]]
                vel = [float(t) for t in token[5:8]]
                data[moon].append(Sample(tt, pos, vel))
    return data


def ParseVector(line):
    m = re.match(r'[VX=\s]+([0-9\.E\-\+]+)[VY=\s]+([0-9\.E\-\+]+)[VZ=\s]+([0-9\.E\-\+]+)', line)
    if not m:
        print('FAILED TO PARSE VECTOR FROM LINE:')
        print(line)
        sys.exit(1)
    return [float(m.group(i+1)) for i in range(3)]


def LoadJpl(filename):
    data = []
    with open(filename, 'rt') as infile:
        lnum = 0
        found = False
        for line in infile:
            line = line.strip()
            lnum += 1
            if not found:
                if line == '$$SOE':
                    found = True
                    part = 0
                elif line.startswith('Revised:'):
                    moon = int(line.split()[-1]) - 501
            elif line == '$$EOE':
                break
            else:
                if part == 0:
                    # 2446545.000000000 = A.D. 1986-Apr-24 12:00:00.0000 TDB
                    tt = float(line.split()[0]) - 2451545.0
                elif part == 1:
                    #  X = 1.134408131605554E-03 Y =-2.590904586750408E-03 Z =-7.490427225904720E-05
                    pos = ParseVector(line)
                elif part == 2:
                    #  VX= 9.148038778472862E-03 VY= 3.973823407182510E-03 VZ= 2.765660368640458E-04
                    vel = ParseVector(line)
                    data.append(Sample(tt, pos, vel))
                part = (part + 1) % 3
    return moon, data


def Compare(ref, jpl):
    if len(ref) != len(jpl):
        print('Lists are different lengths.')
        sys.exit(1)

    sum = 0.0
    n = 0
    for (a, b) in zip(ref, jpl):
        if a.tt != b.tt:
            print('Time mismatch: {} != {}'.format(a.tt, b.tt))

        diff = math.sqrt((a.pos[0] - b.pos[0])**2 + (a.pos[1] - b.pos[1])**2 + (a.pos[2] - b.pos[2])**2)
        mag = math.sqrt(a.pos[0]**2 + a.pos[1]**2 + a.pos[2]**2)
        n += 1
        sum += diff / mag

    print('n={}, fit={:0.3e}'.format(n, sum / n))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('USAGE: jmcompare.py ref.txt jpl.txt')
        sys.exit(1)
    ref = LoadRef(sys.argv[1])
    moon, jpl = LoadJpl(sys.argv[2])
    Compare(ref[moon], jpl)
    sys.exit(0)
