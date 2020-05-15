#!/usr/bin/env python3

import sys
import re
import matplotlib.pyplot as plt


def LoadCsv(inFileName):
    xlist = []
    ylist = []
    with open(inFileName, 'rt') as infile:
        lnum = 0
        for line in infile:
            line = line.strip()
            lnum += 1
            if lnum == 1:
                if line != '"year","delta_t"':
                    raise Exception('Invalid CSV file: ' + inFileName)
            else:
                token = line.split(',')
                if len(token) != 2:
                    raise Exception('Invalid CSV line {} in file {}'.format(lnum, inFileName))
                xlist.append(float(token[0]))
                ylist.append(float(token[1]))
    return xlist, ylist


def LoadJplData(inFileName):
    xlist = []
    ylist = []
    with open(inFileName, 'rt') as infile:
        prevYear = None
        lnum = 0
        indata = False
        for line in infile:
            line = line.strip()
            lnum += 1
            if not indata:
                if line == '$$SOE':
                    indata = True
            else:
                if line == '$$EOE':
                    break
                # 2082-Nov-12 00:00     123.45882  21.10519 124.67208  20.84737  29.5918 -34.2445    69.182656
                # 2083-Jan-01 00:00 Am   61.88124  25.69762  63.15127  25.91410  83.9978  51.6684    69.183873
                token = line.split()
                if len(token) < 9:
                    raise Exception('Invalid tokens on line {} of file {}'.format(lnum, inFileName))
                m = re.match(r'^(\d{4})-([A-Z][a-z][a-z])-\d{2}$', token[0])
                if not m:
                    raise Exception('Invalid date format on line {} of file {}'.format(lnum, inFileName))
                year = int(m.group(1))
                month = m.group(2)
                dt = float(token[-1])
                if year != prevYear:
                    prevYear = year
                    if month != 'Jan':
                        raise Exception('Missed January for year {}'.format(year))
                    xlist.append(year)
                    ylist.append(dt)
    return xlist, ylist


def PlotData(extrapFileName, actualFileName, jplFileName):
    xlist, ylist = LoadCsv(extrapFileName)
    plt.plot(xlist, ylist, 'b-')

    xlist, ylist = LoadCsv(actualFileName)
    plt.plot(xlist, ylist, 'r-')

    xlist, ylist = LoadJplData(jplFileName)
    plt.plot(xlist, ylist, 'g-')

    plt.show()
    plt.close('all')
    return 0


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('USAGE: plot_delta_t.py extrap.csv actual.csv jpldata.txt')
        sys.exit(1)
    sys.exit(PlotData(sys.argv[1], sys.argv[2], sys.argv[3]))
