#!/usr/bin/env python3

import sys
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


def PlotData(extrapFileName, actualFileName):
    xlist, ylist = LoadCsv(extrapFileName)
    plt.plot(xlist, ylist, 'b-')
    xlist, ylist = LoadCsv(actualFileName)
    plt.plot(xlist, ylist, 'r.')
    plt.show()
    plt.close('all')
    return 0


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('USAGE: plot_delta_t.py extrap.csv actual.csv')
        sys.exit(1)
    sys.exit(PlotData(sys.argv[1], sys.argv[2]))
