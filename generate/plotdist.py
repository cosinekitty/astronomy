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
            if lnum > 1:
                token = line.split(',')
                if len(token) < 2:
                    raise Exception('Invalid CSV line {} in file {}'.format(lnum, inFileName))
                xlist.append(float(token[0]))
                ylist.append(float(token[1]))
    return xlist, ylist


def PlotData(novasFileName, topFileName):
    novas_tlist, novas_rlist = LoadCsv(novasFileName)
    top_tlist, top_rlist = LoadCsv(topFileName)
    plt.plot(novas_tlist, novas_rlist, 'b.')
    plt.plot(top_tlist, top_rlist, 'r.')
    plt.show()
    plt.close('all')
    return 0


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('USAGE: plotdist.py novas_dist.csv top_dist.csv')
        sys.exit(1)
    sys.exit(PlotData(sys.argv[1], sys.argv[2]))
