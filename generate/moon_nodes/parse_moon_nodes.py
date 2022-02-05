#!/usr/bin/env python3
import sys
import re

class Node:
    def __init__(self, kind, year, month, day, hour, minute, ra, dec):
        self.kind = kind
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.minute = minute
        self.ra = ra
        self.dec = dec

    def __str__(self):
        return '{} {:4d}-{:02d}-{:02d}T{:02d}:{:02d}Z {:9.4f} {:9.4f}'.format(
            self.kind, self.year, self.month, self.day, self.hour, self.minute, self.ra, self.dec
        )

MonthTable = {
    'Jan':  1,
    'Feb':  2,
    'Mar':  3,
    'Apr':  4,
    'May':  5,
    'Jun':  6,
    'Jul':  7,
    'Aug':  8,
    'Sep':  9,
    'Oct': 10,
    'Nov': 11,
    'Dec': 12
}

def HourMin(htext, mtext):
    if htext[0] in '+-':
        x = float(htext[1:]) + float(mtext)/60.0
        if htext[0] == '+':
            return x
        return -x
    return float(htext) + float(mtext)/60.0

def ParseNode(kind, year, text):
    #            1         2         3
    #  01234567890123456789012345678901234
    # "Jun 10  20:00 A  05h07.2m  +22:53.1"
    if text.strip() == '':
        return None
    month = MonthTable[text[0:3]]
    day = int(text[4:6], 10)
    hour = int(text[8:10], 10)
    minute = int(text[11:13], 10)
    ra = HourMin(text[17:19], text[20:24])
    dec = HourMin(text[27:30], text[31:])
    return Node(kind, year, month, day, hour, minute, ra, dec)

if __name__ == '__main__':
    inFileName = 'espenak_nodes.txt'
    outFileName = 'moon_nodes.txt'
    with open(outFileName, 'wt') as outfile:
        with open(inFileName, 'rt') as infile:
            year = 0
            lnum = 0
            for line in infile:
                lnum += 1
                line = line.rstrip()
                # 2001    Jan 09  13:53 t  07h07.4m  +22:32.1       Jan 22  22:22    19h07.5m  -22:31.9
                if len(line) in [44, 86]:
                    if re.match(r'[0-9]{4}', line[1:5]):
                        year = int(line[1:5])
                    elif line[1:5] != '    ':
                        print('Syntax error in file {} line {}'.format(inFileName, lnum))
                        sys.exit(1)
                    asc = ParseNode('A', year, line[9:44])
                    if asc:
                        outfile.write(str(asc))
                        outfile.write('\n')
                    dsc = ParseNode('D', year, line[51:86])
                    if dsc:
                        outfile.write(str(dsc))
                        outfile.write('\n')
    sys.exit(0)
