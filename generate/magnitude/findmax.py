#!/usr/bin/env python3
import sys
import re
from datetime import datetime

class Item:
    def __init__(self, mag, date, elong, rlon):
        self.mag = mag
        self.date = date
        self.elong = elong
        self.rlon = rlon

def DateStr(date):
    return date.isoformat()[:-3] + 'Z'

def FindMaxMagEvents(filename):
    with open(filename, 'rt') as infile:
        lnum = 0
        prev_item = None
        first_brightest_item = None
        seen_dimming = False
        min_rlon = None
        max_rlon = None
        for line in infile:
            lnum += 1
            line = line.strip()
            # [ Date__(UT)__HR:MN      APmag  S-brt    S-O-T /r    S-T-O    O-P-T]
            # [ 2001-Apr-28 08:00      -4.73   1.43  35.4650 /L 126.2646  18.2704]
            m = re.match(r'^(\d{4}-[A-Z][a-z]{2}-\d{2}\s+\d{2}:\d{2})\s+(-?\d+\.\d{2})\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+/([LT])\s+(\d+\.\d+)\s+(\d+\.\d+)$', line)
            if m:
                date = datetime.strptime(m.group(1), '%Y-%b-%d %H:%M')
                mag = float(m.group(2))
                elong = float(m.group(4))
                rlon = float(m.group(7))
                item = Item(mag, date, elong, rlon)
                if prev_item is not None:
                    if item.mag < prev_item.mag:
                        # Getting brighter
                        first_brightest_item = item
                    elif item.mag > prev_item.mag:
                        # Getting dimmer ... did we just find a peak?
                        # Ignore the minor peaks that happen near inferior conjunction... they are too hard to see anyway.
                        if first_brightest_item and seen_dimming and (first_brightest_item.elong > 10):
                            print('{} {} {:5.2f} {:5.2f} {:4.2f}'.format(DateStr(first_brightest_item.date), DateStr(prev_item.date), first_brightest_item.rlon, prev_item.rlon, first_brightest_item.mag))
                            if min_rlon is None:
                                min_rlon = min(first_brightest_item.rlon, prev_item.rlon)
                                max_rlon = max(first_brightest_item.rlon, prev_item.rlon)
                            else:
                                min_rlon = min(min_rlon, first_brightest_item.rlon, prev_item.rlon)
                                max_rlon = max(max_rlon, first_brightest_item.rlon, prev_item.rlon)
                        seen_dimming = True
                        first_brightest_item = None
                    else:
                        # Staying the same magnitude
                        if first_brightest_item and (item.mag != first_brightest_item.mag):
                            first_brightest_item = None
                prev_item = item
            else:
                # There is a break in the data... "n.a." or end of data
                seen_dimming = False
                prev_item = None
                first_brightest_item = None
        #print('min_rlon={:0.2f}, max_rlon={:0.2f}'.format(min_rlon, max_rlon))

if __name__ == '__main__':
    FindMaxMagEvents('Venus2.txt')
    sys.exit(0)
