#!/usr/bin/env python3
import sys
import re

#--------------------------------------------------------------------------------------

def ParseDuration(s):
    if s == '-':
        return 0
    m = re.match(r'(\d+)m', s)
    if m:
        return int(m.group(1))
    raise Exception('Invalid duration string "{}"'.format(s))

def ParseMonth(s):
    return 1 + ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'].index(s)

def ParseGeo(text, suffixes):
    polarity = +1
    for s in suffixes:
        if text.endswith(s):
            return polarity * float(text[:-1])
        polarity *= -1
    raise Exception('Invalid geographic coordinate "{}"'.format(text))

#--------------------------------------------------------------------------------------

def FixLunarEclipseData():
    # '  2022 May 16  04:11   T+ 131  -0.253  2.397  1.419 104m  43m  15.6  15.52 -19.3'
    r = re.compile(r'''^\s{2}
        (\d{4})\s+                                                  # [1] year
        (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+        # [2] month name
        (\d{2})\s+                                                  # [3] day
        (\d{2}):(\d{2})\s+                                          # [4] hour, [5] minute
        \S+\s+                                                      #     (ignore type)
        \d+\s+                                                      #     (ignore Saros number)
        \S+\s+                                                      #     (ignore gamma)
        \S+\s+                                                      #     (ignore pen mag)
        \S+\s+                                                      #     (ignore umb mag)
        (\d+m|-)\s+                                                 # [6] semi-duration of partial, minutes
        (\d+m|-)\s+                                                 # [7] semi-duration of total, minutes
    ''',
        re.VERBOSE)

    with open('lunar_eclipse.txt', 'wt') as outfile:
        for fn in ['le1701.html', 'le1801.html', 'le1901.html', 'le2001.html', 'le2101.html']:
            with open(fn, 'rt') as infile:
                for line in infile:
                    m = r.match(line)
                    if m:
                        year = int(m.group(1))
                        month = ParseMonth(m.group(2))
                        day = int(m.group(3))
                        hour = int(m.group(4))
                        minute = int(m.group(5))
                        partial = ParseDuration(m.group(6))
                        total = ParseDuration(m.group(7))
                        outfile.write('{:04d}-{:02d}-{:02d}T{:02d}:{:02d}Z {:3d} {:3d}\n'.format(year, month, day, hour, minute, partial, total))

    return 0

#--------------------------------------------------------------------------------------

def FixSolarEclipseData():
    # <a href="/web/20080228210550/http://sunearth.gsfc.nasa.gov/eclipse/5MCSEmap/1701-1800/1719-02-19.gif">08835</a>
    # 1719 Feb 19  06:52:57
    # 10  -3474
    # <a href="/web/20080228210550/http://sunearth.gsfc.nasa.gov/eclipse/SEsaros/SEsaros116.html">116</a>
    # A    0.6856  0.9250  30.5N  68.6E  47 163  384  09m01s
    r = re.compile(r'''^
        <a\s+href="[^"]+">(\d+)</a>\s+                              # [ 1] eclipse number
        (\d{4})\s+                                                  # [ 2] year
        (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+        # [ 3] month name
        (\d{2})\s+                                                  # [ 4] day
        (\d{2}):                                                    # [ 5] hour
        (\d{2}):                                                    # [ 6] minute
        (\d{2})\s+                                                  # [ 7] second
        (-?\d+)\s+                                                  # [ 8] Delta T
        \S+\s+                                                      #      (ignore lunation number)
        <a\s+href="[^"]+">\d+</a>\s+                                #      (ignore Saros number)
        ([PATH])\S?\s+                                              # [ 9] eclipse type
        \S+\s+                                                      #      (ignore gamma number)
        \S+\s+                                                      #      (ignore eclipse mag)
        (\d+\.\d[NS])\s+                                            # [10] latitude of greatest eclipse
        (\d+\.\d[EW])\s+                                            # [11] longitude of greatest eclipse
    ''',
        re.VERBOSE)

    prev_eclnum = None
    with open('solar_eclipse.txt', 'wt') as outfile:
        for fn in['se1701.html', 'se1801.html', 'se1901.html', 'se2001.html', 'se2101.html']:
            with open(fn, 'rt') as infile:
                lnum = 0
                for line in infile:
                    lnum += 1
                    m = r.match(line)
                    if m:
                        #print(line, end='')
                        eclnum = int(m.group(1))
                        if (prev_eclnum is not None) and (prev_eclnum + 1 != eclnum):
                            print('norm.py FATAL(FixSolarEclipseData): Unexpected eclipse number {} (prev = {}) at {} line {}'.format(eclnum, prev_eclnum, fn, lnum))
                            return 1
                        prev_eclnum = eclnum
                        year = int(m.group(2))
                        month = ParseMonth(m.group(3))
                        day = int(m.group(4))
                        hour = int(m.group(5))
                        minute = int(m.group(6))
                        second = int(m.group(7))
                        delta_t = int(m.group(8))
                        ecltype = m.group(9)
                        lat = ParseGeo(m.group(10), 'NS')
                        lon = ParseGeo(m.group(11), 'EW')
                        outfile.write('{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}Z {:4d} {:1s} {:7.1f} {:7.1f}\n'.format(year, month, day, hour, minute, second, delta_t, ecltype, lat, lon))
    return 0

#--------------------------------------------------------------------------------------

if __name__ == '__main__':
    sys.exit(FixLunarEclipseData() or FixSolarEclipseData())
