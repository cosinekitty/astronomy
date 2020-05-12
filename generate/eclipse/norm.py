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
                        month = 1 + ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'].index(m.group(2))
                        day = int(m.group(3))
                        hour = int(m.group(4))
                        minute = int(m.group(5))
                        partial = ParseDuration(m.group(6))
                        total = ParseDuration(m.group(7))
                        outfile.write('{:04d}-{:02d}-{:02d}T{:02d}:{:02d}Z {:3d} {:3d}\n'.format(year, month, day, hour, minute, partial, total))

    return 0

#--------------------------------------------------------------------------------------

if __name__ == '__main__':
    sys.exit(FixLunarEclipseData())
