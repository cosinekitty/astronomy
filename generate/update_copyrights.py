#!/usr/bin/env python
# Copyright (c) 2019-2021 Don Cross <cosinekitty@gmail.com>
import sys
import re
import datetime

CurrentYear = str(datetime.datetime.utcnow().year)

def UpdateCopyrights(fn):
    count = 0
    rx = re.compile(r'Copyright\s+\(c\)\s+2019(-(\d+))?\s+Don\s+Cross')
    update = 'Copyright (c) 2019-' + CurrentYear + ' Don Cross'
    newlines = []
    with open(fn, 'rt') as infile:
        for line in infile:
            m = rx.search(line)
            if m:
                if m.group(2) != CurrentYear:
                    count += 1
                    line = rx.sub(update, line)
            newlines.append(line)
    if count > 0:
        print('update_copyrights.py: {}'.format(fn))
        with open(fn, 'wt') as outfile:
            for line in newlines:
                outfile.write(line)
        return 1
    return 0

def main(fnlist):
    nfiles = 0
    for fn in fnlist:
        nfiles += UpdateCopyrights(fn)
    print('update_copyrights: updated {:d} files'.format(nfiles))
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
