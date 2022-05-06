#!/usr/bin/env python3
#
#   checksum.py  -  Don Cross <cosinekitty@gmail.com>
#
#   This script is a substitute for sha256sum/md5sum
#   on Mac OS and Windows, where it is not available by default.
#   I just needed a common way to verify correct downloads
#   for Astronomy Engine across Linux, Mac OS, and Windows.
#
import sys
import hashlib
import re

def Checksum(algorithm, filename):
    if algorithm == 'sha256':
        factory = hashlib.sha256
    elif algorithm == 'md5':
        factory = hashlib.md5
    else:
        print('ERROR: invalid algorithm: {}'.format(algorithm))
        return 1
    rc = 0
    with open(filename, 'rt') as infile:
        lnum = 0
        for line in infile:
            lnum += 1
            # 039fdcdcfc31968c6938863ac1d293854ba810bbfa0bcd72b1f4cc2d544f3d08  hygdata_v3.csv
            m = re.match(r'^([a-f0-9]+)\s+(.*)$', line)
            if not m:
                print('ERROR({} line {}): Invalid format'.format(filename, lnum))
                return 1
            correctHash = m.group(1)
            verifyFileName = m.group(2)
            BLOCKSIZE = 0x10000
            h = factory()
            # https://www.quickprogrammingtips.com/python/how-to-calculate-sha256-hash-of-a-file-in-python.html
            with open(verifyFileName, 'rb') as binfile:
               for block in iter(lambda: binfile.read(BLOCKSIZE), b""):
                    h.update(block)
            if correctHash == h.hexdigest():
                print('OK:  {}'.format(verifyFileName))
            else:
                print('BAD: {}'.format(verifyFileName))
                rc = 1
    return rc

if __name__ == '__main__':
    if len(sys.argv) == 3:
        algorithm = sys.argv[1]
        filename = sys.argv[2]
        sys.exit(Checksum(algorithm, filename))

    print('USAGE: checksum.py algorithm checkfile')
    sys.exit(1)
