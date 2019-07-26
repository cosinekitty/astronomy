#!/usr/bin/env python3
#
#   astro_demo_common.py  -  Don Cross  -  2019-07-26
#
#   Utility functions shared by Python demo programs.
#
import sys
import re
import astronomy

def ParseTime(text):
    m = re.match(r'^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2}\.?\d*)Z$', text)
    if not m:
        raise Exception('Invalid date/time format in string: "{0}"'.format(text))
    year = int(m.group(1))
    month = int(m.group(2))
    day = int(m.group(3))
    hour = int(m.group(4))
    minute = int(m.group(5))
    second = float(m.group(6))
    return astronomy.Time.Make(year, month, day, hour, minute, second)

def ParseArgs(args):
    # [demo].py latitude longitude [utc]
    if len(args) not in [3, 4]:
        print('USAGE: {} latitude longitude [yyyy-mm-ddThh:mm:ssZ]'.format(args[0]))
        sys.exit(1)
    latitude = float(args[1])
    longitude = float(args[2])
    time = ParseTime(args[3])
    observer = astronomy.Observer(latitude, longitude)
    return (observer, time)
