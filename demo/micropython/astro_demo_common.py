#!/usr/bin/env python3
#
#   astro_demo_common.py  -  Don Cross  -  2019-07-26
#
#   Utility functions shared by Python demo programs.
#
import sys
import astronomy

def ParseArgs(args: list[str]) -> tuple[astronomy.Observer, astronomy.Time]:
    if len(args) not in [3, 4]:
        print('USAGE: {} latitude longitude [yyyy-mm-ddThh:mm:ssZ]'.format(args[0]))
        sys.exit(1)
    latitude = float(args[1])
    longitude = float(args[2])
    if len(args) == 4:
        time = astronomy.Time.Parse(args[3])
    else:
        time = astronomy.Time.Now()
    observer = astronomy.Observer(latitude, longitude)
    return (observer, time)
