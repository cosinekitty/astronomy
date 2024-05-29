#!/usr/bin/env python3
import sys
import math
import dontrig


def Pass(message:str) -> int:
    print('trig.py PASS:', message)
    return


def SineTest() -> int:
    # See how long it takes to converge for a problematic value near 1.
    angle = -0.9535784309745123
    correct = math.sin(angle)
    sum = dontrig.xsin(angle)
    diff = sum - correct
    print('SineTest: correct={:0.16f}, sum={:0.16f}, diff={:g}'.format(correct, sum, diff))
    tolerance = 1.2e-16
    if abs(diff) > tolerance:
        print('SineTest FAIL: EXCESSIVE ERROR')
        return 1
    return Pass('SineTest')


def CosineTest() -> int:
    # See how long it takes to converge for a problematic value near 1.
    angle = -0.9535784309745123
    correct = math.cos(angle)
    sum = dontrig.xcos(angle)
    diff = sum - correct
    print('CosineTest: correct={:0.16f}, sum={:0.16f}, diff={:g}'.format(correct, sum, diff))
    tolerance = 1.2e-16
    if abs(diff) > tolerance:
        print('CosineTest FAIL: EXCESSIVE ERROR')
        return 1
    return Pass('CosineTest')


def UnitTest() -> int:
    return (
        SineTest() or
        CosineTest() or
        Pass('TrigTest')
    )


def GenerateOutputTable(filename: str) -> int:
    # Now that we have verified xsin and xcos are consistent
    # with math.sin and math.cos, let's generate a table of correct values
    # using math.sin and math.cos, to automate downstream unit testing.
    deg = -720.0
    incr = 0.1
    maxDeg = +720.01
    with open(filename, 'wt') as outfile:
        while deg <= maxDeg:
            rad = math.radians(deg)
            c = math.cos(rad)
            s = math.sin(rad)
            outfile.write('{:6.1f} {:22.16f} {:22.16f}\n'.format(deg, c, s))
            deg += incr
    return Pass('GenerateOutputTable')


if __name__ == '__main__':
    dontrig.Debug = True
    sys.exit(
        UnitTest() or
        GenerateOutputTable('trigonometry/trig.txt')
    )
