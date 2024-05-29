#!/usr/bin/env python3
import sys
import math

Debug = False

def xsin(angle:float) -> float:
    angleSquared = angle * angle
    sum = 0.0
    fact = 1
    pow = angle
    n = 1
    while n <= 99:
        # n = 1, 3, 5, ...
        term = pow / fact
        prev = sum
        sum += term
        if Debug:
            print('xsin: n={:02d}  term={:24.16e}  sum={:20.16f}   fact={:d}'.format(n, term, sum, fact))
        if prev == sum:
            return sum
        pow *= angleSquared
        fact *= -((n+1) * (n+2))
        n += 2
    raise Exception('xsin: failure to converge')


def xcos(angle:float) -> float:
    angleSquared = angle * angle
    sum = 0.0
    fact = 1
    pow = 1.0
    n = 0
    while n <= 99:
        # n = 0, 2, 4, ...
        term = pow / fact
        prev = sum
        sum += term
        if Debug:
            print('xcos: n={:02d}  term={:24.16e}  sum={:20.16f}   fact={:d}'.format(n, term, sum, fact))
        if prev == sum:
            return sum
        pow *= angleSquared
        fact *= -((n+1) * (n+2))
        n += 2
    raise Exception('xcos: failure to converge')


def _Pass(message:str) -> int:
    print('trig.py PASS:', message)
    return


def _SineTest() -> int:
    # See how long it takes to converge for a problematic value near 1.
    angle = -0.9535784309745123
    correct = math.sin(angle)
    sum = xsin(angle)
    diff = sum - correct
    print('SineTest: correct={:0.16f}, sum={:0.16f}, diff={:g}'.format(correct, sum, diff))
    tolerance = 1.2e-16
    if abs(diff) > tolerance:
        print('SineTest FAIL: EXCESSIVE ERROR')
        return 1
    return _Pass('SineTest')


def _CosineTest() -> int:
    # See how long it takes to converge for a problematic value near 1.
    angle = -0.9535784309745123
    correct = math.cos(angle)
    sum = xcos(angle)
    diff = sum - correct
    print('CosineTest: correct={:0.16f}, sum={:0.16f}, diff={:g}'.format(correct, sum, diff))
    tolerance = 1.2e-16
    if abs(diff) > tolerance:
        print('CosineTest FAIL: EXCESSIVE ERROR')
        return 1
    return _Pass('CosineTest')


def _UnitTest() -> int:
    return (
        _SineTest() or
        _CosineTest() or
        _Pass('TrigTest')
    )


if __name__ == '__main__':
    Debug = True
    sys.exit(_UnitTest())
