#!/usr/bin/env python3
import sys
import math

def main() -> int:
    # See how long it takes to converge for a problematic value near 1.
    angle = -0.9535784309745123
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
        print('n={:02d}  term={:24.16e}  sum={:20.16f}   fact={:d}'.format(n, term, sum, fact))
        if prev == sum:
            break
        pow *= angleSquared
        fact *= -((n+1) * (n+2))
        n += 2
    correct = math.sin(angle)
    diff = sum - correct
    print('correct={:0.16f}, sum={:0.16f}, diff={:g}'.format(correct, sum, diff))
    return 0

if __name__ == '__main__':
    sys.exit(main())
